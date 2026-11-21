"""
Linear Control Systems Toolkit (Tkinter)
Tabs: Time Response | Frequency Response | Root Locus | Routh-Hurwitz | Nichols Chart | PID Tuning | Compensator
Requires: customtkinter, numpy, scipy, matplotlib
          ->  pip install customtkinter numpy scipy matplotlib
Run:      python control_gui.py
"""
import tkinter as tk
from tkinter import messagebox
import customtkinter as ctk
import numpy as np
from scipy import signal
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk


# ---------------------------------------------------------------- helpers
def parse(text):
    """'1 2, 5' -> array([1., 2., 5.]) (highest power first)."""
    try:
        v = np.array([float(x) for x in text.replace(",", " ").split()])
    except ValueError:
        raise ValueError("Coefficients must be numbers, e.g.  1 2 5")
    v = np.trim_zeros(v, "f")
    if v.size == 0:
        raise ValueError("Polynomial cannot be empty / all zeros")
    return v


def closed_loop(n, d, unity):
    return (n, np.trim_zeros(np.polyadd(d, n), "f")) if unity else (n, d)


def crossings(x, y, level):
    """Log-interpolated x where y crosses `level`."""
    s = y - level
    idx = np.flatnonzero(s[:-1] * s[1:] < 0)
    lx = np.log10(x)
    return [10 ** (lx[i] - s[i] * (lx[i + 1] - lx[i]) / (s[i + 1] - s[i])) for i in idx]


def error_constants(n, d):
    """System type and Kp, Kv, Ka for unity feedback with open-loop G = n/d."""
    N = 0
    dd = d.copy()
    while dd.size > 1 and abs(dd[-1]) < 1e-12:
        dd = dd[:-1]
        N += 1
    c = n[-1] / dd[-1] if abs(n[-1]) > 1e-12 else 0.0
    res = []
    for k in range(3):
        res.append(0.0 if N < k else (c if N == k else np.inf))
    return N, res


def fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "inf" if v is not None and np.isinf(v) else "n/a"
    return f"{v:.4g}"

SUP = str.maketrans("0123456789", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def poly_str(c, var="s"):
    """[1, 2, 0] -> 's\u00b2 + 2s'."""
    c = np.trim_zeros(np.asarray(c, float), "f")
    if c.size == 0:
        return "0"
    N, out = len(c) - 1, ""
    for i, a in enumerate(c):
        p = N - i
        if abs(a) < 1e-12:
            continue
        coef = "" if (abs(abs(a) - 1) < 1e-12 and p > 0) else f"{abs(a):g}"
        term = coef + (var if p >= 1 else "") + (str(p).translate(SUP) if p > 1 else "")
        out += ("-" if a < 0 else "") + term if not out else (" - " if a < 0 else " + ") + term
    return out or "0"


def tf_block(label, n, d):
    ns, ds = poly_str(n), poly_str(d)
    w = max(len(ns), len(ds)) + 2
    pre = f"{label} = "
    return "\n".join([" " * len(pre) + ns.center(w), pre + "\u2500" * w, " " * len(pre) + ds.center(w)])


# ---------------------------------------------------------------- base tab
class Tab(ctk.CTkFrame):
    def __init__(self, master, get_system):
        super().__init__(master, fg_color="transparent")
        self.get_system = get_system
        self.side = ctk.CTkFrame(self, width=350, corner_radius=12)
        self.side.pack(side="left", fill="y", padx=(0, 8), pady=4)
        self.side.grid_propagate(False)
        right = ctk.CTkFrame(self, corner_radius=12)
        right.pack(side="right", fill="both", expand=True, pady=4)
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, right)
        NavigationToolbar2Tk(self.canvas, right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=6, pady=6)
        self.out = ctk.CTkTextbox(self.side, width=330, font=("Consolas", 12), wrap="none", state="disabled")

    # --- small widget helpers (grid inside the side panel)
    def field(self, label, default, row):
        v = tk.StringVar(value=default)
        ctk.CTkLabel(self.side, text=label).grid(row=row, column=0, sticky="w", padx=(12, 4), pady=4)
        ctk.CTkEntry(self.side, textvariable=v, width=100).grid(row=row, column=1, padx=(4, 12), pady=4)
        return v

    def option(self, label, values, row):
        v = tk.StringVar(value=values[0])
        ctk.CTkLabel(self.side, text=label).grid(row=row, column=0, sticky="w", padx=(12, 4), pady=4)
        ctk.CTkOptionMenu(self.side, values=values, variable=v, width=100).grid(row=row, column=1, padx=(4, 12), pady=4)
        return v

    def button(self, text, cmd, row):
        ctk.CTkButton(self.side, text=text, command=cmd).grid(row=row, column=0, columnspan=2,
                                                            sticky="ew", padx=12, pady=4)

    def place_out(self, row):
        self.out.grid(row=row, column=0, columnspan=2, sticky="nsew", padx=8, pady=8)
        self.side.grid_rowconfigure(row, weight=1)

    def show(self, lines):
        self.out.configure(state="normal")
        self.out.delete("1.0", "end")
        self.out.insert("end", "\n".join(lines))
        self.out.configure(state="disabled")
        self.canvas.draw()

    def guard(self, fn):
        try:
            fn()
        except Exception as e:  # show friendly error
            messagebox.showerror("Error", str(e))

    def calculate(self):
        self.guard(self.run)


# ---------------------------------------------------------------- time response
class TimeTab(Tab):
    def __init__(self, master, get_system, set_tf):
        super().__init__(master, get_system)
        self.set_tf = set_tf
        self.kind = self.option("Input signal", ["Step", "Ramp", "Impulse"], 0)
        self.tf = self.field("t final (blank=auto)", "", 1)
        self.button("Calculate & Plot", self.calculate, 2)
        self.zeta = self.field("2nd order: \u03b6", "0.5", 3)
        self.wn = self.field("wn (rad/s)", "5", 4)
        self.button("Load standard 2nd-order TF", self.load2, 5)
        self.place_out(6)

    def load2(self):
        z, w = float(self.zeta.get()), float(self.wn.get())
        self.set_tf(f"{w**2:g}", f"1 {2*z*w:g} {w**2:g}", unity=False)
        self.calculate()

    def run(self):
        n, d, unity = self.get_system()
        cn, cd = closed_loop(n, d, unity)
        sys = signal.lti(cn, cd)
        poles, zeros = np.roots(cd), np.roots(cn)
        stable = bool(np.all(poles.real < 0))
        T = np.linspace(0, float(self.tf.get()), 2000) if self.tf.get().strip() else None

        t1, y1 = signal.step(sys, T=T)
        kind = self.kind.get()
        if kind == "Step":
            t, y = t1, y1
        elif kind == "Impulse":
            t, y = signal.impulse(sys, T=T)
        else:
            t = t1
            _, y, _ = signal.lsim(sys, t, t)

        # ---- plots
        self.fig.clear()
        a1, a2 = self.fig.subplots(1, 2, gridspec_kw={"width_ratios": [2, 1]})
        a1.plot(t, y, lw=2, label="Output")
        if kind == "Ramp":
            a1.plot(t, t, "--", label="Input")
        elif kind == "Step":
            a1.axhline(1, ls="--", c="gray", lw=1, label="Input")
        a1.set(title=f"{kind} response", xlabel="Time (s)", ylabel="c(t)")
        a1.grid(True)

        lines = [f"Closed-loop poles: {np.round(poles, 4).tolist()}",
                 f"Zeros: {np.round(zeros, 4).tolist() if zeros.size else 'none'}",
                 f"Stable: {'YES' if stable else 'NO'}", ""]

        # ---- step-response specs
        if stable:
            ss = cn[-1] / cd[-1]
            lines.append(f"Final value (step): {fmt(ss)}")
            if abs(ss) > 1e-9:
                yn = y1 / ss
                first = lambda c: t1[np.flatnonzero(c)[0]] if np.any(c) else np.nan
                tr = first(yn >= 0.9) - first(yn >= 0.1)
                ip = int(np.argmax(yn))
                tp, Mp = t1[ip], max(0.0, (yn[ip] - 1) * 100)
                out = np.flatnonzero(np.abs(yn - 1) > 0.02)
                ts = 0.0 if out.size == 0 else (np.nan if out[-1] == len(t1) - 1 else t1[out[-1] + 1])
                lines += ["--- Step specs (numerical) ---",
                          f"Rise time tr (10-90%): {fmt(tr)} s",
                          f"Peak time tp: {fmt(tp)} s",
                          f"Peak overshoot Mp: {fmt(Mp)} %",
                          f"Settling time ts (2%): {fmt(ts)} s"]
                if kind == "Step":
                    a1.plot(tp, yn[ip] * ss, "ro")
                    if np.isfinite(ts):
                        a1.axvline(ts, ls=":", c="r", label=f"ts={ts:.3g}s")
        a1.legend()

        # ---- analytic 2nd-order formulas (no zeros)
        if cn.size == 1 and cd.size == 3:
            a = cd / cd[0]
            wn = np.sqrt(a[2])
            z = a[1] / (2 * wn)
            lines += ["", "--- 2nd-order formulas ---", f"wn = {fmt(wn)} rad/s", f"\u03b6 = {fmt(z)}"]
            if 0 < z < 1:
                wd = wn * np.sqrt(1 - z**2)
                lines += [f"wd = {fmt(wd)} rad/s",
                          f"tr = {fmt((np.pi - np.arccos(z)) / wd)} s",
                          f"tp = {fmt(np.pi / wd)} s",
                          f"Mp = {fmt(100 * np.exp(-np.pi * z / np.sqrt(1 - z**2)))} %",
                          f"ts(2%) = {fmt(4 / (z * wn))} s"]
            else:
                lines.append("(Undamped/critically/over-damped: no overshoot formulas)" if z >= 1 else "(\u03b6<=0)")

        # ---- steady-state error (unity feedback)
        if unity:
            N, (Kp, Kv, Ka) = error_constants(n, d)
            ess = [1 / (1 + Kp), np.inf if Kv == 0 else 1 / Kv, np.inf if Ka == 0 else 1 / Ka]
            lines += ["", f"--- Steady-state error (Type {N}) ---",
                      f"Kp={fmt(Kp)}  Kv={fmt(Kv)}  Ka={fmt(Ka)}",
                      f"ess step = {fmt(ess[0])}", f"ess ramp = {fmt(ess[1])}", f"ess parabolic = {fmt(ess[2])}"]

        # ---- pole-zero map
        a2.plot(poles.real, poles.imag, "bx", ms=9, mew=2, label="poles")
        if zeros.size:
            a2.plot(zeros.real, zeros.imag, "go", mfc="none", ms=8, label="zeros")
        a2.axhline(0, c="k", lw=.5)
        a2.axvline(0, c="k", lw=.5)
        a2.set(title="Pole-zero map", xlabel="Real", ylabel="Imag")
        a2.grid(True)
        a2.legend()
        self.fig.tight_layout()
        self.show(lines)


# ---------------------------------------------------------------- frequency response
class FreqTab(Tab):
    def __init__(self, master, get_system):
        super().__init__(master, get_system)
        self.button("Calculate & Plot", self.calculate, 0)
        self.place_out(1)

    def run(self):
        n, d, unity = self.get_system()
        r = np.abs(np.concatenate([np.roots(n), np.roots(d)]))
        r = r[r > 1e-9]
        lo = 10 ** (np.floor(np.log10(r.min())) - 2) if r.size else 1e-2
        hi = 10 ** (np.ceil(np.log10(r.max())) + 2) if r.size else 1e2
        w = np.logspace(np.log10(lo), np.log10(hi), 4000)
        G = np.polyval(n, 1j * w) / np.polyval(d, 1j * w)
        mag = 20 * np.log10(np.abs(G))
        ph = np.degrees(np.unwrap(np.angle(G)))

        wgc, wpc = crossings(w, mag, 0), crossings(w, ph, -180)
        lw = np.log10(w)
        PM = 180 + np.interp(np.log10(wgc[0]), lw, ph) if wgc else None
        GM = -np.interp(np.log10(wpc[0]), lw, mag) if wpc else None

        # ---- plots
        self.fig.clear()
        gs = self.fig.add_gridspec(2, 2)
        am, ap, an = self.fig.add_subplot(gs[0, 0]), self.fig.add_subplot(gs[1, 0]), self.fig.add_subplot(gs[:, 1])
        am.semilogx(w, mag)
        am.axhline(0, c="k", lw=.8)
        am.set(title="Bode magnitude", ylabel="dB")
        ap.semilogx(w, ph)
        ap.axhline(-180, c="k", lw=.8)
        ap.set(title="Bode phase", ylabel="deg", xlabel="w (rad/s)")
        for ax in (am, ap):
            ax.grid(True, which="both", alpha=.4)
            if wgc:
                ax.axvline(wgc[0], ls="--", c="r")
            if wpc:
                ax.axvline(wpc[0], ls="--", c="g")
        Gp = np.where(np.abs(G) < 20, G, np.nan)
        an.plot(Gp.real, Gp.imag, "b", label="w>0")
        an.plot(Gp.real, -Gp.imag, "b--", label="w<0")
        an.plot(-1, 0, "rx", ms=10, mew=2, label="-1+j0")
        an.set(title="Nyquist plot", xlabel="Real", ylabel="Imag")
        an.grid(True)
        an.legend()
        self.fig.tight_layout()

        # ---- numbers
        lines = ["--- Open-loop G(jw) ---",
                 f"Gain crossover wgc: {fmt(wgc[0]) if wgc else 'none'} rad/s",
                 f"Phase margin PM: {fmt(PM) + ' deg' if PM is not None else 'inf'}",
                 f"Phase crossover wpc: {fmt(wpc[0]) if wpc else 'none'} rad/s",
                 f"Gain margin GM: {fmt(GM) + ' dB' if GM is not None else 'inf'}"]
        cn, cd = closed_loop(n, d, unity)
        pl = np.roots(cd)
        lines.append(f"Closed-loop stable: {'YES' if np.all(pl.real < 0) else 'NO'}")
        if wgc and wpc:
            lines.append("Margins say: " + ("stable" if GM > 0 and PM > 0 else "unstable"))

        T = G / (1 + G) if unity else G
        Td = 20 * np.log10(np.abs(T))
        dc = Td[0]
        i = int(np.argmax(Td))
        low = np.flatnonzero(Td < dc - 3)
        lines += ["", "--- Closed-loop response ---" if unity else "--- System response ---",
                  f"Low-freq gain: {fmt(dc)} dB",
                  f"Bandwidth (-3 dB): {fmt(w[low[0]]) + ' rad/s' if low.size else '> plotted range'}"]
        if Td[i] - dc > 0.01:
            Mr = 10 ** ((Td[i] - dc) / 20)
            lines += [f"Resonant peak Mr: {fmt(Mr)} ({fmt(Td[i]-dc)} dB)", f"Resonant freq wr: {fmt(w[i])} rad/s"]
            if Mr > 1:
                z = np.sqrt((1 - np.sqrt(1 - 1 / Mr**2)) / 2)
                lines.append(f"Equiv. 2nd-order \u03b6 ~ {fmt(z)}")
        else:
            lines.append("No resonant peak (Mr = 1)")
        self.show(lines)


# ---------------------------------------------------------------- extra helpers
def stable(n, d, K=1.0):
    return bool(np.all(np.roots(np.polyadd(d, K * n)).real < -1e-9))


def marginal_gain(n, d, kmax=1e4):
    """(Ku, wu): gain and frequency where the locus first crosses the jw axis."""
    Ks = np.logspace(-4, np.log10(kmax), 1500)
    mr = lambda K: np.roots(np.polyadd(d, K * n)).real.max()
    v = np.array([mr(K) for K in Ks])
    i = np.flatnonzero((v[:-1] < 0) & (v[1:] >= 0))
    if not i.size:
        return None
    lo, hi = Ks[i[0]], Ks[i[0] + 1]
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if mr(mid) < 0 else (lo, mid)
    r = np.roots(np.polyadd(d, hi * n))
    return hi, abs(r[np.argmax(r.real)].imag)


def stable_ranges(n, d, kmax=1e4):
    Ks = np.concatenate([[0], np.logspace(-4, np.log10(kmax), 2000)])
    st = np.array([stable(n, d, K) for K in Ks])
    edges = []
    for i in np.flatnonzero(st[:-1] != st[1:]):
        lo, hi, s0 = Ks[i], Ks[i + 1], st[i]
        for _ in range(50):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if stable(n, d, mid) == s0 else (lo, mid)
        edges.append(hi)
    out, start = [], (0.0 if st[0] else None)
    for e in edges:
        if start is None:
            start = e
        else:
            out.append((start, e))
            start = None
    if start is not None:
        out.append((start, np.inf))
    return out


def routh(p, eps=1e-6):
    p = np.asarray(p, float)
    n = len(p)
    if n < 2:
        raise ValueError("Need a polynomial of order >= 1")
    w = int(np.ceil(n / 2))
    R = np.zeros((n, w))
    R[0, :len(p[0::2])] = p[0::2]
    R[1, :len(p[1::2])] = p[1::2]
    notes = []
    for i in range(2, n):
        if np.allclose(R[i - 1], 0):
            deg = n - 1 - (i - 2)
            pw = deg - 2 * np.arange(w)
            R[i - 1] = R[i - 2] * np.where(pw > 0, pw, 0)
            notes.append(f"Row s^{n - i} was all zeros: auxiliary polynomial A(s) from row s^{deg} differentiated "
                         "(roots symmetric about origin -> jw-axis or mirrored RHP roots).")
        if abs(R[i - 1, 0]) < 1e-12:
            R[i - 1, 0] = eps
            notes.append(f"Zero first element in row s^{n - i}: replaced by epsilon={eps:g}.")
        for j in range(w - 1):
            R[i, j] = (R[i - 1, 0] * R[i - 2, j + 1] - R[i - 2, 0] * R[i - 1, j + 1]) / R[i - 1, 0]
    return R, notes


def freq_grid(n, d, pts=4000):
    r = np.abs(np.concatenate([np.roots(n), np.roots(d)]))
    r = r[r > 1e-9]
    lo = 10 ** (np.floor(np.log10(r.min())) - 2) if r.size else 1e-2
    hi = 10 ** (np.ceil(np.log10(r.max())) + 2) if r.size else 1e2
    return np.logspace(np.log10(lo), np.log10(hi), pts)


def step_specs(cn, cd):
    poles = np.roots(cd)
    t, y = signal.step(signal.lti(cn, cd))
    r = {"t": t, "y": y, "stable": bool(np.all(poles.real < 0)), "poles": poles}
    if r["stable"]:
        ss = cn[-1] / cd[-1]
        r["ss"] = ss
        if abs(ss) > 1e-9:
            yn = y / ss
            first = lambda c: t[np.flatnonzero(c)[0]] if np.any(c) else np.nan
            ip = int(np.argmax(yn))
            out = np.flatnonzero(np.abs(yn - 1) > .02)
            r.update(tr=first(yn >= .9) - first(yn >= .1), tp=t[ip], Mp=max(0., (yn[ip] - 1) * 100),
                     ts=0.0 if not out.size else (np.nan if out[-1] == len(t) - 1 else t[out[-1] + 1]))
    return r


def spec_lines(name, r):
    L = [f"[{name}] stable: {'YES' if r['stable'] else 'NO'}"]
    if r["stable"]:
        L.append(f"  final value {fmt(r['ss'])}   ess {fmt(1 - r['ss'])}")
        if "tr" in r:
            L.append(f"  tr {fmt(r['tr'])}s  tp {fmt(r['tp'])}s  Mp {fmt(r['Mp'])}%  ts {fmt(r['ts'])}s")
    return L


# ---------------------------------------------------------------- root locus
class RootLocusTab(Tab):
    def __init__(self, master, get_system):
        super().__init__(master, get_system)
        self.kmax = self.field("K max (blank=auto)", "", 0)
        self.K = self.field("Evaluate at K", "", 1)
        self.zt = self.field("Target \u03b6", "0.5", 2)
        self.button("Plot root locus", self.calculate, 3)
        self.place_out(4)

    def run(self):
        n, d, _ = self.get_system()   # open-loop G(s) = n/d, characteristic eq: d + K n = 0
        P, Z = len(d) - 1, len(n) - 1
        if P < 1:
            raise ValueError("Denominator must have order >= 1")
        poles, zeros = np.roots(d), np.roots(n)
        mg = marginal_gain(n, d)
        ext = max(1.0, np.abs(np.concatenate([poles, zeros])).max()) * 2.5
        if self.kmax.get().strip():
            kmax = float(self.kmax.get())
        elif mg:
            kmax = mg[0] * 4
        else:
            kmax = 10 * abs(np.polyval(d, 2j * ext) / np.polyval(n, 2j * ext))
        Ks = np.concatenate([[0], np.logspace(np.log10(kmax) - 5, np.log10(kmax), 1500)])
        R = np.zeros((len(Ks), P), complex)
        R[0] = prev = poles
        for i in range(1, len(Ks)):
            r = np.roots(np.polyadd(d, Ks[i] * n))
            if len(r) != P:
                raise ValueError("Root count changed - use a strictly proper G(s)")
            ri, ci = linear_sum_assignment(np.abs(prev[:, None] - r[None, :]))
            R[i, ri] = r[ci]
            prev = R[i]

        self.fig.clear()
        ax = self.fig.subplots()
        for b in range(P):
            ax.plot(R[:, b].real, R[:, b].imag, lw=1.3)
        ax.plot(poles.real, poles.imag, "kx", ms=9, mew=2, label="open-loop poles")
        if zeros.size:
            ax.plot(zeros.real, zeros.imag, "ko", mfc="none", ms=8, label="zeros")
        ax.axhline(0, c="k", lw=.5)
        ax.axvline(0, c="k", lw=.5)
        lines = [f"Open-loop poles: {np.round(poles, 4).tolist()}",
                 f"Zeros: {np.round(zeros, 4).tolist() if zeros.size else 'none'}",
                 f"Branches: {P}   Asymptotes: {max(P - Z, 0)}"]
        if P > Z:
            ang = (2 * np.arange(P - Z) + 1) * 180 / (P - Z)
            cen = (poles.real.sum() - zeros.real.sum()) / (P - Z)
            lines += [f"Asymptote angles: {np.round(ang, 2).tolist()} deg", f"Centroid: {fmt(cen)}"]
            if P - Z >= 2:
                for a in np.radians(ang):
                    ax.plot([cen, cen + 3 * ext * np.cos(a)], [0, 3 * ext * np.sin(a)], "k:", lw=.8)
        dp = np.trim_zeros(np.polysub(np.polymul(np.polyder(d), n), np.polymul(d, np.polyder(n))), "f")
        bk = []
        for s in (np.roots(dp) if dp.size > 1 else []):
            if abs(s.imag) < 1e-7:
                K = -np.polyval(d, s.real) / np.polyval(n, s.real)
                if K > 0:
                    bk.append((s.real, K))
        for s, K in bk:
            ax.plot(s, 0, "gD", ms=7)
            lines.append(f"Breakaway/break-in: s = {fmt(s)}  (K = {fmt(K)})")
        if mg:
            ax.plot([0, 0], [mg[1], -mg[1]], "r*", ms=11, label="jw crossing")
            lines += [f"jw-axis crossing: K = {fmt(mg[0])}, w = {fmt(mg[1])} rad/s",
                      f"Stable for 0 < K < {fmt(mg[0])}"]
        else:
            lines.append("No jw-axis crossing found (stable for all K tested)")

        zt = float(self.zt.get()) if self.zt.get().strip() else None
        if zt and 0 < zt < 1:
            a = np.arccos(zt)
            for sg in (1, -1):
                ax.plot([0, -ext * 3 * np.cos(a)], [0, sg * ext * 3 * np.sin(a)], "m--", lw=.9)
            with np.errstate(all="ignore"):
                z = np.where(R.imag > 1e-6, -R.real / np.abs(R), np.inf)
            zm = np.minimum(z.min(axis=1), 1.0)
            for K in crossings(Ks[1:], zm[1:], zt)[:3]:
                r = np.roots(np.polyadd(d, K * n))
                c = r[np.argmax(r.imag)]
                ax.plot(c.real, c.imag, "ms", ms=8)
                ax.plot(c.real, -c.imag, "ms", ms=8)
                lines.append(f"\u03b6={zt:g}: K = {fmt(K)}, dominant poles {fmt(c.real)} +/- j{fmt(c.imag)}")

        if self.K.get().strip():
            K = float(self.K.get())
            r = np.roots(np.polyadd(d, K * n))
            ax.plot(r.real, r.imag, "rs", mfc="none", ms=11, mew=2, label=f"poles at K={K:g}")
            lines += ["", f"--- Closed-loop poles at K = {K:g} ---"]
            for p in r:
                wn = abs(p)
                lines.append(f"  {np.round(p, 4)}   wn={fmt(wn)}  \u03b6={fmt(-p.real / wn)}")
            lines.append(f"Stable: {'YES' if np.all(r.real < 0) else 'NO'}")
        ax.set(title="Root locus (K from 0 upward)", xlabel="Real", ylabel="Imag",
               xlim=(-ext, ext), ylim=(-ext, ext))
        ax.grid(True)
        ax.legend(loc="upper left", fontsize=8)
        self.fig.tight_layout()
        self.show(lines)


# ---------------------------------------------------------------- Routh-Hurwitz
class RouthTab(Tab):
    def __init__(self, master, get_system):
        super().__init__(master, get_system)
        self.K = self.field("Gain K (unity fb)", "1", 0)
        self.button("Build Routh table", self.calculate, 1)
        self.place_out(2)

    def run(self):
        n, d, unity = self.get_system()
        K = float(self.K.get())
        p = np.trim_zeros(np.polyadd(d, K * n) if unity else d, "f")
        R, notes = routh(p)
        deg = len(p) - 1
        lines = [("Characteristic eq: den + K*num = 0" if unity else "Characteristic eq: denominator = 0")
                 + (f"  (K={K:g})" if unity else ""),
                 f"Coefficients: {np.round(p, 5).tolist()}", "", "Routh array:"]
        for i, row in enumerate(R):
            lines.append(f"s^{deg - i}: " + " ".join(f"{x:10.4g}" for x in row))
        roots = np.roots(p)
        rhp, jw = int(np.sum(roots.real > 1e-9)), int(np.sum(np.abs(roots.real) <= 1e-9))
        sc = int(np.sum(np.sign(R[:-1, 0]) * np.sign(R[1:, 0]) < 0))
        lines += ["", f"Sign changes in 1st column: {sc}", f"Roots in RHP: {rhp}   on jw-axis: {jw}",
                  "Verdict: " + ("STABLE" if rhp == 0 and jw == 0 else
                                 "MARGINALLY STABLE" if rhp == 0 else "UNSTABLE")]
        lines += [""] + notes if notes else []
        if unity:
            rg = stable_ranges(n, d)
            lines += ["", "Stable gain range (scan 0..1e4):"]
            lines += [f"  {('0' if a < 2e-4 else fmt(a))} < K < {fmt(b)}" for a, b in rg] or ["  none"]
        lines += ["", f"Actual roots: {np.round(roots, 4).tolist()}"]
        self.fig.clear()
        ax = self.fig.subplots()
        ax.axvspan(0, max(1, np.abs(roots).max() * 1.5), color="red", alpha=.08)
        ax.plot(roots.real, roots.imag, "bx", ms=10, mew=2)
        ax.axhline(0, c="k", lw=.5)
        ax.axvline(0, c="k", lw=.5)
        ax.set(title="Roots of characteristic equation (shaded = RHP)", xlabel="Real", ylabel="Imag")
        ax.grid(True)
        self.fig.tight_layout()
        self.show(lines)


# ---------------------------------------------------------------- Nichols chart
class NicholsTab(Tab):
    def __init__(self, master, get_system):
        super().__init__(master, get_system)
        self.tgt = self.field("Target Mr (dB)", "3", 0)
        self.button("Plot Nichols chart", self.calculate, 1)
        self.place_out(2)

    def run(self):
        n, d, _ = self.get_system()   # loop transfer function L(s) = n/d, unity feedback
        w = freq_grid(n, d)
        G = np.polyval(n, 1j * w) / np.polyval(d, 1j * w)
        mag = 20 * np.log10(np.abs(G))
        ph = np.degrees(np.unwrap(np.angle(G)))
        xp = (ph % 360) - 360
        j = np.flatnonzero(np.abs(np.diff(xp)) > 180) + 1
        xp, yp = np.insert(xp, j, np.nan), np.insert(mag, j, np.nan)

        self.fig.clear()
        ax = self.fig.subplots()
        thd = np.linspace(-359.5, -0.5, 720)
        th = np.radians(thd)
        for Mdb in (-20, -12, -6, -3, -1, 0, 1, 3, 6, 12):
            M2 = 10 ** (Mdb / 10)
            with np.errstate(invalid="ignore", divide="ignore"):
                if Mdb == 0:
                    rs = [np.where(np.cos(th) < 0, -1 / (2 * np.cos(th)), np.nan)]
                else:
                    a, b = M2 - 1, 2 * M2 * np.cos(th)
                    disc = b * b - 4 * a * M2
                    rs = [np.where((disc >= 0) & ((r := (-b + s * np.sqrt(disc)) / (2 * a)) > 0), r, np.nan)
                          for s in (1, -1)]
                for r in rs:
                    y = 20 * np.log10(r)
                    ax.plot(thd, y, c="tab:red", lw=.7, alpha=.6)
                    k = np.flatnonzero(~np.isnan(y))
                    if k.size and r is rs[0]:
                        ax.text(thd[k[0]], y[k[0]], f"{Mdb}dB", fontsize=7, color="tab:red")
        for al in (-5, -15, -30, -45, -60, -90, -120, -150):
            a = np.radians(al)
            r = np.sin(th) * np.cos(a) / np.sin(a) - np.cos(th)
            with np.errstate(invalid="ignore", divide="ignore"):
                y = 20 * np.log10(np.where(r > 0, r, np.nan))
            ax.plot(thd, y, c="tab:green", lw=.7, alpha=.6)
            k = np.flatnonzero(~np.isnan(y))
            if k.size:
                ax.text(thd[k[-1]], y[k[-1]], f"{al}\u00b0", fontsize=7, color="tab:green")
        ax.plot(xp, yp, "b", lw=2, label="L(jw)")
        ax.plot(-180, 0, "rx", ms=10, mew=2, label="critical point")
        ax.set(title="Nichols chart (red: M-contours, green: N-contours)", xlabel="Phase (deg)",
               ylabel="Gain (dB)", xlim=(-360, 0),
               ylim=(min(-40, max(-100, mag.min())), max(20, min(60, mag.max()))))
        ax.grid(True, alpha=.4)
        ax.legend(loc="lower left")
        self.fig.tight_layout()

        wgc, wpc = crossings(w, mag, 0), crossings(w, ph, -180)
        lw_ = np.log10(w)
        PM = 180 + np.interp(np.log10(wgc[0]), lw_, ph) if wgc else None
        GM = -np.interp(np.log10(wpc[0]), lw_, mag) if wpc else None
        Td = 20 * np.log10(np.abs(G / (1 + G)))
        i = int(np.argmax(Td))
        low = np.flatnonzero(Td < Td[0] - 3)
        peak = (f"Closed-loop peak Mr: {fmt(10 ** (Td[i] / 20))} ({fmt(Td[i])} dB) at wr = {fmt(w[i])} rad/s"
                if i > 0 and Td[i] > Td[0] + .01 else "Closed-loop: no resonant peak (Mr = 1)")
        lines = [f"Gain margin: {fmt(GM) + ' dB' if GM is not None else 'inf'}",
                 f"Phase margin: {fmt(PM) + ' deg' if PM is not None else 'inf'}",
                 f"Closed-loop stable: {'YES' if stable(n, d) else 'NO'}", "",
                 peak,
                 f"Bandwidth (-3 dB): {fmt(w[low[0]]) + ' rad/s' if low.size else '> range'}"]
        if self.tgt.get().strip():
            tg = float(self.tgt.get())
            Ks = np.logspace(-3, 3, 500)
            with np.errstate(all="ignore"):
                mr = np.array([(20 * np.log10(np.abs(K * G / (1 + K * G)))).max() for K in Ks])
            cr = crossings(Ks, mr, tg)
            cr = [k for k in cr if stable(n, d, k)]
            lines += ["", f"--- Gain adjustment for Mr = {tg:g} dB (stable loop) ---"]
            lines += ([f"Multiply loop gain by K = {fmt(k)}" for k in cr[:3]] or
                      ["No stable gain in 0.001..1000 gives that Mr"])
        self.show(lines)


# ---------------------------------------------------------------- PID tuning
class PIDTab(Tab):
    def __init__(self, master, get_system):
        super().__init__(master, get_system)
        self.Kp, self.Ki = self.field("Kp", "1", 0), self.field("Ki", "0", 1)
        self.Kd, self.N = self.field("Kd", "0", 2), self.field("Deriv. filter N", "20", 3)
        self.button("Simulate PID", self.calculate, 4)
        for r, k in enumerate(("P", "PI", "PID"), 5):
            self.button(f"Ziegler-Nichols {k} (auto)", lambda k=k: self.guard(lambda: self.zn(k)), r)
        self.place_out(8)

    def zn(self, kind):
        n, d, _ = self.get_system()
        mg = marginal_gain(n, d)
        if not mg:
            raise ValueError("Plant never reaches marginal stability with P control - Z-N not applicable")
        Ku, wu = mg
        Pu = 2 * np.pi / wu
        kp, ti, td = {"P": (.5 * Ku, np.inf, 0), "PI": (.45 * Ku, Pu / 1.2, 0), "PID": (.6 * Ku, Pu / 2, Pu / 8)}[kind]
        self.Kp.set(f"{kp:.4g}")
        self.Ki.set(f"{kp / ti:.4g}")
        self.Kd.set(f"{kp * td:.4g}")
        self.run(extra=[f"Z-N: Ku={fmt(Ku)}, Pu={fmt(Pu)} s"])

    def run(self, extra=()):
        n, d, _ = self.get_system()   # plant G(s) = n/d
        Kp, Ki, Kd, N = (float(v.get()) for v in (self.Kp, self.Ki, self.Kd, self.N))
        tau = Kd / N if Kd > 0 else 0.0
        f = np.trim_zeros(np.array([tau, 1.0]), "f")           # (tau s + 1)
        if Ki != 0:
            dc = np.polymul(f, [1.0, 0.0])
            nc = np.polyadd(np.polyadd(Kp * dc, Ki * f), Kd * np.array([1.0, 0, 0]))
        else:
            dc, nc = f, np.polyadd(Kp * f, Kd * np.array([1.0, 0]))
        nc = np.trim_zeros(nc, "f")
        if nc.size == 0:
            raise ValueError("All controller gains are zero")
        Ln, Ld = np.polymul(nc, n), np.polymul(dc, d)
        ru = step_specs(n, np.trim_zeros(np.polyadd(d, n), "f"))
        rp = step_specs(Ln, np.trim_zeros(np.polyadd(Ld, Ln), "f"))
        self.fig.clear()
        ax = self.fig.subplots()
        ax.plot(ru["t"], ru["y"], "--", label="P-only, K=1 (plant + unity fb)")
        ax.plot(rp["t"], rp["y"], lw=2, label=f"PID  Kp={Kp:g} Ki={Ki:g} Kd={Kd:g}")
        ax.axhline(1, c="gray", lw=.8)
        ax.set(title="Closed-loop step response", xlabel="Time (s)", ylabel="y(t)")
        if not (ru["stable"] and rp["stable"]):
            ax.set_ylim(-1, 3)
        ax.grid(True)
        ax.legend(fontsize=8)
        self.fig.tight_layout()
        lines = list(extra) + [f"C(s) = Kp + Ki/s + Kd*s/(Kd/N*s+1)", ""]
        lines += spec_lines("Without PID", ru) + [""] + spec_lines("With PID", rp)
        lines += ["", f"PID closed-loop poles: {np.round(rp['poles'], 3).tolist()}"]
        self.show(lines)


# ---------------------------------------------------------------- compensator design
def margins(w, G):
    mag = 20 * np.log10(np.abs(G))
    ph = np.degrees(np.unwrap(np.angle(G)))
    lw = np.log10(w)
    gc, pc = crossings(w, mag, 0), crossings(w, ph, -180)
    PM = 180 + np.interp(np.log10(gc[0]), lw, ph) if gc else None
    GM = -np.interp(np.log10(pc[0]), lw, mag) if pc else None
    return mag, ph, (gc[0] if gc else None), PM, (pc[0] if pc else None), GM


def _pm_after(w, KG, nc, dc):
    return margins(w, KG * np.polyval(nc, 1j * w) / np.polyval(dc, 1j * w))[3]


def design_lead(w, KG, pm_des, margin):
    mag, ph, wgc, PM, _, _ = margins(w, KG)
    if PM is None:
        raise ValueError("K*G(s) has no gain crossover - adjust K")
    if PM >= pm_des:
        raise ValueError(f"Phase margin is already {PM:.1f} deg (>= {pm_des:g}); no lead needed")
    for extra in np.arange(0, 40, 1.0):          # raise the phase allowance until the PM is met
        phim = np.radians(pm_des - PM + margin + extra)
        if phim >= np.radians(70):
            raise ValueError("Lead needs more than 70 deg: use a lag-lead or two lead stages")
        a = (1 - np.sin(phim)) / (1 + np.sin(phim))
        cr = crossings(w, mag, 10 * np.log10(a))
        if not cr:
            continue
        wm = cr[0]
        T = 1 / (wm * np.sqrt(a))
        nc, dc = np.array([T, 1.0]), np.array([a * T, 1.0])
        if (_pm_after(w, KG, nc, dc) or 0) >= pm_des:
            break
    else:
        raise ValueError("Could not reach the phase margin with a single lead stage")
    return nc, dc, [f"alpha = {fmt(a)}   T = {fmt(T)} s", f"phi_max = {fmt(np.degrees(phim))} deg at wm = {fmt(wm)} rad/s",
                    f"zero at s = -{fmt(1 / T)},  pole at s = -{fmt(1 / (a * T))}"]


def design_lag(w, KG, pm_des, margin):
    mag, ph, wgc, PM, _, _ = margins(w, KG)
    cr = crossings(w, ph, -180 + pm_des + margin)
    if not cr:
        raise ValueError("Phase never reaches the required value - a lag compensator cannot give this PM")
    wc, lw = cr[0], np.log10(w)
    g = np.interp(np.log10(wc), lw, mag)
    if g <= 0:
        raise ValueError(f"Gain at the new crossover is already {g:.1f} dB - no lag needed")
    beta, T = 10 ** (g / 20), 10 / wc
    return np.array([T, 1.0]), np.array([beta * T, 1.0]), [
        f"new crossover wc' = {fmt(wc)} rad/s  (|KG| = {fmt(g)} dB there)",
        f"beta = {fmt(beta)}   T = {fmt(T)} s",
        f"zero at s = -{fmt(1 / T)},  pole at s = -{fmt(1 / (beta * T))}"]


def design_laglead(w, KG, pm_des, margin, w1=None):
    mag, ph, *_ = margins(w, KG)
    lw = np.log10(w)
    if w1 is None:
        cr = crossings(w, ph, -180)
        if not cr:
            raise ValueError("No phase crossover found - enter the new crossover w1 manually")
        w1 = cr[0]
    phi1, g1 = np.interp(np.log10(w1), lw, ph), np.interp(np.log10(w1), lw, mag)
    lag_ph = 0.0
    for _ in range(40):                              # lead phase must also cancel the lag's phase
        phil = pm_des + margin - 180 - phi1 + abs(lag_ph)
        if not 0 < phil < 70:
            raise ValueError(f"Lead section would need {phil:.1f} deg at w1={w1:.3g}: pick another w1")
        sn = np.sin(np.radians(phil))
        gam = (1 + sn) / (1 - sn)
        X = 101 * 10 ** ((g1 + 10 * np.log10(gam)) / 10)
        beta = np.sqrt(max((X - 1) / 100, 1.0))
        lag_ph = np.degrees(np.arctan(10) - np.arctan(10 * beta))
    T1, T2 = np.sqrt(gam) / w1, 10 / w1
    nc = np.polymul([T1, 1.0], [T2, 1.0])
    dc = np.polymul([T1 / gam, 1.0], [beta * T2, 1.0])
    return nc, dc, [f"new crossover w1 = {fmt(w1)} rad/s  (|KG| = {fmt(g1)} dB, phase {fmt(phi1)} deg)",
                    f"Lead part: gamma = {fmt(gam)}  T1 = {fmt(T1)} s  (phase {fmt(phil)} deg at w1)",
                    f"Lag part:  beta = {fmt(beta)}  T2 = {fmt(T2)} s  (phase {fmt(lag_ph)} deg at w1)"]


class CompensatorTab(Tab):
    def __init__(self, master, get_system):
        super().__init__(master, get_system)
        self.kind = self.option("Type", ["Lead", "Lag", "Lag-Lead"], 0)
        self.K = self.field("Gain K (for error spec)", "1", 1)
        self.pm = self.field("Desired PM (deg)", "45", 2)
        self.mg = self.field("Safety margin (deg)", "5", 3)
        self.w1 = self.field("w1 (lag-lead, blank=auto)", "", 4)
        self.button("Design compensator", self.calculate, 5)
        self.place_out(6)

    def run(self):
        n, d, _ = self.get_system()                 # plant G(s) = n/d, unity feedback
        K, pm, mg = float(self.K.get()), float(self.pm.get()), float(self.mg.get())
        w = freq_grid(K * n, d, 6000)
        KG = K * np.polyval(n, 1j * w) / np.polyval(d, 1j * w)
        kind = self.kind.get()
        if kind == "Lead":
            nc, dc, info = design_lead(w, KG, pm, mg)
        elif kind == "Lag":
            nc, dc, info = design_lag(w, KG, pm, mg)
        else:
            nc, dc, info = design_laglead(w, KG, pm, mg, float(self.w1.get()) if self.w1.get().strip() else None)

        Gc = np.polyval(nc, 1j * w) / np.polyval(dc, 1j * w)
        m0, p0, wg0, PM0, wp0, GM0 = margins(w, KG)
        m1, p1, wg1, PM1, wp1, GM1 = margins(w, KG * Gc)
        cn0, cd0 = K * n, np.trim_zeros(np.polyadd(d, K * n), "f")
        Ln, Ld = K * np.polymul(n, nc), np.polymul(d, dc)
        cd1 = np.trim_zeros(np.polyadd(Ld, Ln), "f")
        r0, r1 = step_specs(cn0, cd0), step_specs(Ln, cd1)

        self.fig.clear()
        gs = self.fig.add_gridspec(2, 2)
        am, ap, ast = self.fig.add_subplot(gs[0, 0]), self.fig.add_subplot(gs[1, 0]), self.fig.add_subplot(gs[:, 1])
        for ax, a0, a1, ac, ttl, yl in ((am, m0, m1, 20 * np.log10(np.abs(Gc)), "Bode magnitude", "dB"),
                                        (ap, p0, p1, np.degrees(np.unwrap(np.angle(Gc))), "Bode phase", "deg")):
            ax.semilogx(w, a0, "--", label="uncompensated")
            ax.semilogx(w, a1, lw=2, label="compensated")
            ax.semilogx(w, ac, ":", c="g", label="Gc(jw)")
            ax.set(title=ttl, ylabel=yl)
            ax.grid(True, which="both", alpha=.4)
        am.axhline(0, c="k", lw=.8)
        ap.axhline(-180, c="k", lw=.8)
        ap.set_xlabel("w (rad/s)")
        am.legend(fontsize=7)
        ast.plot(r0["t"], r0["y"], "--", label="uncompensated")
        ast.plot(r1["t"], r1["y"], lw=2, label="compensated")
        ast.axhline(1, c="gray", lw=.8)
        ast.set(title="Closed-loop step response", xlabel="Time (s)")
        if not (r0["stable"] and r1["stable"]):
            ast.set_ylim(-1, 3)
        ast.grid(True)
        ast.legend(fontsize=8)
        self.fig.tight_layout()

        N, (Kp, Kv, Ka) = error_constants(K * n, d)
        f = lambda v, u: (fmt(v) + u) if v is not None else "inf"
        lines = [f"--- {kind} compensator ---", "Gc(s) = (num)/(den), DC gain 1:", tf_block("Gc(s)", nc, dc), ""] + info
        lines += ["", f"Total gain applied: K = {fmt(K)}   (type {N}, Kp={fmt(Kp)}, Kv={fmt(Kv)}, Ka={fmt(Ka)})",
                  "", "            before     after",
                  f"PM        {f(PM0, ' deg'):>10} {f(PM1, ' deg'):>10}",
                  f"GM        {f(GM0, ' dB'):>10} {f(GM1, ' dB'):>10}",
                  f"w_gc      {f(wg0, ''):>10} {f(wg1, ''):>10}", ""]
        lines += spec_lines("Before", r0) + spec_lines("After", r1)
        if PM1 is not None and PM1 < pm - 0.5:
            lines += ["", f"Note: PM {PM1:.1f} < target {pm:g}; increase the safety margin."]
        self.show(lines)


# ---------------------------------------------------------------- main app
ctk.set_appearance_mode("System")
ctk.set_default_color_theme("blue")


class App(ctk.CTk):
    TABS = ["Time Response", "Frequency Response", "Root Locus", "Routh-Hurwitz", "Nichols Chart", "PID Tuning", "Compensator"]

    def __init__(self):
        super().__init__()
        self.title("Linear Control Systems Toolkit")
        self.geometry("1400x880")
        self.num, self.den = tk.StringVar(value="10"), tk.StringVar(value="1 2 0")
        self.unity = tk.BooleanVar(value=True)

        head = ctk.CTkFrame(self, corner_radius=12)
        head.pack(fill="x", padx=12, pady=(12, 6))
        head.grid_columnconfigure(1, weight=1)
        box = ctk.CTkFrame(head, fg_color="transparent")
        box.grid(row=0, column=0, padx=12, pady=10, sticky="n")
        ctk.CTkLabel(box, text="Linear Control Systems", font=ctk.CTkFont(size=20, weight="bold")
                     ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 6))
        for r, (lab, var) in enumerate((("Numerator", self.num), ("Denominator", self.den)), 1):
            ctk.CTkLabel(box, text=lab).grid(row=r, column=0, sticky="w", pady=3)
            ctk.CTkEntry(box, textvariable=var, width=230).grid(row=r, column=1, padx=(8, 0), pady=3)
            var.trace_add("write", lambda *a: self.update_tf())
        ctk.CTkCheckBox(box, text="Entered G(s) is open-loop, unity feedback", variable=self.unity,
                        command=self.update_tf).grid(row=3, column=0, columnspan=2, sticky="w", pady=6)
        ctk.CTkLabel(box, text="Coefficients, highest power first:  1 2 0  =  s\u00b2 + 2s",
                     text_color="gray").grid(row=4, column=0, columnspan=2, sticky="w")
        ctk.CTkOptionMenu(box, values=["System", "Light", "Dark"], width=100,
                          command=ctk.set_appearance_mode).grid(row=5, column=0, sticky="w", pady=(8, 0))

        self.display = ctk.CTkTextbox(head, height=190, font=("Consolas", 15), wrap="none", state="disabled")
        self.display.grid(row=0, column=1, padx=12, pady=12, sticky="nsew")

        self.tabs = ctk.CTkTabview(self)
        self.tabs.pack(fill="both", expand=True, padx=12, pady=(6, 12))
        makers = [lambda m: TimeTab(m, self.system, self.set_tf), lambda m: FreqTab(m, self.system),
                  lambda m: RootLocusTab(m, self.system), lambda m: RouthTab(m, self.system),
                  lambda m: NicholsTab(m, self.system), lambda m: PIDTab(m, self.system),
                  lambda m: CompensatorTab(m, self.system)]
        for name, mk in zip(self.TABS, makers):
            mk(self.tabs.add(name)).pack(fill="both", expand=True)
        self.update_tf()

    def set_tf(self, n, d, unity):
        self.num.set(n)
        self.den.set(d)
        self.unity.set(unity)
        self.update_tf()

    def system(self):
        n, d = parse(self.num.get()), parse(self.den.get())
        if n.size > d.size:
            raise ValueError("Improper transfer function (numerator order > denominator order)")
        return n, d, self.unity.get()

    def tf_text(self):
        n, d, unity = self.system()
        if unity:
            cn, cd = closed_loop(n, d, True)
            L = [tf_block("G(s)", n, d), "", "Closed loop  T(s) = G/(1+G):", tf_block("T(s)", cn, cd), "",
                 f"Characteristic eq:  {poly_str(cd)} = 0",
                 f"Open-loop poles:   {np.round(np.roots(d), 4).tolist()}",
                 f"Closed-loop poles: {np.round(np.roots(cd), 4).tolist()}"]
        else:
            L = [tf_block("T(s)", n, d), "", f"Poles: {np.round(np.roots(d), 4).tolist()}"]
            ol = np.trim_zeros(np.polysub(d, n), "f")
            if ol.size and ol.size >= n.size:
                L += ["", "Equivalent unity-feedback open loop  G = T/(1-T):", tf_block("G(s)", n, ol)]
        return "\n".join(L)

    def update_tf(self):
        try:
            txt = self.tf_text()
        except Exception as e:
            txt = f"  (waiting for valid input: {e})"
        self.display.configure(state="normal")
        self.display.delete("1.0", "end")
        self.display.insert("end", txt)
        self.display.configure(state="disabled")


if __name__ == "__main__":
    App().mainloop()