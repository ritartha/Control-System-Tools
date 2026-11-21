# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 14:14:14 2026

@author: ritar
"""

"""
Linear Control Systems Toolkit (Tkinter)
Tabs: Time Response | Frequency Response | Root Locus | Routh-Hurwitz | Nichols Chart | PID Tuning
Requires: numpy, scipy, matplotlib      ->  pip install numpy scipy matplotlib
Run:      python control_gui.py
"""
import tkinter as tk
from tkinter import ttk, messagebox
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


# ---------------------------------------------------------------- base tab
class Tab(ttk.Frame):
    def __init__(self, master, get_system):
        super().__init__(master)
        self.get_system = get_system
        self.side = ttk.Frame(self)
        self.side.pack(side="left", fill="y", padx=6, pady=6)
        right = ttk.Frame(self)
        right.pack(side="right", fill="both", expand=True)
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, right)
        NavigationToolbar2Tk(self.canvas, right)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.out = tk.Text(self.side, width=46, height=28, font=("Consolas", 9))

    def show(self, lines):
        self.out.config(state="normal")
        self.out.delete("1.0", "end")
        self.out.insert("end", "\n".join(lines))
        self.out.config(state="disabled")
        self.canvas.draw()

    def field(self, label, default, row):
        v = tk.StringVar(value=default)
        ttk.Label(self.side, text=label).grid(row=row, column=0, sticky="w")
        ttk.Entry(self.side, textvariable=v, width=10).grid(row=row, column=1)
        return v

    def button(self, text, cmd, row):
        ttk.Button(self.side, text=text, command=cmd).grid(row=row, column=0, columnspan=2, pady=3)

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
        s = self.side
        self.kind = tk.StringVar(value="Step")
        self.tf = tk.StringVar()
        self.zeta = tk.StringVar(value="0.5")
        self.wn = tk.StringVar(value="5")
        ttk.Label(s, text="Input signal").grid(row=0, column=0, sticky="w")
        ttk.Combobox(s, textvariable=self.kind, values=["Step", "Ramp", "Impulse"],
                     width=10, state="readonly").grid(row=0, column=1)
        ttk.Label(s, text="t final (blank=auto)").grid(row=1, column=0, sticky="w")
        ttk.Entry(s, textvariable=self.tf, width=10).grid(row=1, column=1)
        ttk.Button(s, text="Calculate & Plot", command=self.calculate).grid(row=2, column=0, columnspan=2, pady=4)
        ttk.Separator(s).grid(row=3, column=0, columnspan=2, sticky="ew", pady=3)
        ttk.Label(s, text="Standard 2nd order: zeta").grid(row=4, column=0, sticky="w")
        ttk.Entry(s, textvariable=self.zeta, width=10).grid(row=4, column=1)
        ttk.Label(s, text="wn (rad/s)").grid(row=5, column=0, sticky="w")
        ttk.Entry(s, textvariable=self.wn, width=10).grid(row=5, column=1)
        ttk.Button(s, text="Load wn^2/(s^2+2*zeta*wn*s+wn^2)", command=self.load2).grid(row=6, column=0, columnspan=2, pady=4)
        self.out.grid(row=7, column=0, columnspan=2)

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
            lines += ["", "--- 2nd-order formulas ---", f"wn = {fmt(wn)} rad/s", f"zeta = {fmt(z)}"]
            if 0 < z < 1:
                wd = wn * np.sqrt(1 - z**2)
                lines += [f"wd = {fmt(wd)} rad/s",
                          f"tr = {fmt((np.pi - np.arccos(z)) / wd)} s",
                          f"tp = {fmt(np.pi / wd)} s",
                          f"Mp = {fmt(100 * np.exp(-np.pi * z / np.sqrt(1 - z**2)))} %",
                          f"ts(2%) = {fmt(4 / (z * wn))} s"]
            else:
                lines.append("(Undamped/critically/over-damped: no overshoot formulas)" if z >= 1 else "(zeta<=0)")

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
        ttk.Button(self.side, text="Calculate & Plot", command=self.calculate).pack(pady=4)
        self.out.pack()

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
                lines.append(f"Equiv. 2nd-order zeta ~ {fmt(z)}")
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
        self.zt = self.field("Target zeta", "0.5", 2)
        self.button("Plot root locus", self.calculate, 3)
        self.out.grid(row=4, column=0, columnspan=2)

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
                lines.append(f"zeta={zt:g}: K = {fmt(K)}, dominant poles {fmt(c.real)} +/- j{fmt(c.imag)}")

        if self.K.get().strip():
            K = float(self.K.get())
            r = np.roots(np.polyadd(d, K * n))
            ax.plot(r.real, r.imag, "rs", mfc="none", ms=11, mew=2, label=f"poles at K={K:g}")
            lines += ["", f"--- Closed-loop poles at K = {K:g} ---"]
            for p in r:
                wn = abs(p)
                lines.append(f"  {np.round(p, 4)}   wn={fmt(wn)}  zeta={fmt(-p.real / wn)}")
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
        self.out.grid(row=2, column=0, columnspan=2)

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
        self.out.grid(row=2, column=0, columnspan=2)

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
        self.out.grid(row=8, column=0, columnspan=2)

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


# ---------------------------------------------------------------- main app
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Linear Control Systems Toolkit")
        self.geometry("1350x760")
        top = ttk.LabelFrame(self, text="Transfer function (coefficients, highest power first)")
        top.pack(fill="x", padx=8, pady=6)
        self.num, self.den = tk.StringVar(value="10"), tk.StringVar(value="1 2 0")
        self.unity = tk.BooleanVar(value=True)
        ttk.Label(top, text="Numerator:").pack(side="left", padx=4)
        ttk.Entry(top, textvariable=self.num, width=24).pack(side="left")
        ttk.Label(top, text="Denominator:").pack(side="left", padx=4)
        ttk.Entry(top, textvariable=self.den, width=24).pack(side="left")
        ttk.Checkbutton(top, text="Entered G(s) is open-loop with unity feedback",
                        variable=self.unity).pack(side="left", padx=10)
        ttk.Label(top, text="e.g. 10/(s^2+2s):  num = 10,  den = 1 2 0").pack(side="left", padx=6)
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=4)
        nb.add(TimeTab(nb, self.system, self.set_tf), text="Time Response")
        nb.add(FreqTab(nb, self.system), text="Frequency Response")
        nb.add(RootLocusTab(nb, self.system), text="Root Locus")
        nb.add(RouthTab(nb, self.system), text="Routh-Hurwitz")
        nb.add(NicholsTab(nb, self.system), text="Nichols Chart")
        nb.add(PIDTab(nb, self.system), text="PID Tuning")

    def set_tf(self, n, d, unity):
        self.num.set(n)
        self.den.set(d)
        self.unity.set(unity)

    def system(self):
        n, d = parse(self.num.get()), parse(self.den.get())
        if n.size > d.size:
            raise ValueError("Improper transfer function (numerator order > denominator order)")
        return n, d, self.unity.get()


if __name__ == "__main__":
    App().mainloop()