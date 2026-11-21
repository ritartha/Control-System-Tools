"""
Linear Control Systems Toolkit (Tkinter)
Tabs: Time Response (Ch. 6)  |  Frequency Response (Ch. 7)
Requires: numpy, scipy, matplotlib      ->  pip install numpy scipy matplotlib
Run:      python control_gui.py
"""
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from scipy import signal
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

    def calculate(self):
        try:
            self.run()
        except Exception as e:  # show friendly error
            messagebox.showerror("Error", str(e))


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
