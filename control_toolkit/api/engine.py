"""
Linear Control Systems Toolkit - Computation Engine
Ports all math from the original Tkinter app to return JSON-serializable dicts.
Each compute_* function returns {plots: {...}, results: {...}}.
Plot data uses Plotly.js trace format.
"""
import numpy as np
from scipy import signal
from scipy.optimize import linear_sum_assignment


# ──────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────

def _safe(v):
    """Convert numpy scalar / NaN / Inf to JSON-safe Python type."""
    if v is None:
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        if np.isnan(v):
            return None
        if np.isinf(v):
            return "Infinity" if v > 0 else "-Infinity"
        return float(v)
    if isinstance(v, (complex, np.complexfloating)):
        return {"real": _safe(v.real), "imag": _safe(v.imag)}
    return v


def _list(a):
    """numpy array → Python list of safe values."""
    if a is None:
        return []
    arr = np.asarray(a)
    result = []
    for x in arr.flat:
        if np.isnan(x) if np.isrealobj(arr) else False:
            result.append(None)
        elif isinstance(x, (complex, np.complexfloating)):
            result.append(float(x.real) if abs(x.imag) < 1e-15 else None)
        else:
            result.append(float(x))
    return result


def _clist(arr):
    """Complex array → list of {real, imag} dicts."""
    return [{"real": float(c.real), "imag": float(c.imag)} for c in np.asarray(arr)]


def fmt(v):
    """Format number to 4 significant figures."""
    if v is None:
        return "n/a"
    if isinstance(v, str):
        return v
    if isinstance(v, (float, np.floating)):
        if np.isnan(v):
            return "n/a"
        if np.isinf(v):
            return "∞" if v > 0 else "-∞"
    return f"{v:.4g}"


def parse(text):
    """Parse coefficient string: '1 2, 5' → array([1., 2., 5.])."""
    try:
        v = np.array([float(x) for x in text.replace(",", " ").split()])
    except ValueError:
        raise ValueError("Coefficients must be numbers, e.g. 1 2 5")
    v = np.trim_zeros(v, "f")
    if v.size == 0:
        raise ValueError("Polynomial cannot be empty / all zeros")
    return v


def closed_loop(n, d, unity):
    """Return (cl_num, cl_den)."""
    if unity:
        return (n, np.trim_zeros(np.polyadd(d, n), "f"))
    return (n, d)


def crossings(x, y, level):
    """Log-interpolated x where y crosses level."""
    s = y - level
    idx = np.flatnonzero(s[:-1] * s[1:] < 0)
    lx = np.log10(np.maximum(x, 1e-30))
    result = []
    for i in idx:
        denom = s[i + 1] - s[i]
        if abs(denom) < 1e-30:
            continue
        result.append(10 ** (lx[i] - s[i] * (lx[i + 1] - lx[i]) / denom))
    return result


def error_constants(n, d):
    """System type N and [Kp, Kv, Ka]."""
    N = 0
    dd = d.copy()
    while dd.size > 1 and abs(dd[-1]) < 1e-12:
        dd = dd[:-1]
        N += 1
    c = n[-1] / dd[-1] if abs(n[-1]) > 1e-12 else 0.0
    res = []
    for k in range(3):
        if N < k:
            res.append(0.0)
        elif N == k:
            res.append(c)
        else:
            res.append(float('inf'))
    return N, res


SUP = str.maketrans("0123456789", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")


def poly_str(c, var="s"):
    """[1, 2, 0] → 's² + 2s'."""
    c = np.trim_zeros(np.asarray(c, float), "f")
    if c.size == 0:
        return "0"
    deg, out = len(c) - 1, ""
    for i, a in enumerate(c):
        p = deg - i
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
    return "\n".join([
        " " * len(pre) + ns.center(w),
        pre + "─" * w,
        " " * len(pre) + ds.center(w)
    ])


def _stable(n, d, K=1.0):
    return bool(np.all(np.roots(np.polyadd(d, K * n)).real < -1e-9))


def marginal_gain(n, d, kmax=1e4):
    """(Ku, wu) where the locus crosses jw axis, or None."""
    Ks = np.logspace(-4, np.log10(kmax), 1500)
    def mr(K):
        return np.roots(np.polyadd(d, K * n)).real.max()
    v = np.array([mr(K) for K in Ks])
    idx = np.flatnonzero((v[:-1] < 0) & (v[1:] >= 0))
    if not idx.size:
        return None
    lo, hi = Ks[idx[0]], Ks[idx[0] + 1]
    for _ in range(60):
        mid = (lo + hi) / 2
        if mr(mid) < 0:
            lo = mid
        else:
            hi = mid
    r = np.roots(np.polyadd(d, hi * n))
    return hi, abs(r[np.argmax(r.real)].imag)


def stable_ranges(n, d, kmax=1e4):
    Ks = np.concatenate([[0], np.logspace(-4, np.log10(kmax), 2000)])
    st = np.array([_stable(n, d, K) for K in Ks])
    edges = []
    for i in np.flatnonzero(st[:-1] != st[1:]):
        lo, hi, s0 = Ks[i], Ks[i + 1], st[i]
        for _ in range(50):
            mid = (lo + hi) / 2
            if _stable(n, d, mid) == s0:
                lo = mid
            else:
                hi = mid
        edges.append(hi)
    out, start = [], (0.0 if st[0] else None)
    for e in edges:
        if start is None:
            start = e
        else:
            out.append((start, e))
            start = None
    if start is not None:
        out.append((start, float('inf')))
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
            notes.append(
                f"Row s^{n - i} was all zeros: auxiliary polynomial from row s^{deg} "
                "differentiated (roots symmetric about origin)."
            )
        if abs(R[i - 1, 0]) < 1e-12:
            R[i - 1, 0] = eps
            notes.append(f"Zero first element in row s^{n - i}: replaced by ε={eps:g}.")
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
    sys_lti = signal.lti(cn, cd)
    t, y = signal.step(sys_lti)
    r = {"t": t, "y": y, "stable": bool(np.all(poles.real < 0)), "poles": poles}
    if r["stable"]:
        ss = cn[-1] / cd[-1]
        r["ss"] = ss
        if abs(ss) > 1e-9:
            yn = y / ss
            first = lambda c: t[np.flatnonzero(c)[0]] if np.any(c) else np.nan
            ip = int(np.argmax(yn))
            out = np.flatnonzero(np.abs(yn - 1) > 0.02)
            r.update(
                tr=first(yn >= 0.9) - first(yn >= 0.1),
                tp=t[ip],
                Mp=max(0.0, (yn[ip] - 1) * 100),
                ts=0.0 if not out.size else (np.nan if out[-1] == len(t) - 1 else t[out[-1] + 1])
            )
    return r


def _spec_dict(r):
    """Convert step_specs result to JSON-safe dict."""
    d = {"stable": r["stable"]}
    if r["stable"]:
        d["final_value"] = _safe(r.get("ss"))
        d["ess"] = _safe(1 - r.get("ss", 1)) if r.get("ss") is not None else None
        d["tr"] = _safe(r.get("tr"))
        d["tp"] = _safe(r.get("tp"))
        d["Mp"] = _safe(r.get("Mp"))
        d["ts"] = _safe(r.get("ts"))
    return d


def margins_calc(w, G):
    mag = 20 * np.log10(np.abs(G) + 1e-30)
    ph = np.degrees(np.unwrap(np.angle(G)))
    lw = np.log10(w)
    gc = crossings(w, mag, 0)
    pc = crossings(w, ph, -180)
    PM = 180 + np.interp(np.log10(gc[0]), lw, ph) if gc else None
    GM = -np.interp(np.log10(pc[0]), lw, mag) if pc else None
    return mag, ph, (gc[0] if gc else None), PM, (pc[0] if pc else None), GM


# ──────────────────────────────────────────────────────────
# Compensator design helpers
# ──────────────────────────────────────────────────────────

def _pm_after(w, KG, nc, dc):
    _, _, _, pm, _, _ = margins_calc(w, KG * np.polyval(nc, 1j * w) / np.polyval(dc, 1j * w))
    return pm


def design_lead(w, KG, pm_des, margin):
    mag, ph, wgc, PM, _, _ = margins_calc(w, KG)
    if PM is None:
        raise ValueError("K*G(s) has no gain crossover - adjust K")
    if PM >= pm_des:
        raise ValueError(f"Phase margin is already {PM:.1f}° (≥ {pm_des:g}); no lead needed")
    for extra in np.arange(0, 40, 1.0):
        phim = np.radians(pm_des - PM + margin + extra)
        if phim >= np.radians(70):
            raise ValueError("Lead needs more than 70°: use lag-lead or two stages")
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
    return nc, dc, [
        f"α = {fmt(a)}   T = {fmt(T)} s",
        f"φ_max = {fmt(np.degrees(phim))}° at ωm = {fmt(wm)} rad/s",
        f"zero at s = -{fmt(1/T)},  pole at s = -{fmt(1/(a*T))}"
    ]


def design_lag(w, KG, pm_des, margin):
    mag, ph, wgc, PM, _, _ = margins_calc(w, KG)
    cr = crossings(w, ph, -180 + pm_des + margin)
    if not cr:
        raise ValueError("Phase never reaches the required value - lag compensator cannot give this PM")
    wc = cr[0]
    lw = np.log10(w)
    g = np.interp(np.log10(wc), lw, mag)
    if g <= 0:
        raise ValueError(f"Gain at the new crossover is already {g:.1f} dB - no lag needed")
    beta = 10 ** (g / 20)
    T = 10 / wc
    return np.array([T, 1.0]), np.array([beta * T, 1.0]), [
        f"new crossover ωc' = {fmt(wc)} rad/s  (|KG| = {fmt(g)} dB)",
        f"β = {fmt(beta)}   T = {fmt(T)} s",
        f"zero at s = -{fmt(1/T)},  pole at s = -{fmt(1/(beta*T))}"
    ]


def design_laglead(w, KG, pm_des, margin, w1=None):
    mag, ph = margins_calc(w, KG)[:2]
    lw = np.log10(w)
    if w1 is None:
        cr = crossings(w, ph, -180)
        if not cr:
            raise ValueError("No phase crossover found - enter w1 manually")
        w1 = cr[0]
    phi1 = np.interp(np.log10(w1), lw, ph)
    g1 = np.interp(np.log10(w1), lw, mag)
    lag_ph = 0.0
    for _ in range(40):
        phil = pm_des + margin - 180 - phi1 + abs(lag_ph)
        if not 0 < phil < 70:
            raise ValueError(f"Lead section would need {phil:.1f}° at w1={w1:.3g}: pick another w1")
        sn = np.sin(np.radians(phil))
        gam = (1 + sn) / (1 - sn)
        X = 101 * 10 ** ((g1 + 10 * np.log10(gam)) / 10)
        beta = np.sqrt(max((X - 1) / 100, 1.0))
        lag_ph = np.degrees(np.arctan(10) - np.arctan(10 * beta))
    T1 = np.sqrt(gam) / w1
    T2 = 10 / w1
    nc = np.polymul([T1, 1.0], [T2, 1.0])
    dc = np.polymul([T1 / gam, 1.0], [beta * T2, 1.0])
    return nc, dc, [
        f"new crossover ω1 = {fmt(w1)} rad/s  (|KG| = {fmt(g1)} dB, phase {fmt(phi1)}°)",
        f"Lead: γ = {fmt(gam)},  T1 = {fmt(T1)} s  (phase {fmt(phil)}° at ω1)",
        f"Lag:  β = {fmt(beta)},  T2 = {fmt(T2)} s  (phase {fmt(lag_ph)}° at ω1)"
    ]


# ──────────────────────────────────────────────────────────
# API: Transfer Function Display
# ──────────────────────────────────────────────────────────

def compute_transfer_function(num_str, den_str, unity):
    n, d = parse(num_str), parse(den_str)
    if n.size > d.size:
        raise ValueError("Improper TF (numerator order > denominator order)")

    poles = np.roots(d)
    zeros = np.roots(n)
    lines = []

    if unity:
        cn, cd = closed_loop(n, d, True)
        lines.append(tf_block("G(s)", n, d))
        lines.append("")
        lines.append("Closed loop  T(s) = G/(1+G):")
        lines.append(tf_block("T(s)", cn, cd))
        lines.append("")
        lines.append(f"Characteristic eq:  {poly_str(cd)} = 0")
        lines.append(f"Open-loop poles:   {np.round(poles, 4).tolist()}")
        cl_poles = np.roots(cd)
        lines.append(f"Closed-loop poles: {np.round(cl_poles, 4).tolist()}")
    else:
        lines.append(tf_block("T(s)", n, d))
        lines.append("")
        lines.append(f"Poles: {np.round(poles, 4).tolist()}")
        cl_poles = poles

    return {
        "display": "\n".join(lines),
        "poles": _clist(poles),
        "zeros": _clist(zeros),
        "closed_loop_poles": _clist(cl_poles),
    }


# ──────────────────────────────────────────────────────────
# API: Time Response
# ──────────────────────────────────────────────────────────

def compute_time_response(num_str, den_str, unity, signal_type="Step", t_final=None):
    n, d = parse(num_str), parse(den_str)
    if n.size > d.size:
        raise ValueError("Improper TF")
    cn, cd = closed_loop(n, d, unity)
    sys_lti = signal.lti(cn, cd)
    poles = np.roots(cd)
    zeros = np.roots(cn)
    is_stable = bool(np.all(poles.real < 0))

    T = np.linspace(0, t_final, 2000) if t_final else None
    t1, y1 = signal.step(sys_lti, T=T)

    if signal_type == "Step":
        t, y = t1, y1
    elif signal_type == "Impulse":
        t, y = signal.impulse(sys_lti, T=T)
    else:  # Ramp
        t = t1
        _, y, _ = signal.lsim(sys_lti, t, t)

    # ── Response plot ──
    response_traces = [
        {"x": _list(t), "y": _list(y), "mode": "lines", "name": "Output",
         "line": {"width": 2.5, "color": "#06b6d4"}}
    ]
    if signal_type == "Ramp":
        response_traces.append(
            {"x": _list(t), "y": _list(t), "mode": "lines", "name": "Input",
             "line": {"dash": "dash", "color": "#94a3b8"}}
        )
    elif signal_type == "Step":
        response_traces.append(
            {"x": [float(t[0]), float(t[-1])], "y": [1, 1], "mode": "lines", "name": "Input",
             "line": {"dash": "dash", "color": "#94a3b8", "width": 1}}
        )

    response_layout = {
        "title": f"{signal_type} Response",
        "xaxis": {"title": "Time (s)"},
        "yaxis": {"title": "c(t)"}
    }

    # ── Step specs ──
    step_specs_data = None
    if is_stable:
        ss = cn[-1] / cd[-1]
        step_specs_data = {"final_value": _safe(ss)}
        if abs(ss) > 1e-9:
            yn = y1 / ss
            first = lambda c: t1[np.flatnonzero(c)[0]] if np.any(c) else np.nan
            tr = first(yn >= 0.9) - first(yn >= 0.1)
            ip = int(np.argmax(yn))
            tp, Mp = t1[ip], max(0.0, (yn[ip] - 1) * 100)
            out = np.flatnonzero(np.abs(yn - 1) > 0.02)
            ts = 0.0 if out.size == 0 else (np.nan if out[-1] == len(t1) - 1 else t1[out[-1] + 1])
            step_specs_data.update({
                "rise_time": _safe(tr), "peak_time": _safe(tp),
                "overshoot": _safe(Mp), "settling_time": _safe(ts)
            })
            if signal_type == "Step":
                response_traces.append(
                    {"x": [float(tp)], "y": [float(yn[ip] * ss)], "mode": "markers",
                     "name": "Peak", "marker": {"symbol": "circle", "size": 10, "color": "#f87171"}}
                )
                if np.isfinite(ts):
                    response_layout.setdefault("shapes", []).append({
                        "type": "line", "x0": float(ts), "x1": float(ts),
                        "y0": 0, "y1": float(np.max(y) * 1.1),
                        "line": {"color": "#f87171", "dash": "dot", "width": 1}
                    })

    # ── 2nd order formulas ──
    second_order = None
    if cn.size == 1 and cd.size == 3:
        a = cd / cd[0]
        wn = np.sqrt(a[2])
        z = a[1] / (2 * wn)
        second_order = {"wn": _safe(wn), "zeta": _safe(z)}
        if 0 < z < 1:
            wd = wn * np.sqrt(1 - z ** 2)
            second_order.update({
                "wd": _safe(wd),
                "tr": _safe((np.pi - np.arccos(z)) / wd),
                "tp": _safe(np.pi / wd),
                "Mp": _safe(100 * np.exp(-np.pi * z / np.sqrt(1 - z ** 2))),
                "ts": _safe(4 / (z * wn)),
                "type_desc": "Underdamped"
            })
        elif z >= 1:
            second_order["type_desc"] = "Critically damped" if abs(z - 1) < 1e-6 else "Overdamped"
        else:
            second_order["type_desc"] = "Unstable (ζ < 0)"

    # ── Steady-state error ──
    steady_state = None
    if unity:
        N, (Kp, Kv, Ka) = error_constants(n, d)
        ess_step = 1 / (1 + Kp) if np.isfinite(Kp) else 0.0
        ess_ramp = float('inf') if Kv == 0 else 1 / Kv
        ess_par = float('inf') if Ka == 0 else 1 / Ka
        steady_state = {
            "type": int(N), "Kp": _safe(Kp), "Kv": _safe(Kv), "Ka": _safe(Ka),
            "ess_step": _safe(ess_step), "ess_ramp": _safe(ess_ramp), "ess_parabolic": _safe(ess_par)
        }

    # ── Pole-Zero map plot ──
    all_pz = np.concatenate([poles, zeros]) if zeros.size else poles
    pz_ext = max(1.0, np.abs(all_pz).max() * 1.8)
    pz_traces = [
        {"x": _list(poles.real), "y": _list(poles.imag), "mode": "markers",
         "name": "Poles", "marker": {"symbol": "x", "size": 12, "color": "#3b82f6", "line": {"width": 2}}}
    ]
    if zeros.size:
        pz_traces.append(
            {"x": _list(zeros.real), "y": _list(zeros.imag), "mode": "markers",
             "name": "Zeros", "marker": {"symbol": "circle-open", "size": 10, "color": "#4ade80", "line": {"width": 2}}}
        )
    pz_layout = {
        "title": "Pole-Zero Map",
        "xaxis": {"title": "Real", "zeroline": True, "range": [float(-pz_ext), float(pz_ext)]},
        "yaxis": {"title": "Imag", "zeroline": True, "scaleanchor": "x", "range": [float(-pz_ext), float(pz_ext)]},
        "shapes": [
            # LHS (stable) green shading
            {"type": "rect", "x0": float(-pz_ext), "x1": 0, "y0": float(-pz_ext), "y1": float(pz_ext),
             "fillcolor": "rgba(74,222,128,0.07)", "line": {"width": 0}, "layer": "below"},
            # RHS (unstable) red shading
            {"type": "rect", "x0": 0, "x1": float(pz_ext), "y0": float(-pz_ext), "y1": float(pz_ext),
             "fillcolor": "rgba(248,113,113,0.07)", "line": {"width": 0}, "layer": "below"},
            # jω axis
            {"type": "line", "x0": 0, "x1": 0, "y0": float(-pz_ext), "y1": float(pz_ext),
             "line": {"color": "rgba(100,116,139,0.5)", "width": 1.5, "dash": "dot"}}
        ]
    }

    return {
        "plots": {
            "response": {"traces": response_traces, "layout": response_layout},
            "pole_zero": {"traces": pz_traces, "layout": pz_layout}
        },
        "results": {
            "poles": _clist(poles), "zeros": _clist(zeros),
            "stable": is_stable,
            "step_specs": step_specs_data,
            "second_order": second_order,
            "steady_state": steady_state
        }
    }


# ──────────────────────────────────────────────────────────
# API: Frequency Response
# ──────────────────────────────────────────────────────────

def compute_frequency_response(num_str, den_str, unity):
    n, d = parse(num_str), parse(den_str)
    if n.size > d.size:
        raise ValueError("Improper TF")

    w = freq_grid(n, d)
    G = np.polyval(n, 1j * w) / np.polyval(d, 1j * w)
    mag = 20 * np.log10(np.abs(G) + 1e-30)
    ph = np.degrees(np.unwrap(np.angle(G)))

    wgc_list = crossings(w, mag, 0)
    wpc_list = crossings(w, ph, -180)
    lw = np.log10(w)
    PM = 180 + np.interp(np.log10(wgc_list[0]), lw, ph) if wgc_list else None
    GM = -np.interp(np.log10(wpc_list[0]), lw, mag) if wpc_list else None

    # ── Bode magnitude ──
    bode_mag_traces = [
        {"x": _list(w), "y": _list(mag), "mode": "lines", "name": "|G(jω)|",
         "line": {"width": 2, "color": "#06b6d4"}}
    ]
    bode_mag_layout = {
        "title": "Bode Magnitude",
        "xaxis": {"title": "ω (rad/s)", "type": "log"},
        "yaxis": {"title": "Magnitude (dB)"},
        "shapes": [{"type": "line", "x0": float(w[0]), "x1": float(w[-1]),
                     "y0": 0, "y1": 0, "line": {"color": "#64748b", "width": 1}}]
    }
    if wgc_list:
        bode_mag_traces.append(
            {"x": [wgc_list[0], wgc_list[0]], "y": [float(mag.min()), float(mag.max())],
             "mode": "lines", "name": f"ωgc={fmt(wgc_list[0])}",
             "line": {"dash": "dash", "color": "#f87171", "width": 1}}
        )
    if wpc_list:
        bode_mag_traces.append(
            {"x": [wpc_list[0], wpc_list[0]], "y": [float(mag.min()), float(mag.max())],
             "mode": "lines", "name": f"ωpc={fmt(wpc_list[0])}",
             "line": {"dash": "dash", "color": "#4ade80", "width": 1}}
        )

    # ── Bode phase ──
    bode_phase_traces = [
        {"x": _list(w), "y": _list(ph), "mode": "lines", "name": "∠G(jω)",
         "line": {"width": 2, "color": "#a78bfa"}}
    ]
    bode_phase_layout = {
        "title": "Bode Phase",
        "xaxis": {"title": "ω (rad/s)", "type": "log"},
        "yaxis": {"title": "Phase (deg)"},
        "shapes": [{"type": "line", "x0": float(w[0]), "x1": float(w[-1]),
                     "y0": -180, "y1": -180, "line": {"color": "#64748b", "width": 1}}]
    }
    if wgc_list:
        bode_phase_traces.append(
            {"x": [wgc_list[0], wgc_list[0]], "y": [float(ph.min()), float(ph.max())],
             "mode": "lines", "name": "ωgc", "showlegend": False,
             "line": {"dash": "dash", "color": "#f87171", "width": 1}}
        )
    if wpc_list:
        bode_phase_traces.append(
            {"x": [wpc_list[0], wpc_list[0]], "y": [float(ph.min()), float(ph.max())],
             "mode": "lines", "name": "ωpc", "showlegend": False,
             "line": {"dash": "dash", "color": "#4ade80", "width": 1}}
        )

    # ── Nyquist ──
    Gp = np.where(np.abs(G) < 20, G, np.nan + 0j)
    nyq_traces = [
        {"x": _list(np.real(Gp)), "y": _list(np.imag(Gp)),
         "mode": "lines", "name": "ω > 0", "line": {"width": 2, "color": "#3b82f6"}},
        {"x": _list(np.real(Gp)), "y": _list(-np.imag(Gp)),
         "mode": "lines", "name": "ω < 0", "line": {"dash": "dash", "width": 1.5, "color": "#3b82f6"}},
        {"x": [-1], "y": [0], "mode": "markers", "name": "-1+j0",
         "marker": {"symbol": "x", "size": 12, "color": "#f87171", "line": {"width": 2}}}
    ]
    nyq_layout = {
        "title": "Nyquist Plot",
        "xaxis": {"title": "Real", "zeroline": True},
        "yaxis": {"title": "Imag", "zeroline": True, "scaleanchor": "x"}
    }

    # ── Closed-loop info ──
    cn, cd = closed_loop(n, d, unity)
    cl_poles = np.roots(cd)
    cl_stable = bool(np.all(cl_poles.real < 0))

    T = G / (1 + G) if unity else G
    Td = 20 * np.log10(np.abs(T) + 1e-30)
    dc_gain = Td[0]
    i_peak = int(np.argmax(Td))
    low = np.flatnonzero(Td < dc_gain - 3)
    bw = float(w[low[0]]) if low.size else None

    Mr = None
    wr = None
    equiv_zeta = None
    if Td[i_peak] - dc_gain > 0.01:
        Mr_val = 10 ** ((Td[i_peak] - dc_gain) / 20)
        Mr = _safe(Mr_val)
        wr = _safe(w[i_peak])
        if Mr_val > 1:
            eq_z = np.sqrt((1 - np.sqrt(1 - 1 / Mr_val ** 2)) / 2)
            equiv_zeta = _safe(eq_z)

    margins_verdict = None
    if wgc_list and wpc_list and GM is not None and PM is not None:
        margins_verdict = "stable" if GM > 0 and PM > 0 else "unstable"

    return {
        "plots": {
            "bode_mag": {"traces": bode_mag_traces, "layout": bode_mag_layout},
            "bode_phase": {"traces": bode_phase_traces, "layout": bode_phase_layout},
            "nyquist": {"traces": nyq_traces, "layout": nyq_layout}
        },
        "results": {
            "wgc": _safe(wgc_list[0]) if wgc_list else None,
            "PM": _safe(PM), "wpc": _safe(wpc_list[0]) if wpc_list else None,
            "GM": _safe(GM),
            "closed_loop_stable": cl_stable,
            "bandwidth": _safe(bw), "dc_gain": _safe(dc_gain),
            "Mr": Mr, "wr": wr, "equiv_zeta": equiv_zeta,
            "margins_verdict": margins_verdict
        }
    }


# ──────────────────────────────────────────────────────────
# API: Root Locus
# ──────────────────────────────────────────────────────────

def compute_root_locus(num_str, den_str, k_max=None, k_eval=None, zeta_target=None):
    n, d = parse(num_str), parse(den_str)
    P = len(d) - 1
    Z = len(n) - 1
    if P < 1:
        raise ValueError("Denominator must have order >= 1")

    poles = np.roots(d)
    zeros = np.roots(n)
    mg = marginal_gain(n, d)
    ext = max(1.0, np.abs(np.concatenate([poles, zeros])).max()) * 2.5

    if k_max is not None and k_max > 0:
        km = k_max
    elif mg:
        km = mg[0] * 4
    else:
        km = 10 * abs(np.polyval(d, 2j * ext) / (np.polyval(n, 2j * ext) + 1e-30))

    Ks = np.concatenate([[0], np.logspace(np.log10(km) - 5, np.log10(km), 1500)])
    R = np.zeros((len(Ks), P), complex)
    R[0] = prev = poles
    for i in range(1, len(Ks)):
        r = np.roots(np.polyadd(d, Ks[i] * n))
        if len(r) != P:
            r = np.pad(r, (0, P - len(r)), constant_values=0)[:P]
        ri, ci = linear_sum_assignment(np.abs(prev[:, None] - r[None, :]))
        R[i, ri] = r[ci]
        prev = R[i]

    # ── Plot traces ──
    colors = ["#06b6d4", "#a78bfa", "#f97316", "#4ade80", "#f87171",
              "#fbbf24", "#ec4899", "#14b8a6", "#8b5cf6", "#64748b"]
    traces = []
    for b in range(P):
        traces.append({
            "x": _list(R[:, b].real), "y": _list(R[:, b].imag),
            "mode": "lines", "name": f"Branch {b+1}",
            "line": {"width": 1.5, "color": colors[b % len(colors)]},
            "showlegend": False
        })
    traces.append({
        "x": _list(poles.real), "y": _list(poles.imag),
        "mode": "markers", "name": "Open-loop poles",
        "marker": {"symbol": "x", "size": 12, "color": "#dc2626", "line": {"width": 2}}
    })
    if zeros.size:
        traces.append({
            "x": _list(zeros.real), "y": _list(zeros.imag),
            "mode": "markers", "name": "Zeros",
            "marker": {"symbol": "circle-open", "size": 10, "color": "#16a34a", "line": {"width": 2}}
        })

    # ── Results ──
    results = {
        "poles": _clist(poles), "zeros": _clist(zeros) if zeros.size else [],
        "branches": P, "num_asymptotes": max(P - Z, 0)
    }

    if P > Z:
        ang = (2 * np.arange(P - Z) + 1) * 180 / (P - Z)
        cen = (poles.real.sum() - zeros.real.sum()) / (P - Z)
        results["asymptote_angles"] = _list(np.round(ang, 2))
        results["centroid"] = _safe(cen)
        # Draw asymptotes
        for a_deg in ang:
            a_rad = np.radians(a_deg)
            traces.append({
                "x": [float(cen), float(cen + 3 * ext * np.cos(a_rad))],
                "y": [0, float(3 * ext * np.sin(a_rad))],
                "mode": "lines", "showlegend": False,
                "line": {"dash": "dot", "width": 0.8, "color": "#64748b"}
            })
    else:
        results["asymptote_angles"] = []
        results["centroid"] = None

    # ── Breakpoints ──
    dp = np.trim_zeros(
        np.polysub(np.polymul(np.polyder(d), n), np.polymul(d, np.polyder(n))), "f"
    )
    bk = []
    if dp.size > 1:
        for s in np.roots(dp):
            if abs(s.imag) < 1e-7:
                nv = np.polyval(n, s.real)
                if abs(nv) > 1e-15:
                    K_val = -np.polyval(d, s.real) / nv
                    if K_val > 0:
                        bk.append({"s": _safe(s.real), "K": _safe(float(K_val))})
                        traces.append({
                            "x": [float(s.real)], "y": [0],
                            "mode": "markers", "name": f"Break s={fmt(s.real)}",
                            "marker": {"symbol": "diamond", "size": 9, "color": "#4ade80"},
                            "showlegend": False
                        })
    results["breakpoints"] = bk

    # ── jω crossing ──
    if mg:
        results["jw_crossing"] = {"K": _safe(mg[0]), "w": _safe(mg[1])}
        traces.append({
            "x": [0, 0], "y": [float(mg[1]), float(-mg[1])],
            "mode": "markers", "name": f"jω crossing K={fmt(mg[0])}",
            "marker": {"symbol": "star", "size": 12, "color": "#f87171"}
        })
        results["stable_range"] = f"0 < K < {fmt(mg[0])}"
    else:
        results["jw_crossing"] = None
        results["stable_range"] = "Stable for all K tested"

    # ── Zeta target ──
    zeta_points = []
    if zeta_target is not None and 0 < zeta_target < 1:
        a_ang = np.arccos(zeta_target)
        for sg in (1, -1):
            traces.append({
                "x": [0, float(-ext * 3 * np.cos(a_ang))],
                "y": [0, float(sg * ext * 3 * np.sin(a_ang))],
                "mode": "lines", "showlegend": False,
                "line": {"dash": "dash", "width": 1, "color": "#a78bfa"}
            })
        with np.errstate(all="ignore"):
            z_arr = np.where(R.imag > 1e-6, -R.real / np.abs(R), np.inf)
        zm = np.minimum(z_arr.min(axis=1), 1.0)
        for K_val in crossings(Ks[1:], zm[1:], zeta_target)[:3]:
            r = np.roots(np.polyadd(d, K_val * n))
            c = r[np.argmax(r.imag)]
            zeta_points.append({"K": _safe(K_val), "real": _safe(c.real), "imag": _safe(c.imag)})
            traces.append({
                "x": [float(c.real), float(c.real)],
                "y": [float(c.imag), float(-c.imag)],
                "mode": "markers", "showlegend": False,
                "marker": {"symbol": "square", "size": 10, "color": "#a78bfa"}
            })
    results["zeta_points"] = zeta_points

    # ── Evaluate at K ──
    if k_eval is not None:
        r = np.roots(np.polyadd(d, k_eval * n))
        kep = []
        for p in r:
            wn = abs(p)
            kep.append({
                "value": f"{np.round(p, 4)}", "wn": _safe(wn),
                "zeta": _safe(-p.real / wn) if wn > 1e-12 else None
            })
        results["k_eval_poles"] = kep
        results["k_eval_stable"] = bool(np.all(r.real < 0))
        traces.append({
            "x": _list(r.real), "y": _list(r.imag),
            "mode": "markers", "name": f"Poles at K={k_eval:g}",
            "marker": {"symbol": "square-open", "size": 13, "color": "#fbbf24", "line": {"width": 2}}
        })
    else:
        results["k_eval_poles"] = None

    layout = {
        "title": "Root Locus (K from 0 upward)",
        "xaxis": {"title": "Real", "range": [float(-ext), float(ext)], "zeroline": True},
        "yaxis": {"title": "Imag", "range": [float(-ext), float(ext)], "zeroline": True, "scaleanchor": "x"},
        "shapes": [
            # LHS (stable) green shading
            {"type": "rect", "x0": float(-ext), "x1": 0, "y0": float(-ext), "y1": float(ext),
             "fillcolor": "rgba(74,222,128,0.07)", "line": {"width": 0}, "layer": "below"},
            # RHS (unstable) red shading
            {"type": "rect", "x0": 0, "x1": float(ext), "y0": float(-ext), "y1": float(ext),
             "fillcolor": "rgba(248,113,113,0.07)", "line": {"width": 0}, "layer": "below"},
            # jω axis emphasis
            {"type": "line", "x0": 0, "x1": 0, "y0": float(-ext), "y1": float(ext),
             "line": {"color": "rgba(100,116,139,0.5)", "width": 1.5, "dash": "dot"}}
        ]
    }

    return {
        "plots": {"root_locus": {"traces": traces, "layout": layout}},
        "results": results
    }


# ──────────────────────────────────────────────────────────
# API: Routh-Hurwitz
# ──────────────────────────────────────────────────────────

def compute_routh_hurwitz(num_str, den_str, unity, K=1.0):
    n, d = parse(num_str), parse(den_str)
    if n.size > d.size:
        raise ValueError("Improper TF")

    p = np.trim_zeros(np.polyadd(d, K * n) if unity else d, "f")
    R, notes = routh(p)
    deg = len(p) - 1

    eq_desc = ("den + K·num = 0" if unity else "den = 0") + (f"  (K={K:g})" if unity else "")

    row_labels = [f"s^{deg - i}" for i in range(len(R))]
    routh_table = []
    for row in R:
        routh_table.append([_safe(x) for x in row])

    roots_arr = np.roots(p)
    rhp = int(np.sum(roots_arr.real > 1e-9))
    jw = int(np.sum(np.abs(roots_arr.real) <= 1e-9))
    sc = int(np.sum(np.sign(R[:-1, 0]) * np.sign(R[1:, 0]) < 0))

    if rhp == 0 and jw == 0:
        verdict = "STABLE"
    elif rhp == 0:
        verdict = "MARGINALLY STABLE"
    else:
        verdict = "UNSTABLE"

    sr = []
    if unity:
        sr_raw = stable_ranges(n, d)
        for a, b in sr_raw:
            sr.append({"lo": _safe(a), "hi": _safe(b)})

    # ── Root plot ──
    traces = [
        {"x": _list(roots_arr.real), "y": _list(roots_arr.imag),
         "mode": "markers", "name": "Roots",
         "marker": {"symbol": "x", "size": 12, "color": "#3b82f6", "line": {"width": 2}}}
    ]
    max_r = float(max(1, np.abs(roots_arr).max() * 1.5)) if roots_arr.size else 2.0
    layout = {
        "title": "Roots of Characteristic Equation",
        "xaxis": {"title": "Real", "zeroline": True, "range": [-max_r, max_r]},
        "yaxis": {"title": "Imag", "zeroline": True, "scaleanchor": "x", "range": [-max_r, max_r]},
        "shapes": [
            # LHS (stable) green shading
            {"type": "rect", "x0": float(-max_r), "x1": 0, "y0": float(-max_r), "y1": float(max_r),
             "fillcolor": "rgba(74,222,128,0.07)", "line": {"width": 0}, "layer": "below"},
            # RHS (unstable) red shading
            {"type": "rect", "x0": 0, "x1": float(max_r), "y0": float(-max_r), "y1": float(max_r),
             "fillcolor": "rgba(248,113,113,0.07)", "line": {"width": 0}, "layer": "below"},
            # jω axis
            {"type": "line", "x0": 0, "x1": 0, "y0": float(-max_r), "y1": float(max_r),
             "line": {"color": "rgba(100,116,139,0.5)", "width": 1.5, "dash": "dot"}}
        ]
    }

    return {
        "plots": {"roots": {"traces": traces, "layout": layout}},
        "results": {
            "equation_desc": eq_desc,
            "coefficients": _list(np.round(p, 5)),
            "routh_table": routh_table, "row_labels": row_labels,
            "sign_changes": sc, "rhp_roots": rhp, "jw_roots": jw,
            "verdict": verdict, "notes": notes,
            "stable_ranges": sr,
            "actual_roots": _clist(roots_arr)
        }
    }


# ──────────────────────────────────────────────────────────
# API: Nichols Chart
# ──────────────────────────────────────────────────────────

def compute_nichols(num_str, den_str, target_mr=None):
    n, d = parse(num_str), parse(den_str)
    if n.size > d.size:
        raise ValueError("Improper TF")

    w = freq_grid(n, d)
    G = np.polyval(n, 1j * w) / np.polyval(d, 1j * w)
    mag = 20 * np.log10(np.abs(G) + 1e-30)
    ph = np.degrees(np.unwrap(np.angle(G)))

    # Wrapped phase for Nichols
    xp = (ph % 360) - 360
    j_breaks = np.flatnonzero(np.abs(np.diff(xp)) > 180) + 1
    xp_plot = np.insert(xp, j_breaks, np.nan)
    yp_plot = np.insert(mag, j_breaks, np.nan)

    traces = []

    # ── M-contours and N-contours (simplified as reference grid) ──
    thd = np.linspace(-359.5, -0.5, 720)
    th = np.radians(thd)
    for Mdb in (-20, -12, -6, -3, -1, 0, 1, 3, 6, 12):
        M2 = 10 ** (Mdb / 10)
        with np.errstate(invalid="ignore", divide="ignore"):
            if Mdb == 0:
                r_vals = np.where(np.cos(th) < 0, -1 / (2 * np.cos(th)), np.nan)
                y_vals = 20 * np.log10(np.where(r_vals > 0, r_vals, np.nan))
                traces.append({
                    "x": _list(thd), "y": _list(y_vals),
                    "mode": "lines", "showlegend": False,
                    "line": {"width": 0.7, "color": "rgba(248,113,113,0.5)"},
                    "hoverinfo": "skip"
                })
            else:
                a_coef = M2 - 1
                b_coef = 2 * M2 * np.cos(th)
                disc = b_coef ** 2 - 4 * a_coef * M2
                for sign in (1, -1):
                    r_vals = np.where(disc >= 0, (-b_coef + sign * np.sqrt(np.maximum(disc, 0))) / (2 * a_coef), np.nan)
                    r_vals = np.where(r_vals > 0, r_vals, np.nan)
                    y_vals = 20 * np.log10(np.where(r_vals > 0, r_vals, np.nan))
                    traces.append({
                        "x": _list(thd), "y": _list(y_vals),
                        "mode": "lines", "showlegend": False,
                        "line": {"width": 0.6, "color": "rgba(248,113,113,0.4)"},
                        "hoverinfo": "skip"
                    })

    # N-contours
    for al in (-5, -15, -30, -45, -60, -90, -120, -150):
        a_rad = np.radians(al)
        r_vals = np.sin(th) * np.cos(a_rad) / np.sin(a_rad) - np.cos(th)
        with np.errstate(invalid="ignore", divide="ignore"):
            y_vals = 20 * np.log10(np.where(r_vals > 0, r_vals, np.nan))
        traces.append({
            "x": _list(thd), "y": _list(y_vals),
            "mode": "lines", "showlegend": False,
            "line": {"width": 0.6, "color": "rgba(74,222,128,0.4)"},
            "hoverinfo": "skip"
        })

    # ── Main curve ──
    traces.append({
        "x": _list(xp_plot), "y": _list(yp_plot),
        "mode": "lines", "name": "L(jω)",
        "line": {"width": 2.5, "color": "#06b6d4"}
    })
    traces.append({
        "x": [-180], "y": [0], "mode": "markers", "name": "Critical point",
        "marker": {"symbol": "x", "size": 12, "color": "#f87171", "line": {"width": 2}}
    })

    layout = {
        "title": "Nichols Chart",
        "xaxis": {"title": "Phase (deg)", "range": [-360, 0]},
        "yaxis": {
            "title": "Gain (dB)",
            "range": [min(-40, max(-100, float(mag.min()))), max(20, min(60, float(mag.max())))]
        }
    }

    # ── Margins ──
    wgc_list = crossings(w, mag, 0)
    wpc_list = crossings(w, ph, -180)
    lw = np.log10(w)
    PM = 180 + np.interp(np.log10(wgc_list[0]), lw, ph) if wgc_list else None
    GM_val = -np.interp(np.log10(wpc_list[0]), lw, mag) if wpc_list else None
    is_stable = _stable(n, d)

    Td = 20 * np.log10(np.abs(G / (1 + G)) + 1e-30)
    i_peak = int(np.argmax(Td))
    low = np.flatnonzero(Td < Td[0] - 3)
    bw = float(w[low[0]]) if low.size else None

    if i_peak > 0 and Td[i_peak] > Td[0] + 0.01:
        peak_info = f"Mr = {fmt(10 ** (Td[i_peak] / 20))} ({fmt(Td[i_peak])} dB) at ωr = {fmt(w[i_peak])} rad/s"
    else:
        peak_info = "No resonant peak (Mr = 1)"

    # ── Gain adjustment for target Mr ──
    gain_adjustments = []
    if target_mr is not None:
        Ks_scan = np.logspace(-3, 3, 500)
        with np.errstate(all="ignore"):
            mr_scan = np.array([(20 * np.log10(np.abs(K * G / (1 + K * G)) + 1e-30)).max() for K in Ks_scan])
        cr = crossings(Ks_scan, mr_scan, target_mr)
        cr = [k for k in cr if _stable(n, d, k)]
        if cr:
            gain_adjustments = [f"Multiply loop gain by K = {fmt(k)}" for k in cr[:3]]
        else:
            gain_adjustments = ["No stable gain in 0.001..1000 gives that Mr"]

    return {
        "plots": {"nichols": {"traces": traces, "layout": layout}},
        "results": {
            "GM": _safe(GM_val), "PM": _safe(PM), "stable": is_stable,
            "peak_info": peak_info, "bandwidth": _safe(bw),
            "gain_adjustments": gain_adjustments
        }
    }


# ──────────────────────────────────────────────────────────
# API: PID Tuning
# ──────────────────────────────────────────────────────────

def compute_pid(num_str, den_str, Kp=1.0, Ki=0.0, Kd=0.0, N=20.0, zn_type=None):
    n, d = parse(num_str), parse(den_str)
    if n.size > d.size:
        raise ValueError("Improper TF")

    zn_info = None

    # ── Ziegler-Nichols auto-tune ──
    if zn_type:
        mg = marginal_gain(n, d)
        if not mg:
            raise ValueError("Plant never reaches marginal stability with P control — Z-N not applicable")
        Ku, wu = mg
        Pu = 2 * np.pi / wu
        if zn_type == "P":
            Kp, ti, td = 0.5 * Ku, float('inf'), 0
        elif zn_type == "PI":
            Kp, ti, td = 0.45 * Ku, Pu / 1.2, 0
        else:  # PID
            Kp, ti, td = 0.6 * Ku, Pu / 2, Pu / 8
        Ki = Kp / ti if np.isfinite(ti) else 0.0
        Kd = Kp * td
        zn_info = f"Ziegler-Nichols {zn_type}: Ku={fmt(Ku)}, Pu={fmt(Pu)} s"

    # ── Build PID controller TF ──
    tau = Kd / N if Kd > 0 else 0.0
    f = np.trim_zeros(np.array([tau, 1.0]), "f")
    if Ki != 0:
        dc = np.polymul(f, [1.0, 0.0])
        nc = np.polyadd(np.polyadd(Kp * dc, Ki * f), Kd * np.array([1.0, 0, 0]))
    else:
        dc = f
        nc = np.polyadd(Kp * f, Kd * np.array([1.0, 0]))
    nc = np.trim_zeros(nc, "f")
    if nc.size == 0:
        raise ValueError("All controller gains are zero")

    Ln = np.polymul(nc, n)
    Ld = np.polymul(dc, d)

    ru = step_specs(n, np.trim_zeros(np.polyadd(d, n), "f"))
    rp = step_specs(Ln, np.trim_zeros(np.polyadd(Ld, Ln), "f"))

    traces = [
        {"x": _list(ru["t"]), "y": _list(ru["y"]),
         "mode": "lines", "name": "P-only, K=1 (plant + unity fb)",
         "line": {"dash": "dash", "width": 1.5, "color": "#94a3b8"}},
        {"x": _list(rp["t"]), "y": _list(rp["y"]),
         "mode": "lines", "name": f"PID  Kp={Kp:g} Ki={Ki:g} Kd={Kd:g}",
         "line": {"width": 2.5, "color": "#06b6d4"}},
    ]
    layout = {
        "title": "Closed-Loop Step Response",
        "xaxis": {"title": "Time (s)"},
        "yaxis": {"title": "y(t)"},
        "shapes": [{"type": "line", "x0": 0, "x1": float(max(rp["t"][-1], ru["t"][-1])),
                     "y0": 1, "y1": 1, "line": {"color": "#64748b", "width": 1}}]
    }
    if not (ru["stable"] and rp["stable"]):
        layout["yaxis"]["range"] = [-1, 3]

    with_pid_data = _spec_dict(rp)
    with_pid_data["poles"] = _clist(rp["poles"])

    return {
        "plots": {"step_response": {"traces": traces, "layout": layout}},
        "results": {
            "zn_info": zn_info,
            "Kp": _safe(Kp), "Ki": _safe(Ki), "Kd": _safe(Kd),
            "without_pid": _spec_dict(ru),
            "with_pid": with_pid_data
        }
    }


# ──────────────────────────────────────────────────────────
# API: Compensator Design
# ──────────────────────────────────────────────────────────

def compute_compensator(num_str, den_str, comp_type="Lead", K=1.0, pm=45.0, margin=5.0, w1=None):
    n, d = parse(num_str), parse(den_str)
    if n.size > d.size:
        raise ValueError("Improper TF")

    w = freq_grid(K * n, d, 6000)
    KG = K * np.polyval(n, 1j * w) / np.polyval(d, 1j * w)

    if comp_type == "Lead":
        nc, dc, info = design_lead(w, KG, pm, margin)
    elif comp_type == "Lag":
        nc, dc, info = design_lag(w, KG, pm, margin)
    else:
        nc, dc, info = design_laglead(w, KG, pm, margin, w1)

    Gc = np.polyval(nc, 1j * w) / np.polyval(dc, 1j * w)
    m0, p0, wg0, PM0, wp0, GM0 = margins_calc(w, KG)
    m1, p1, wg1, PM1, wp1, GM1 = margins_calc(w, KG * Gc)

    cn0 = K * n
    cd0 = np.trim_zeros(np.polyadd(d, K * n), "f")
    Ln = K * np.polymul(n, nc)
    Ld = np.polymul(d, dc)
    cd1 = np.trim_zeros(np.polyadd(Ld, Ln), "f")
    r0 = step_specs(cn0, cd0)
    r1 = step_specs(Ln, cd1)

    Gc_mag = 20 * np.log10(np.abs(Gc) + 1e-30)
    Gc_ph = np.degrees(np.unwrap(np.angle(Gc)))

    # ── Bode magnitude ──
    bode_mag_traces = [
        {"x": _list(w), "y": _list(m0), "mode": "lines", "name": "Uncompensated",
         "line": {"dash": "dash", "width": 1.5, "color": "#94a3b8"}},
        {"x": _list(w), "y": _list(m1), "mode": "lines", "name": "Compensated",
         "line": {"width": 2.5, "color": "#06b6d4"}},
        {"x": _list(w), "y": _list(Gc_mag), "mode": "lines", "name": "Gc(jω)",
         "line": {"dash": "dot", "width": 1.5, "color": "#4ade80"}}
    ]
    bode_mag_layout = {
        "title": "Bode Magnitude",
        "xaxis": {"title": "ω (rad/s)", "type": "log"},
        "yaxis": {"title": "dB"},
        "shapes": [{"type": "line", "x0": float(w[0]), "x1": float(w[-1]),
                     "y0": 0, "y1": 0, "line": {"color": "#64748b", "width": 1}}]
    }

    # ── Bode phase ──
    bode_phase_traces = [
        {"x": _list(w), "y": _list(p0), "mode": "lines", "name": "Uncompensated",
         "line": {"dash": "dash", "width": 1.5, "color": "#94a3b8"}},
        {"x": _list(w), "y": _list(p1), "mode": "lines", "name": "Compensated",
         "line": {"width": 2.5, "color": "#a78bfa"}},
        {"x": _list(w), "y": _list(Gc_ph), "mode": "lines", "name": "Gc(jω)",
         "line": {"dash": "dot", "width": 1.5, "color": "#4ade80"}}
    ]
    bode_phase_layout = {
        "title": "Bode Phase",
        "xaxis": {"title": "ω (rad/s)", "type": "log"},
        "yaxis": {"title": "deg"},
        "shapes": [{"type": "line", "x0": float(w[0]), "x1": float(w[-1]),
                     "y0": -180, "y1": -180, "line": {"color": "#64748b", "width": 1}}]
    }

    # ── Step response ──
    step_traces = [
        {"x": _list(r0["t"]), "y": _list(r0["y"]), "mode": "lines", "name": "Uncompensated",
         "line": {"dash": "dash", "width": 1.5, "color": "#94a3b8"}},
        {"x": _list(r1["t"]), "y": _list(r1["y"]), "mode": "lines", "name": "Compensated",
         "line": {"width": 2.5, "color": "#06b6d4"}},
    ]
    step_layout = {
        "title": "Closed-Loop Step Response",
        "xaxis": {"title": "Time (s)"},
        "yaxis": {"title": "y(t)"},
        "shapes": [{"type": "line", "x0": 0, "x1": float(max(r0["t"][-1], r1["t"][-1])),
                     "y0": 1, "y1": 1, "line": {"color": "#64748b", "width": 1}}]
    }
    if not (r0["stable"] and r1["stable"]):
        step_layout["yaxis"]["range"] = [-1, 3]

    # ── Error constants ──
    N_type, (Kp_v, Kv_v, Ka_v) = error_constants(K * n, d)
    gain_info = f"K = {fmt(K)}  (Type {N_type}, Kp={fmt(Kp_v)}, Kv={fmt(Kv_v)}, Ka={fmt(Ka_v)})"

    comp_tf = tf_block("Gc(s)", nc, dc)

    before_data = _spec_dict(r0)
    before_data.update({"PM": _safe(PM0), "GM": _safe(GM0), "wgc": _safe(wg0)})
    after_data = _spec_dict(r1)
    after_data.update({"PM": _safe(PM1), "GM": _safe(GM1), "wgc": _safe(wg1)})

    note = None
    if PM1 is not None and PM1 < pm - 0.5:
        note = f"PM {PM1:.1f}° < target {pm:g}°; increase the safety margin."

    return {
        "plots": {
            "bode_mag": {"traces": bode_mag_traces, "layout": bode_mag_layout},
            "bode_phase": {"traces": bode_phase_traces, "layout": bode_phase_layout},
            "step_response": {"traces": step_traces, "layout": step_layout}
        },
        "results": {
            "comp_type": comp_type,
            "compensator_tf": comp_tf,
            "design_info": info,
            "gain_info": gain_info,
            "before": before_data,
            "after": after_data,
            "note": note
        }
    }
