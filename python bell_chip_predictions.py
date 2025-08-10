#!/usr/bin/env python3
"""
Bell-chip ARP Predictions
=========================
Run falsifiable predictions (P1–P3) for the Bloch–ARP model from the README.

Usage examples:
  python bell_chip_predictions.py --which P1
  python bell_chip_predictions.py --which P2 --Gamma2 0.03 --eta 0.8
  python bell_chip_predictions.py --which P3 --light
  python bell_chip_predictions.py --which ALL --out ./outputs

Notes:
- Uses RK4 integration.
- Each prediction renders exactly one figure.
- Requires: numpy, matplotlib, dataclasses (py3.7+).
"""

import argparse
import os
from dataclasses import dataclass
from typing import Callable, Tuple

import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# Core Bloch–ARP dynamics
# ------------------------

@dataclass
class Params:
    Delta: float
    Omega: float
    alpha: float
    mu: float
    lam: float
    Gamma2: float = 0.0
    eta: float = 0.0

def deriv_bloch_arp(state: np.ndarray, t: float, p: Params,
                    Omega_of_t: Callable[[float], float] = None) -> np.ndarray:
    """
    state = [x, y, z, g_x, g_y]
    dx/dt = -Δ y - 2 λ g_y z - Gamma2_eff * x
    dy/dt =  Δ x - Ω z + 2 λ g_x z - Gamma2_eff * y
    dz/dt =  Ω y + 2 λ (g_y x - g_x y)
    dgx/dt = 0.5 α x - μ g_x
    dgy/dt = 0.5 α y - μ g_y
    where Gamma2_eff = p.Gamma2 + p.eta * g_y^2
    """
    x, y, z, g_x, g_y = state
    Delta = p.Delta
    Omega = p.Omega if Omega_of_t is None else Omega_of_t(t)
    Gamma2_eff = p.Gamma2 + p.eta * (g_y ** 2)

    dx = -Delta * y - 2 * p.lam * g_y * z - Gamma2_eff * x
    dy = Delta * x - Omega * z + 2 * p.lam * g_x * z - Gamma2_eff * y
    dz = Omega * y + 2 * p.lam * (g_y * x - g_x * y)
    dgx = 0.5 * p.alpha * x - p.mu * g_x
    dgy = 0.5 * p.alpha * y - p.mu * g_y
    return np.array([dx, dy, dz, dgx, dgy], dtype=float)

def rk4_step(f, y, t, dt, *args, **kwargs):
    k1 = f(y, t, *args, **kwargs)
    k2 = f(y + 0.5 * dt * k1, t + 0.5 * dt, *args, **kwargs)
    k3 = f(y + 0.5 * dt * k2, t + 0.5 * dt, *args, **kwargs)
    k4 = f(y + dt * k3, t + dt, *args, **kwargs)
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

def integrate(f, y0, t0, t1, dt, *args, **kwargs):
    n = int(np.ceil((t1 - t0) / dt))
    y = np.array(y0, dtype=float)
    t = t0
    traj = np.zeros((n + 1, len(y0)))
    ts = np.zeros(n + 1)
    traj[0] = y
    ts[0] = t
    for i in range(n):
        y = rk4_step(f, y, t, dt, *args, **kwargs)
        t = t0 + (i + 1) * dt
        traj[i + 1] = y
        ts[i + 1] = t
    return ts, traj

# ------------------------
# P1: Hysteresis in z*(Δ)
# ------------------------

def run_P1_hysteresis(outdir: str, base: Params, dt: float = 0.002, T_settle: float = 15.0):
    Deltas_up = np.linspace(-2.0, 2.0, 41)
    Deltas_dn = Deltas_up[::-1]

    def settle_to_ss(Delta, state_init):
        p = Params(**{**base.__dict__, "Delta": float(Delta)})
        _, tr = integrate(deriv_bloch_arp, state_init, 0.0, T_settle, dt, p)
        return tr[-1]

    state0 = np.array([0.0, 0.0, 1.0, 0.0, 0.0], dtype=float)
    z_up = []
    s = state0.copy()
    for D in Deltas_up:
        s = settle_to_ss(D, s)
        z_up.append(s[2])

    z_dn = []
    for D in Deltas_dn:
        s = settle_to_ss(D, s)
        z_dn.append(s[2])

    plt.figure()
    plt.plot(Deltas_up, z_up, label="Up sweep")
    plt.plot(Deltas_dn, z_dn, label="Down sweep")
    plt.axhline(0, linestyle="--", linewidth=0.8)
    plt.xlabel("Detuning Δ")
    plt.ylabel("Steady-state z*")
    plt.title("P1: Hysteresis in z*(Δ) for Λ=λα/μ>0")
    plt.legend()
    plt.tight_layout()
    path = os.path.join(outdir, "P1_hysteresis.png")
    plt.savefig(path, dpi=150)
    plt.close()
    return path

# ------------------------------------------------------------
# P2: Controllable decoherence floor — T2 vs drive amplitude
# ------------------------------------------------------------

def estimate_T2(times, y_signal):
    eps = 1e-12
    y_abs = np.abs(y_signal) + eps
    logy = np.log(y_abs)
    A = np.vstack([np.ones_like(times), times]).T
    coeff, *_ = np.linalg.lstsq(A, logy, rcond=None)
    b = coeff[1]
    if b >= 0:
        return np.inf
    return -1.0 / b

def run_P2_decoherence(outdir: str, base: Params, dt: float = 0.001,
                       T_build: float = 20.0, T_decay: float = 8.0):
    amps = np.linspace(0.2, 2.0, 10)
    T2_est = []

    for A in amps:
        p_build = Params(**{**base.__dict__, "Omega": float(A)})
        state0 = np.array([0.0, 0.0, 1.0, 0.0, 0.0], dtype=float)
        _, tr1 = integrate(deriv_bloch_arp, state0, 0.0, T_build, dt, p_build)
        s_switch = tr1[-1]
        p_decay = Params(**{**base.__dict__, "Omega": 0.0})
        t2, tr2 = integrate(deriv_bloch_arp, s_switch, 0.0, T_decay, dt, p_decay)
        k0 = int(0.05 * len(t2))
        T2_val = estimate_T2(t2[k0:], tr2[k0:, 1])  # use y component
        T2_est.append(T2_val)

    plt.figure()
    plt.plot(amps, T2_est, marker="o")
    plt.xlabel("Drive amplitude A")
    plt.ylabel("Estimated T2 (from FID of y)")
    plt.title("P2: Controllable Decoherence Floor via η g_y^2")
    plt.tight_layout()
    path = os.path.join(outdir, "P2_T2_vs_amplitude.png")
    plt.savefig(path, dpi=150)
    plt.close()
    return path

# ------------------------------------------------------------
# P3: Line-shape skew — lock-in response vs detuning
# ------------------------------------------------------------

def lockin_amplitude(times, signal, omega):
    cos_ref = np.cos(omega * times)
    sin_ref = np.sin(omega * times)
    sig = signal - np.mean(signal)
    X = 2.0 / len(times) * np.sum(sig * cos_ref)
    Y = 2.0 / len(times) * np.sum(sig * sin_ref)
    R = np.sqrt(X ** 2 + Y ** 2)
    return R, X, Y

def run_P3_lineshape(outdir: str, base: Params, A: float = 0.6, omega: float = 1.0,
                     dt: float = 0.001, T_total: float = 80.0, T_skip: float = 30.0,
                     light: bool = False):
    if light:
        dt = 0.002
        T_total = 50.0
        T_skip = 20.0
        Deltas = np.linspace(-2.0, 2.0, 21)
    else:
        Deltas = np.linspace(-2.0, 2.0, 41)

    amps_noARP, amps_ARP = [], []
    base_noARP = Params(**{**base.__dict__, "lam": 0.0})
    base_ARP = base

    def Omega_of_t(t):
        return A * np.cos(omega * t)

    for D in Deltas:
        # No ARP (baseline symmetric line)
        p = Params(**{**base_noARP.__dict__, "Delta": float(D)})
        s0 = np.array([0.0, 0.0, 1.0, 0.0, 0.0], dtype=float)
        t, tr = integrate(deriv_bloch_arp, s0, 0.0, T_total, dt, p, Omega_of_t)
        mask = t >= T_skip
        R, _, _ = lockin_amplitude(t[mask], tr[mask, 2], omega)
        amps_noARP.append(R)

        # With ARP
        p = Params(**{**base_ARP.__dict__, "Delta": float(D)})
        s0 = np.array([0.0, 0.0, 1.0, 0.0, 0.0], dtype=float)
        t, tr = integrate(deriv_bloch_arp, s0, 0.0, T_total, dt, p, Omega_of_t)
        mask = t >= T_skip
        R, _, _ = lockin_amplitude(t[mask], tr[mask, 2], omega)
        amps_ARP.append(R)

    plt.figure()
    plt.plot(Deltas, amps_noARP, label="No ARP (λ=0)")
    plt.plot(Deltas, amps_ARP, label="With ARP (λ>0)")
    plt.xlabel("Detuning Δ")
    plt.ylabel("|Z(ω)| (lock-in amplitude)")
    title = "P3: Line-shape — ARP-induced dispersive skew"
    if light:
        title += " (light)"
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    path = os.path.join(outdir, "P3_lineshape_skew.png")
    plt.savefig(path, dpi=150)
    plt.close()
    return path

# ------------------
# CLI
# ------------------

def main():
    ap = argparse.ArgumentParser(description="Bell-chip ARP Predictions (P1–P3)")
    ap.add_argument("--which", choices=["P1", "P2", "P3", "ALL"], default="ALL",
                    help="Which prediction(s) to run.")
    ap.add_argument("--out", default="./outputs", help="Output directory for figures.")
    # Common params
    ap.add_argument("--Delta", type=float, default=0.0, help="Base detuning (some runs sweep Δ).")
    ap.add_argument("--Omega", type=float, default=1.0, help="Base Rabi amplitude.")
    ap.add_argument("--alpha", type=float, default=0.6, help="Reinforcement rate.")
    ap.add_argument("--mu", type=float, default=0.2, help="Memory decay rate.")
    ap.add_argument("--lam", type=float, default=0.8, help="Coupling λ.")
    ap.add_argument("--Gamma2", type=float, default=0.02, help="Base dephasing 1/T2 (for P2/P3).")
    ap.add_argument("--eta", type=float, default=0.6, help="Dephasing coefficient for g_y^2 (P2).")
    # P3 tweaks
    ap.add_argument("--light", action="store_true", help="Faster P3 run with coarser grid.")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Base parameters for runs; individual routines override specifics as needed.
    base = Params(
        Delta=args.Delta, Omega=args.Omega, alpha=args.alpha, mu=args.mu,
        lam=args.lam, Gamma2=args.Gamma2, eta=args.eta
    )

    saved = []

    if args.which in ("P1", "ALL"):
        # For P1 we prefer zero dephasing to isolate hysteresis mechanism
        p1_base = Params(**{**base.__dict__, "Gamma2": 0.0, "eta": 0.0, "Omega": 1.0})
        saved.append(run_P1_hysteresis(args.out, p1_base))

    if args.which in ("P2", "ALL"):
        # P2 needs dephasing and eta>0 for amplitude-dependent floor
        p2_base = Params(**{**base.__dict__, "Delta": 0.0})
        saved.append(run_P2_decoherence(args.out, p2_base))

    if args.which in ("P3", "ALL"):
        # P3 compares λ=0 vs λ>0; base Gamma2 used; light mode optional
        p3_base = Params(**{**base.__dict__, "mu": max(0.15, base.mu), "lam": max(0.0, base.lam)})
        saved.append(run_P3_lineshape(args.out, p3_base, light=args.light))

    print("Saved:")
    for p in saved:
        print(" -", p)

if __name__ == "__main__":
    main()
