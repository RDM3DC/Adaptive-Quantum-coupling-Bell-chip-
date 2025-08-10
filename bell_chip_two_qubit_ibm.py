#!/usr/bin/env python3
"""
bell_chip_two_qubit_ibm.py
--------------------------
Two-qubit Bell-chip with shared memory — IBM transmon-style units.

Examples:
  # IBM-like: Δ=5 GHz, Ω=0.1 GHz, λ=0.5, μ=0.1α (rotating-frame; grid-search J, λc)
  python bell_chip_two_qubit_ibm.py --delta-ghz 5 --omega-ghz 0.1 --lam 0.5 --mu-ratio 0.1 --mode rotframe --grid

  # Keep lab detuning (Δ/Ω large) just to illustrate
  python bell_chip_two_qubit_ibm.py --delta-ghz 5 --omega-ghz 0.1 --lam 0.5 --mode labdelta
"""
import argparse, numpy as np, matplotlib.pyplot as plt

# Pauli ops
I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
def kron(a,b): return np.kron(a,b)
sxA = kron(sx, I2); syA = kron(sy, I2); szA = kron(sz, I2)
sxB = kron(I2, sx); syB = kron(I2, sy); szB = kron(I2, sz)
sxsx = kron(sx, sx)

phi_plus = (1/np.sqrt(2))*np.array([1,0,0,1], dtype=complex)
def fidelity_phi_plus(psi): return float(np.abs(np.vdot(phi_plus, psi))**2)
def normalize(psi): return psi/np.linalg.norm(psi)

def rhs(psi, gxA, gyA, gxB, gyB, gc, p):
    H = 0.5*(p["DeltaA"]*szA + p["OmegaA"]*sxA) + 0.5*(p["DeltaB"]*szB + p["OmegaB"]*sxB) \
        + p["J"]*sxsx + p["lam"]*(gxA*sxA + gyA*syA + gxB*sxB + gyB*syB) + p["lamc"]*gc*sxsx
    dpsi = -1j*(H @ psi)
    xA = np.real(np.vdot(psi, sxA @ psi)); yA = np.real(np.vdot(psi, syA @ psi))
    xB = np.real(np.vdot(psi, sxB @ psi)); yB = np.real(np.vdot(psi, syB @ psi))
    xAxB = np.real(np.vdot(psi, sxsx @ psi))
    dgxA = 0.5*p["alpha"]*xA - p["mu"]*gxA; dgyA = 0.5*p["alpha"]*yA - p["mu"]*gyA
    dgxB = 0.5*p["alpha"]*xB - p["mu"]*gxB; dgyB = 0.5*p["alpha"]*yB - p["mu"]*gyB
    dgc  = p["alphac"]*xAxB - p["muc"]*gc
    return dpsi, dgxA, dgyA, dgxB, dgyB, dgc

def rk4_step(state, dt, p):
    psi, gxA, gyA, gxB, gyB, gc = state
    k1 = rhs(psi, gxA, gyA, gxB, gyB, gc, p)
    k2 = rhs(psi + 0.5*dt*k1[0], gxA + 0.5*dt*k1[1], gyA + 0.5*dt*k1[2],
             gxB + 0.5*dt*k1[3], gyB + 0.5*dt*k1[4], gc + 0.5*dt*k1[5], p)
    k3 = rhs(psi + 0.5*dt*k2[0], gxA + 0.5*dt*k2[1], gyA + 0.5*dt*k2[2],
             gxB + 0.5*dt*k2[3], gyB + 0.5*dt*k2[4], gc + 0.5*dt*k2[5], p)
    k4 = rhs(psi + dt*k3[0], gxA + dt*k3[1], gyA + dt*k3[2],
             gxB + dt*k3[3], gyB + dt*k3[4], gc + dt*k3[5], p)
    psi_new = psi + (dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
    gxA_new = gxA + (dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    gyA_new = gyA + (dt/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    gxB_new = gxB + (dt/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
    gyB_new = gyB + (dt/6)*(k1[4]+2*k2[4]+2*k3[4]+k4[4])
    gc_new  = gc  + (dt/6)*(k1[5]+2*k2[5]+2*k3[5]+k4[5])
    return (normalize(psi_new), gxA_new, gyA_new, gxB_new, gyB_new, gc_new)

def simulate(T, dt, p):
    state = (np.array([1,0,0,0], dtype=complex), 0.0,0.0,0.0,0.0,0.0)
    state = (state[0]/np.linalg.norm(state[0]), *state[1:])
    n = int(np.ceil(T/dt)); times = np.linspace(0, T, n+1)
    F = np.zeros(n+1)
    for i in range(n+1):
        F[i] = float(np.abs(np.vdot(phi_plus, state[0]))**2)
        if i < n: state = rk4_step(state, dt, p)
    return times, F

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta-ghz", type=float, default=5.0)
    ap.add_argument("--omega-ghz", type=float, default=0.1)
    ap.add_argument("--lam", type=float, default=0.5, help="local λ")
    ap.add_argument("--alpha", type=float, default=0.6)
    ap.add_argument("--mu-ratio", type=float, default=0.1, help="μ = mu_ratio * α")
    ap.add_argument("--lamc", type=float, default=0.6, help="shared λ_c (used if --grid not set)")
    ap.add_argument("--J", type=float, default=0.15, help="XX coupling (normalized by Ω) if --grid not set")
    ap.add_argument("--mode", choices=["rotframe","labdelta"], default="rotframe")
    ap.add_argument("--grid", action="store_true", help="grid search J, λc to maximize F_max")
    args = ap.parse_args()

    Delta_model = args.delta_ghz / args.omega_ghz
    alpha = args.alpha; mu = args.alpha * args.mu_ratio; alphac = alpha*0.8; muc = mu
    if args.mode == "rotframe":
        p = dict(DeltaA=0.0, DeltaB=0.0, OmegaA=1.0, OmegaB=1.0, J=args.J, lam=args.lam,
                 lamc=args.lamc, alpha=alpha, mu=mu, alphac=alphac, muc=muc)
    else:
        p = dict(DeltaA=Delta_model, DeltaB=Delta_model, OmegaA=1.0, OmegaB=1.0, J=args.J, lam=args.lam,
                 lamc=args.lamc, alpha=alpha, mu=mu, alphac=alphac, muc=muc)

    if args.grid and args.mode=="rotframe":
        J_grid = np.linspace(0.05, 0.35, 7)
        lamc_grid = np.linspace(0.2, 0.9, 8)
        best = (0.0, None, None, None, None)
        for J in J_grid:
            for lamc in lamc_grid:
                p_try = dict(**p, J=float(J), lamc=float(lamc))
                t, F = simulate(40.0, 0.004, p_try)
                Fm = float(F.max())
                if Fm > best[0]:
                    best = (Fm, J, lamc, t, F)
        Fmax, Jopt, Lcopt, t, F = best
        print(f"[rotframe grid] best F_max={Fmax:.3f} at J={Jopt:.2f}, λc={Lcopt:.2f}")
        # First-peak time (Ω=1 units and ns for Ω_GHz)
        idx = np.where(F > 0.9)[0]
        if len(idx)>0:
            t_peak = t[idx[0]]
            print(f"First peak at t≈{t_peak:.2f} (Ω=1 units) ≈ {t_peak*1e1:.1f} ns for Ω={args.omega_ghz} GHz")
        plt.figure(); plt.plot(t, F); plt.ylim(0,1.02)
        plt.xlabel("Time (Ω=1 units)"); plt.ylabel("F_Φ+(t)")
        plt.title(f"Δ_eff=0, best: F_max={Fmax:.3f} [J={Jopt:.2f}, λc={Lcopt:.2f}]")
        plt.tight_layout(); plt.savefig("two_qubit_rotframe_fidelity.png", dpi=150); plt.show()
    else:
        t, F = simulate(40.0, 0.004, p)
        Fmax = float(F.max())
        label = "Δ_eff=0" if args.mode=="rotframe" else f"lab Δ/Ω={Delta_model:.1f}"
        print(f"[{label}] F_max={Fmax:.3f}")
        plt.figure(); plt.plot(t, F); plt.ylim(0,1.02)
        plt.xlabel("Time (Ω=1 units)"); plt.ylabel("F_Φ+(t)")
        plt.title(f"{label}, F_max={Fmax:.3f} (J={p['J']:.2f}, λc={p['lamc']:.2f})")
        plt.tight_layout(); plt.savefig("two_qubit_fidelity.png", dpi=150); plt.show()

if __name__ == "__main__":
    main()
