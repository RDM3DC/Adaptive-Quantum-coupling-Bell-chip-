#!/usr/bin/env python3
# Parameter map for Bell-chip ARP hysteresis vs (Lambda, Delta)
# Saves heatmap to ./outputs/parameter_map_hysteresis.png

import os, numpy as np, matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class Params:
    Delta: float
    Omega: float
    alpha: float
    mu: float
    lam: float

def deriv(state, t, p):
    x, y, z, gx, gy = state
    dx = -p.Delta*y - 2*p.lam*gy*z
    dy =  p.Delta*x - p.Omega*z + 2*p.lam*gx*z
    dz =  p.Omega*y + 2*p.lam*(gy*x - gx*y)
    dgx = 0.5*p.alpha*x - p.mu*gx
    dgy = 0.5*p.alpha*y - p.mu*gy
    return np.array([dx, dy, dz, dgx, dgy])

def rk4_step(f, y, t, dt, p):
    k1 = f(y, t, p)
    k2 = f(y + 0.5*dt*k1, t + 0.5*dt, p)
    k3 = f(y + 0.5*dt*k2, t + 0.5*dt, p)
    k4 = f(y + dt*k3, t + dt, p)
    return y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def settle_ss(p, state_init, T=8.0, dt=0.003):
    n = int(T/dt); y = state_init.copy(); t=0.0
    for _ in range(n):
        y = rk4_step(deriv, y, t, dt, p); t += dt
    return y

def hysteresis_metric(alpha=0.6, mu=0.2, Omega=1.0, lam=0.8, deltas=None):
    if deltas is None:
        deltas = np.linspace(-2,2,31)
    z_up=[]; s = np.array([0,0,1,0,0], float)
    for D in deltas:
        p = Params(Delta=float(D), Omega=Omega, alpha=alpha, mu=mu, lam=lam)
        s = settle_ss(p, s, T=8.0, dt=0.003); z_up.append(s[2])
    z_dn=[] 
    for D in deltas[::-1]:
        p = Params(Delta=float(D), Omega=Omega, alpha=alpha, mu=mu, lam=lam)
        s = settle_ss(p, s, T=8.0, dt=0.003); z_dn.append(s[2])
    z_dn_rev = z_dn[::-1]
    diff = np.abs(np.array(z_up) - np.array(z_dn_rev))
    area = np.trapz(diff, deltas) / (deltas[-1]-deltas[0])
    return area

def main():
    outdir = "./outputs"; os.makedirs(outdir, exist_ok=True)
    alpha, mu, Omega = 0.6, 0.2, 1.0
    lam_vals = np.linspace(0.0, 1.2, 10)
    Delta_offsets = np.linspace(-1.0, 1.0, 21)
    base_deltas = np.linspace(-2,2,31)
    heat = np.zeros((len(lam_vals), len(Delta_offsets)))
    for i, lam in enumerate(lam_vals):
        for j, d0 in enumerate(Delta_offsets):
            deltas = base_deltas + d0
            heat[i,j] = hysteresis_metric(alpha=alpha, mu=mu, Omega=Omega, lam=lam, deltas=deltas)
    Lambda_eff = lam_vals * alpha / mu
    plt.figure()
    extent = [Delta_offsets[0], Delta_offsets[-1], Lambda_eff[0], Lambda_eff[-1]]
    plt.imshow(heat, origin='lower', aspect='auto', extent=extent)
    plt.colorbar(label="Hysteresis metric (area)")
    plt.xlabel("Δ sweep center offset")
    plt.ylabel("Λ = λ α / μ")
    plt.title("Parameter map: Hysteresis strength vs (Λ, Δ offset)")
    plt.tight_layout()
    path = os.path.join(outdir, "parameter_map_hysteresis.png")
    plt.savefig(path, dpi=150); print("Saved:", path)

if __name__ == "__main__":
    main()
