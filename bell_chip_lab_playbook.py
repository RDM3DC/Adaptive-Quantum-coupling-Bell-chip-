#!/usr/bin/env python3
# Bell-chip Lab Playbook Helper
# Choose λ from target Λ, run light checks for P1–P3, export CSV + PNG.

import os, csv, numpy as np, matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Callable

@dataclass
class Params:
    Delta: float
    Omega: float
    alpha: float
    mu: float
    lam: float
    Gamma2: float=0.02
    eta: float=0.6

def deriv(state, t, p, Omega_of_t=None):
    x, y, z, gx, gy = state
    Omega = p.Omega if Omega_of_t is None else Omega_of_t(t)
    Gamma2_eff = p.Gamma2 + p.eta*(gy**2)
    dx = -p.Delta*y - 2*p.lam*gy*z - Gamma2_eff*x
    dy =  p.Delta*x - Omega*z + 2*p.lam*gx*z - Gamma2_eff*y
    dz =  Omega*y + 2*p.lam*(gy*x - gx*y)
    dgx = 0.5*p.alpha*x - p.mu*gx
    dgy = 0.5*p.alpha*y - p.mu*gy
    return np.array([dx, dy, dz, dgx, dgy])

def rk4_step(f, y, t, dt, *args, **kwargs):
    k1 = f(y, t, *args, **kwargs)
    k2 = f(y + 0.5*dt*k1, t + 0.5*dt, *args, **kwargs)
    k3 = f(y + 0.5*dt*k2, t + 0.5*dt, *args, **kwargs)
    k4 = f(y + dt*k3, t + dt, *args, **kwargs)
    return y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def integrate(f, y0, t0, t1, dt, *args, **kwargs):
    n = int((t1 - t0) / dt)
    y = np.array(y0, float); t=t0
    traj = np.zeros((n+1, len(y0))); ts = np.zeros(n+1)
    traj[0]=y; ts[0]=t
    for i in range(n):
        y = rk4_step(f, y, t, dt, *args, **kwargs)
        t += dt; traj[i+1]=y; ts[i+1]=t
    return ts, traj

def lockin_amplitude(times, signal, omega):
    sig = signal - np.mean(signal)
    X = 2/len(times)*np.sum(sig*np.cos(omega*times))
    Y = 2/len(times)*np.sum(sig*np.sin(omega*times))
    R = np.sqrt(X**2+Y**2)
    return R, X, Y

def lambda_from_Lambda(Lambda, alpha, mu):
    return Lambda * mu / alpha

def quick_verify(outdir, target_Lambda=2.4, alpha=0.6, mu=0.2, 
                 Delta_center=0.0, Omega=1.0, Gamma2=0.02, eta=0.6):
    os.makedirs(outdir, exist_ok=True)
    lam = lambda_from_Lambda(target_Lambda, alpha, mu)
    # Save chosen parameters
    with open(os.path.join(outdir, "chosen_params.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["param","value"])
        for k,v in dict(target_Lambda=target_Lambda, alpha=alpha, mu=mu, lam=lam, 
                        Delta_center=Delta_center, Omega=Omega, Gamma2=Gamma2, eta=eta).items():
            w.writerow([k,v])

    # P1 light hysteresis
    p = Params(Delta=0.0, Omega=Omega, alpha=alpha, mu=mu, lam=lam, Gamma2=0.0, eta=0.0)
    deltas = np.linspace(-1.6, 1.6, 33) + Delta_center
    dt=0.003; T=8.0
    def settle(Delta, s0):
        pp = Params(**{**p.__dict__, "Delta": float(Delta)})
        _, tr = integrate(deriv, s0, 0.0, T, dt, pp)
        return tr[-1]
    s = np.array([0,0,1,0,0], float)
    z_up=[]; 
    for D in deltas: s = settle(D, s); z_up.append(s[2])
    z_dn=[]; 
    for D in deltas[::-1]: s = settle(D, s); z_dn.append(s[2])
    z_dn = z_dn[::-1]
    area = np.trapz(np.abs(np.array(z_up)-np.array(z_dn)), deltas)/(deltas[-1]-deltas[0])
    import matplotlib.pyplot as plt
    plt.figure(); plt.plot(deltas, z_up, label="up"); plt.plot(deltas, z_dn, label="down")
    plt.xlabel("Δ"); plt.ylabel("z*"); plt.title(f"P1 light | hysteresis area≈{area:.3f}")
    plt.legend(); plt.tight_layout()
    p1_path = os.path.join(outdir,"P1_light.png"); plt.savefig(p1_path, dpi=140); plt.close()

    # P2 T2 vs amplitude (3 points)
    amps = [0.4, 1.0, 1.6]; T_build=12.0; T_decay=6.0; dt=0.002
    rows=[["A","T2_est"]]
    def estimate_T2(times, y):
        eps=1e-12; logy=np.log(np.abs(y)+eps)
        A=np.vstack([np.ones_like(times), times]).T
        b=np.linalg.lstsq(A, logy, rcond=None)[0][1]
        return np.inf if b>=0 else -1.0/b
    for Aamp in amps:
        pb = Params(Delta=0.0, Omega=Aamp, alpha=alpha, mu=mu, lam=lam, Gamma2=Gamma2, eta=eta)
        s0=np.array([0,0,1,0,0], float)
        _, tr1 = integrate(deriv, s0, 0.0, T_build, dt, pb)
        s1=tr1[-1]
        pd = Params(Delta=0.0, Omega=0.0, alpha=alpha, mu=mu, lam=lam, Gamma2=Gamma2, eta=eta)
        t2, tr2 = integrate(deriv, s1, 0.0, T_decay, dt, pd)
        k0=int(0.1*len(t2)); T2=estimate_T2(t2[k0:], tr2[k0:,1]); rows.append([Aamp, T2])
    with open(os.path.join(outdir,"P2_T2_summary.csv"),"w",newline="") as f: csv.writer(f).writerows(rows)
    plt.figure(); plt.plot([r[0] for r in rows[1:]], [r[1] for r in rows[1:]], marker="o")
    plt.xlabel("Drive amplitude A"); plt.ylabel("Estimated T2"); plt.title("P2 light")
    plt.tight_layout(); p2_path=os.path.join(outdir,"P2_light.png"); plt.savefig(p2_path, dpi=140); plt.close()

    # P3 lineshape (very light)
    base = Params(Delta=0.0, Omega=0.0, alpha=alpha, mu=mu, lam=lam, Gamma2=0.015, eta=0.0)
    Aamp=0.6; omega=1.0; dt=0.002; Ttot=40.0; Tskip=16.0
    D = np.linspace(-1.6, 1.6, 19)
    def Omega_of_t(t): return Aamp*np.cos(omega*t)
    amps_ARP=[]; amps_no=[]

    # no-ARP baseline
    base_no = Params(**{**base.__dict__, "lam": 0.0})
    for d in D:
        s0=np.array([0,0,1,0,0], float)
        t,tr = integrate(deriv, s0, 0.0, Ttot, dt, Params(**{**base_no.__dict__, "Delta":float(d)}), Omega_of_t)
        m=t>=Tskip; R,_X,_Y = lockin_amplitude(t[m], tr[m,2], omega); amps_no.append(R)
    # with ARP
    for d in D:
        s0=np.array([0,0,1,0,0], float)
        t,tr = integrate(deriv, s0, 0.0, Ttot, dt, Params(**{**base.__dict__, "Delta":float(d)}), Omega_of_t)
        m=t>=Tskip; R,_X,_Y = lockin_amplitude(t[m], tr[m,2], omega); amps_ARP.append(R)
    plt.figure(); plt.plot(D, amps_no, label="λ=0"); plt.plot(D, amps_ARP, label="λ>0")
    plt.xlabel("Δ"); plt.ylabel("|Z(ω)|"); plt.title("P3 light"); plt.legend(); plt.tight_layout()
    p3_path=os.path.join(outdir,"P3_light.png"); plt.savefig(p3_path, dpi=140); plt.close()

    return dict(params_csv=os.path.join(outdir,"chosen_params.csv"),
                p1=p1_path, p2=p2_path, p3=p3_path,
                p2_csv=os.path.join(outdir,"P2_T2_summary.csv"),
                lambda_value=lam, hysteresis_area=area)

if __name__ == "__main__":
    out = quick_verify("./bell_chip_lab_playbook")
    print("Saved files:")
    for k,v in out.items():
        if isinstance(v,str): print(" -", k, ":", v)
    print("λ chosen from Λ target:", out["lambda_value"])
    print("Hysteresis metric (area):", out["hysteresis_area"])
