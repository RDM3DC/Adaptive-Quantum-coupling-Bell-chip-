#!/usr/bin/env python3
# bell_chip_three_qubit.py — 3-qubit Bell-chip (A–B–C chain) with shared memory. GHZ fidelity.
import numpy as np, matplotlib.pyplot as plt

# Pauli / helpers
I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
def kron3(a,b,c): return np.kron(np.kron(a,b),c)

# Single-qubit ops
sxA = kron3(sx, I2, I2); syA = kron3(sy, I2, I2); szA = kron3(sz, I2, I2)
sxB = kron3(I2, sx, I2); syB = kron3(I2, sy, I2); szB = kron3(I2, sz, I2)
sxC = kron3(I2, I2, sx); syC = kron3(I2, I2, sy); szC = kron3(I2, I2, sz)

# Couplings
sxsx_AB = kron3(sx, sx, I2)
sxsx_BC = kron3(I2, sx, sx)

# GHZ+
ghz = (1/np.sqrt(2))*np.array([1,0,0,0,0,0,0,1], dtype=complex)
def normalize(psi): return psi/np.linalg.norm(psi)
def F_ghz(psi): return float(np.abs(np.vdot(ghz, psi))**2)

def rhs(psi, g, p):
    H = 0.5*(p["DeltaA"]*szA + p["OmegaA"]*sxA) \
      + 0.5*(p["DeltaB"]*szB + p["OmegaB"]*sxB) \
      + 0.5*(p["DeltaC"]*szC + p["OmegaC"]*sxC) \
      + p["JAB"]*sxsx_AB + p["JBC"]*sxsx_BC \
      + p["lam"]*(g["gxA"]*sxA + g["gyA"]*syA + g["gxB"]*sxB + g["gyB"]*syB + g["gxC"]*sxC + g["gyC"]*syC) \
      + p["lamc"]*(g["gcAB"]*sxsx_AB + g["gcBC"]*sxsx_BC)
    dpsi = -1j*(H @ psi)
    # expectations
    xA = np.real(np.vdot(psi, sxA @ psi)); yA = np.real(np.vdot(psi, syA @ psi))
    xB = np.real(np.vdot(psi, sxB @ psi)); yB = np.real(np.vdot(psi, syB @ psi))
    xC = np.real(np.vdot(psi, sxC @ psi)); yC = np.real(np.vdot(psi, syC @ psi))
    xAB = np.real(np.vdot(psi, sxsx_AB @ psi))
    xBC = np.real(np.vdot(psi, sxsx_BC @ psi))
    dg = {}
    dg["gxA"] = 0.5*p["alpha"]*xA - p["mu"]*g["gxA"]; dg["gyA"] = 0.5*p["alpha"]*yA - p["mu"]*g["gyA"]
    dg["gxB"] = 0.5*p["alpha"]*xB - p["mu"]*g["gxB"]; dg["gyB"] = 0.5*p["alpha"]*yB - p["mu"]*g["gyB"]
    dg["gxC"] = 0.5*p["alpha"]*xC - p["mu"]*g["gxC"]; dg["gyC"] = 0.5*p["alpha"]*yC - p["mu"]*g["gyC"]
    dg["gcAB"] = p["alphac"]*xAB - p["muc"]*g["gcAB"]; dg["gcBC"] = p["alphac"]*xBC - p["muc"]*g["gcBC"]
    return dpsi, dg

def rk4_step(psi, g, dt, p):
    def addg(ga, gb, s): return {k:ga[k]+s*gb[k] for k in ga}
    k1psi, k1g = rhs(psi, g, p)
    k2psi, k2g = rhs(psi + 0.5*dt*k1psi, addg(g,k1g,0.5*dt), p)
    k3psi, k3g = rhs(psi + 0.5*dt*k2psi, addg(g,k2g,0.5*dt), p)
    k4psi, k4g = rhs(psi + dt*k3psi,     addg(g,k3g,dt),     p)
    psi_new = psi + (dt/6)*(k1psi + 2*k2psi + 2*k3psi + k4psi)
    g_new = {k: g[k] + (dt/6)*(k1g[k] + 2*k2g[k] + 2*k3g[k] + k4g[k]) for k in g}
    return normalize(psi_new), g_new

def simulate(T=20.0, dt=0.01, p=None, with_memory=True):
    if p is None:
        p = dict(DeltaA=0.0, DeltaB=0.0, DeltaC=0.0,
                 OmegaA=1.0, OmegaB=1.0, OmegaC=1.0,
                 JAB=0.18, JBC=0.18, lam=0.15, lamc=0.7,
                 alpha=0.6, mu=0.06, alphac=0.48, muc=0.06)
    psi = np.zeros(8, dtype=complex); psi[0]=1.0; psi = normalize(psi)
    g = {k:0.0 for k in ["gxA","gyA","gxB","gyB","gxC","gyC","gcAB","gcBC"]}
    if not with_memory:
        p = dict(**p); p["lam"]=0.0; p["lamc"]=0.0; p["alphac"]=0.0
    n=int(np.ceil(T/dt)); t=np.linspace(0,T,n+1); F=np.zeros(n+1)
    for i in range(n+1):
        F[i] = F_ghz(psi)
        if i<n: psi, g = rk4_step(psi, g, dt, p)
    return t, F

if __name__ == "__main__":
    t1, F0 = simulate(T=20.0, dt=0.01, with_memory=False)
    t2, Fm = simulate(T=20.0, dt=0.01, with_memory=True)
    plt.figure()
    plt.plot(t1, F0, label="No ARP memory"); plt.plot(t2, Fm, label="With ARP memory")
    plt.ylim(0,1.02); plt.xlabel("Time (Ω=1 units)"); plt.ylabel("F_GHZ(t)")
    plt.title("Three-qubit Bell-chip: GHZ fidelity"); plt.legend(); plt.tight_layout()
    plt.savefig("three_qubit_fidelity.png", dpi=160); plt.show()
