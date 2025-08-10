#!/usr/bin/env python3
# bell_chip_two_qubit.py — Simulate two-qubit Bell-chip with shared memory and plot concurrence.
import numpy as np, matplotlib.pyplot as plt

I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
def kron(a,b): return np.kron(a,b)
sxA = kron(sx, I2); syA = kron(sy, I2); szA = kron(sz, I2)
sxB = kron(I2, sx); syB = kron(I2, sy); szB = kron(I2, sz)
sxsx = kron(sx, sx)

def concurrence(psi):
    a,b,c,d = psi
    return float(abs(2*(a*d - b*c)))

def normalize(psi): return psi/np.linalg.norm(psi)

def rhs(psi, gxA, gyA, gxB, gyB, gc, p, t):
    H = 0.5*(p["DeltaA"]*szA + p["OmegaA"]*sxA) \
        + 0.5*(p["DeltaB"]*szB + p["OmegaB"]*sxB) \
        + p["J"]*sxsx \
        + p["lam"]*(gxA*sxA + gyA*syA + gxB*sxB + gyB*syB) \
        + p["lamc"]*gc*sxsx
    dpsi = -1j*(H @ psi)
    xA = np.real(np.vdot(psi, sxA @ psi)); yA = np.real(np.vdot(psi, syA @ psi))
    xB = np.real(np.vdot(psi, sxB @ psi)); yB = np.real(np.vdot(psi, syB @ psi))
    xAxB = np.real(np.vdot(psi, sxsx @ psi))
    dgxA = 0.5*p["alpha"]*xA - p["mu"]*gxA
    dgyA = 0.5*p["alpha"]*yA - p["mu"]*gyA
    dgxB = 0.5*p["alpha"]*xB - p["mu"]*gxB
    dgyB = 0.5*p["alpha"]*yB - p["mu"]*gyB
    dgc  = p["alphac"]*xAxB - p["muc"]*gc
    return dpsi, dgxA, dgyA, dgxB, dgyB, dgc

def rk4_step(state, t, dt, p):
    psi, gxA, gyA, gxB, gyB, gc = state
    def f(s, tt): return rhs(s[0], s[1], s[2], s[3], s[4], s[5], p, tt)
    k1 = f((psi, gxA, gyA, gxB, gyB, gc), t)
    k2 = f((psi + 0.5*dt*k1[0], gxA + 0.5*dt*k1[1], gyA + 0.5*dt*k1[2],
            gxB + 0.5*dt*k1[3], gyB + 0.5*dt*k1[4], gc + 0.5*dt*k1[5]), t+0.5*dt)
    k3 = f((psi + 0.5*dt*k2[0], gxA + 0.5*dt*k2[1], gyA + 0.5*dt*k2[2],
            gxB + 0.5*dt*k2[3], gyB + 0.5*dt*k2[4], gc + 0.5*dt*k2[5]), t+0.5*dt)
    k4 = f((psi + dt*k3[0],    gxA + dt*k3[1],    gyA + dt*k3[2],
            gxB + dt*k3[3],    gyB + dt*k3[4],    gc + dt*k3[5]), t+dt)
    psi_new = psi + (dt/6)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
    gxA_new = gxA + (dt/6)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    gyA_new = gyA + (dt/6)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
    gxB_new = gxB + (dt/6)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])
    gyB_new = gyB + (dt/6)*(k1[4]+2*k2[4]+2*k3[4]+k4[4])
    gc_new  = gc  + (dt/6)*(k1[5]+2*k2[5]+2*k3[5]+k4[5])
    return (normalize(psi_new), gxA_new, gyA_new, gxB_new, gyB_new, gc_new)

def simulate(T=30.0, dt=0.002, p=None, psi0=None):
    if p is None:
        p = dict(DeltaA=0.0, DeltaB=0.0, OmegaA=1.0, OmegaB=1.0, J=0.15,
                 lam=0.8, lamc=0.6, alpha=0.6, mu=0.2, alphac=0.4, muc=0.15)
    if psi0 is None:
        psi0 = np.array([1,0,0,0], dtype=complex)
    state = (psi0/np.linalg.norm(psi0), 0.0,0.0,0.0,0.0,0.0)
    n=int(np.ceil(T/dt)); t=np.linspace(0,T,n+1); C=np.zeros(n+1)
    for i,ti in enumerate(t):
        C[i] = float(abs(2*(state[0][0]*state[0][3]-state[0][1]*state[0][2])))
        if i<n: state = rk4_step(state, ti, dt, p)
    return t, C

if __name__ == "__main__":
    T, dt = 30.0, 0.002
    p_with = dict(DeltaA=0.0, DeltaB=0.0, OmegaA=1.0, OmegaB=1.0, J=0.15,
                  lam=0.8, lamc=0.6, alpha=0.6, mu=0.2, alphac=0.4, muc=0.15)
    p_no   = dict(DeltaA=0.0, DeltaB=0.0, OmegaA=1.0, OmegaB=1.0, J=0.15,
                  lam=0.0, lamc=0.0, alpha=0.6, mu=0.2, alphac=0.0, muc=0.15)
    t, C_with = simulate(T, dt, p_with)
    _, C_no   = simulate(T, dt, p_no)
    plt.figure()
    plt.plot(t, C_no, label="Concurrence (λ=0, no shared)")
    plt.plot(t, C_with, label="Concurrence (λ>0, shared)")
    plt.xlabel("Time"); plt.ylabel("Concurrence")
    plt.title("Two-qubit Bell-chip: Entanglement Gain with Shared Memory")
    plt.legend(); plt.tight_layout()
    plt.savefig("two_qubit_concurrence.png", dpi=150)
    plt.show()
