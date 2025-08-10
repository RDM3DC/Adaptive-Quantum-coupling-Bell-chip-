#!/usr/bin/env python3
"""
Bell‑chip 4‑qubit chain (A‑B‑C‑D): triple‑sequence + phase sweep targeting high W‑state fidelity.

Design goals (mirrors your 2Q/3Q scripts):
- Segment-by-segment evolution with static H per segment (no ODE solver needed for the baseline run).
- Coupling model: XY (exchange) interaction on chosen edge: H_ij = J_eff * (X_i X_j + Y_i Y_j)/2
- Optional local Z‑phase offsets on targeted qubits to implement phase‑steer between segments.
- Hook for ARP scaling (multiplicative gain on J via lambda_c and a custom time‑gain function).
- Sweep utility over (J, lambda_c, phi, tau) with CSV output, plus single‑shot run for quick tests.

Targets from previous thread:
- 2Q: F_{Phi+} ~ 0.992 at t ~ 16.2 (synced in 2Q script)
- 3Q: F_GHZ > 0.85 near (J ~ 0.29, lambda_c ~ 0.88) with sequential AB/BC
- 4Q: Triple‑sequence (AB -> BC -> CD) + phase control; aim W4 fidelity > 0.8

Author: (generated)
"""
from __future__ import annotations
import numpy as np
from numpy import kron
from scipy.linalg import expm
import argparse
import itertools
import csv
from typing import Tuple, Dict, List

# ---------- Basic Pauli / tensor helpers ----------
I = np.eye(2, dtype=complex)
X = np.array([[0, 1],[1, 0]], dtype=complex)
Y = np.array([[0, -1j],[1j, 0]], dtype=complex)
Z = np.array([[1, 0],[0,-1]], dtype=complex)

PAULI = {'I': I, 'X': X, 'Y': Y, 'Z': Z}


def kronN(ops: List[np.ndarray]) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for op in ops:
        out = kron(out, op)
    return out


def op_on(n_qubits: int, which: int, op: np.ndarray) -> np.ndarray:
    """Return operator acting with `op` on qubit index `which` (0-based), identity elsewhere."""
    ops = [I]*n_qubits
    ops[which] = op
    return kronN(ops)


def two_body_xy(n_qubits: int, i: int, j: int) -> np.ndarray:
    """(X_i X_j + Y_i Y_j)/2 on an n-qubit register."""
    return 0.5*(op_on(n_qubits, i, X) @ op_on(n_qubits, j, X) +
                op_on(n_qubits, i, Y) @ op_on(n_qubits, j, Y))


# ---------- States & fidelity ----------

def basis_state(n_qubits: int, one_positions: List[int]) -> np.ndarray:
    """|b> with ones at the given positions (0-based from left=qubit0=A)."""
    N = 2**n_qubits
    idx = 0
    for q in range(n_qubits):
        bit = 1 if q in one_positions else 0
        idx = (idx << 1) | bit
    v = np.zeros((N,), dtype=complex)
    v[idx] = 1.0
    return v


def W4_state() -> np.ndarray:
    # |W4> = (|0001> + |0010> + |0100> + |1000>)/2
    n = 4
    terms = [basis_state(n, [3]),  # |...0001> (D)
             basis_state(n, [2]),  # C
             basis_state(n, [1]),  # B
             basis_state(n, [0])]  # A
    psi = sum(terms)
    psi = psi / np.linalg.norm(psi)
    return psi


def fidelity(psi: np.ndarray, phi: np.ndarray) -> float:
    return float(np.abs(np.vdot(phi, psi))**2)


def density_fidelity(rho: np.ndarray, phi: np.ndarray) -> float:
    return float(np.real(np.vdot(phi, rho @ phi)))


# ---------- ARP scaling hook ----------

def arp_gain(t: float, lambda_c: float, params: Dict) -> float:
    """Hook for ARP scaling. Default: smooth envelope in [0,1] multiplied by (1 + lambda_c).
    Replace with your ARP ODE gain if available.
    """
    # Smooth 0->1 envelope (raised-cosine) within each segment duration if provided
    T = params.get('segment_T', None)
    if T is None:
        env = 1.0
    else:
        x = np.clip(t / max(T, 1e-12), 0.0, 1.0)
        env = 0.5 - 0.5*np.cos(np.pi * x)  # cosine ramp
    return 1.0 + lambda_c * env


# ---------- Hamiltonians per segment ----------

def segment_hamiltonian_AB(J: float, lambda_c: float, phi: float, n=4, params: Dict|None=None) -> np.ndarray:
    # XY coupling on (A,B) = (0,1)
    Hij = two_body_xy(n, 0, 1)
    # local Z phase (steer): apply +phi/2 on A and -phi/2 on B as effective Z detunings
    Hloc = (phi/2.0)*op_on(n, 0, Z) + (-phi/2.0)*op_on(n, 1, Z)
    # Effective coupling
    J_eff = J * arp_gain(0.0, lambda_c, params or {})
    return J_eff*Hij + Hloc


def segment_hamiltonian_BC(J: float, lambda_c: float, phi: float, n=4, params: Dict|None=None) -> np.ndarray:
    Hij = two_body_xy(n, 1, 2)
    Hloc = (phi/2.0)*op_on(n, 1, Z) + (-phi/2.0)*op_on(n, 2, Z)
    J_eff = J * arp_gain(0.0, lambda_c, params or {})
    return J_eff*Hij + Hloc


def segment_hamiltonian_CD(J: float, lambda_c: float, phi: float, n=4, params: Dict|None=None) -> np.ndarray:
    Hij = two_body_xy(n, 2, 3)
    Hloc = (phi/2.0)*op_on(n, 2, Z) + (-phi/2.0)*op_on(n, 3, Z)
    J_eff = J * arp_gain(0.0, lambda_c, params or {})
    return J_eff*Hij + Hloc


# ---------- Evolution utilities ----------

def evolve_segment(psi: np.ndarray, H: np.ndarray, tau: float) -> np.ndarray:
    U = expm(-1j * H * tau)
    return U @ psi


def triple_sequence_run(J: float, lambda_c: float, phi: float, tau: float,
                        init: np.ndarray|None=None,
                        params: Dict|None=None) -> np.ndarray:
    """Run AB -> BC -> CD segments sequentially."""
    if init is None:
        init = basis_state(4, [])  # |0000>
    p = params or {}
    p = {**p, 'segment_T': tau}
    H1 = segment_hamiltonian_AB(J, lambda_c, phi, n=4, params=p)
    H2 = segment_hamiltonian_BC(J, lambda_c, phi, n=4, params=p)
    H3 = segment_hamiltonian_CD(J, lambda_c, phi, n=4, params=p)

    psi = evolve_segment(init, H1, tau)
    psi = evolve_segment(psi,  H2, tau)
    psi = evolve_segment(psi,  H3, tau)
    # Optional global phase cleanup (irrelevant to fidelity)
    return psi / np.linalg.norm(psi)


# ---------- W‑fidelity helper ----------

def compute_W4_fidelity(psi: np.ndarray) -> float:
    return fidelity(psi, W4_state())


# ---------- Sweep ----------

def sweep_grid(J_vals, lc_vals, phi_vals, tau_vals,
               csv_path: str = '4q_W_sweep.csv',
               topk: int = 20) -> List[Tuple[float,float,float,float,float]]:
    """Return top‑k rows and write full CSV: J,lambda_c,phi,tau,F_W."""
    header = ['J', 'lambda_c', 'phi', 'tau', 'F_W']
    rows = []
    for J, lc, phi, tau in itertools.product(J_vals, lc_vals, phi_vals, tau_vals):
        psi = triple_sequence_run(J, lc, phi, tau)
        FW = compute_W4_fidelity(psi)
        rows.append((J, lc, phi, tau, FW))
    # Write CSV
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    # Top‑k by fidelity
    rows.sort(key=lambda r: r[-1], reverse=True)
    best = rows[:topk]
    return best


# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(description='4Q bell‑chip triple‑sequence W‑fidelity sweep')
    sub = parser.add_subparsers(dest='cmd', required=False)

    p_run = sub.add_parser('run', help='single run')
    p_run.add_argument('--J', type=float, default=0.29)
    p_run.add_argument('--lambda_c', type=float, default=0.88)
    p_run.add_argument('--phi', type=float, default=0.0, help='phase steer [rad]')
    p_run.add_argument('--tau', type=float, default=5.4, help='segment duration')

    p_sweep = sub.add_parser('sweep', help='parameter sweep -> CSV')
    p_sweep.add_argument('--Jmin', type=float, default=0.20)
    p_sweep.add_argument('--Jmax', type=float, default=0.40)
    p_sweep.add_argument('--Jnum', type=int, default=9)
    p_sweep.add_argument('--LCmin', type=float, default=0.70)
    p_sweep.add_argument('--LCmax', type=float, default=1.00)
    p_sweep.add_argument('--LCnum', type=int, default=7)
    p_sweep.add_argument('--PHInum', type=int, default=13, help='phases in [0,2pi)')
    p_sweep.add_argument('--TAUmin', type=float, default=2.0)
    p_sweep.add_argument('--TAUmax', type=float, default=8.0)
    p_sweep.add_argument('--TAUnum', type=int, default=7)
    p_sweep.add_argument('--csv', type=str, default='4q_W_sweep.csv')
    p_sweep.add_argument('--topk', type=int, default=20)

    args = parser.parse_args()

    if args.cmd == 'run':
        psi = triple_sequence_run(args.J, args.lambda_c, args.phi, args.tau)
        FW = compute_W4_fidelity(psi)
        print(f"F_W4={FW:.6f}  (J={args.J}, lambda_c={args.lambda_c}, phi={args.phi}, tau={args.tau})")
        return

    # default: sweep
    J_vals  = np.linspace(getattr(args, 'Jmin', 0.20), getattr(args, 'Jmax', 0.40), getattr(args, 'Jnum', 9))
    LC_vals = np.linspace(getattr(args, 'LCmin', 0.70), getattr(args, 'LCmax', 1.00), getattr(args, 'LCnum', 7))
    PHI_vals = np.linspace(0.0, 2*np.pi, getattr(args, 'PHInum', 13), endpoint=False)
    TAU_vals = np.linspace(getattr(args, 'TAUmin', 2.0), getattr(args, 'TAUmax', 8.0), getattr(args, 'TAUnum', 7))

    best = sweep_grid(J_vals, LC_vals, PHI_vals, TAU_vals, csv_path=getattr(args, 'csv', '4q_W_sweep.csv'), topk=getattr(args, 'topk', 20))
    print("Top results (J, lambda_c, phi, tau, F_W):")
    for row in best:
        print(row)


if __name__ == '__main__':
    main()
