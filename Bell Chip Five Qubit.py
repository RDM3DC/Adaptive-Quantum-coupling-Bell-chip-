#!/usr/bin/env python3
"""
Bell‑chip 5‑qubit chain (A‑B‑C‑D‑E): sequential + phase cascade targeting high GHZ5 fidelity.

Design mirrors 2Q/3Q/4Q scaffolds:
- Segment evolution with static H per segment, XY exchange on each active edge.
- Phase steer via local Z terms; ARP scaling hook for time‑envelope gain on J.
- Sweep over (J, lambda_c, phi, tau) with CSV output and top‑k report.

Targets based on thread:
- 2Q: F_{Phi+} ~ 0.992 @ t ~ 16.2 (raised‑cosine gate)
- 3Q: F_GHZ > 0.8 band; peak near (J ~ 0.29, lambda_c ~ 0.88)
- 4Q: W4 > 0.82 at (J ~ 0.275, lambda_c ~ 0.87) [reported]
- 5Q: this file — aim GHZ5 > 0.8 with ARP sequential cascade
"""
from __future__ import annotations
import numpy as np
from numpy import kron
from scipy.linalg import expm
import argparse
import itertools
import csv
from typing import List, Dict, Tuple

# ---------- Pauli / tensor helpers ----------
I = np.eye(2, dtype=complex)
X = np.array([[0, 1],[1, 0]], dtype=complex)
Y = np.array([[0, -1j],[1j, 0]], dtype=complex)
Z = np.array([[1, 0],[0,-1]], dtype=complex)


def kronN(ops: List[np.ndarray]) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for op in ops:
        out = kron(out, op)
    return out


def op_on(n: int, which: int, op: np.ndarray) -> np.ndarray:
    ops = [I]*n
    ops[which] = op
    return kronN(ops)


def two_body_xy(n: int, i: int, j: int) -> np.ndarray:
    return 0.5*(op_on(n, i, X) @ op_on(n, j, X) + op_on(n, i, Y) @ op_on(n, j, Y))

# ---------- States & fidelities ----------

def basis_state(n: int, ones: List[int]) -> np.ndarray:
    N = 2**n
    idx = 0
    for q in range(n):
        bit = 1 if q in ones else 0
        idx = (idx << 1) | bit
    v = np.zeros((N,), dtype=complex)
    v[idx] = 1.0
    return v


def GHZ5_state() -> np.ndarray:
    # |GHZ5> = (|00000> + |11111>)/sqrt(2)
    n = 5
    psi = basis_state(n, []) + basis_state(n, [0,1,2,3,4])
    psi = psi / np.linalg.norm(psi)
    return psi


def ket_fidelity(psi: np.ndarray, phi: np.ndarray) -> float:
    return float(np.abs(np.vdot(phi, psi))**2)

# ---------- ARP scaling hook ----------

def arp_gain(t: float, lambda_c: float, params: Dict) -> float:
    T = params.get('segment_T', None)
    if T is None:
        env = 1.0
    else:
        x = np.clip(t/max(T,1e-12), 0.0, 1.0)
        env = 0.5 - 0.5*np.cos(np.pi*x)  # raised‑cosine
    return 1.0 + lambda_c*env

# ---------- Segment Hamiltonians ----------

def seg_H(n: int, edge: Tuple[int,int], J: float, lambda_c: float, phi: float, params: Dict|None=None) -> np.ndarray:
    i, j = edge
    Hij = two_body_xy(n, i, j)
    Hloc = (phi/2.0)*op_on(n, i, Z) + (-phi/2.0)*op_on(n, j, Z)
    J_eff = J * arp_gain(0.0, lambda_c, params or {})
    return J_eff*Hij + Hloc

# ---------- Evolution ----------

def evolve_seg(psi: np.ndarray, H: np.ndarray, tau: float) -> np.ndarray:
    U = expm(-1j*H*tau)
    return U @ psi


def seq_run_5q(J: float, lambda_c: float, phi: float, tau: float,
               sequence: List[Tuple[int,int]]|None=None,
               init: np.ndarray|None=None,
               params: Dict|None=None) -> np.ndarray:
    """Default sequence: AB->BC->CD->DE on 5‑qubit chain indices (0..4)."""
    n = 5
    if init is None:
        init = basis_state(n, [])  # |00000>
    if sequence is None:
        sequence = [(0,1),(1,2),(2,3),(3,4)]
    p = params or {}
    p = {**p, 'segment_T': tau}

    psi = init
    for (i,j) in sequence:
        H = seg_H(n, (i,j), J, lambda_c, phi, params=p)
        psi = evolve_seg(psi, H, tau)
    return psi/np.linalg.norm(psi)

# ---------- GHZ5 fidelity ----------

def compute_GHZ5_fidelity(psi: np.ndarray) -> float:
    return ket_fidelity(psi, GHZ5_state())

# ---------- Sweep ----------

def sweep_5q(J_vals, lc_vals, phi_vals, tau_vals,
             csv_path: str='5q_GHZ_sweep.csv', topk: int=20,
             sequence: List[Tuple[int,int]]|None=None):
    header = ['J','lambda_c','phi','tau','F_GHZ5']
    rows = []
    for J, lc, phi, tau in itertools.product(J_vals, lc_vals, phi_vals, tau_vals):
        psi = seq_run_5q(J, lc, phi, tau, sequence=sequence)
        FG = compute_GHZ5_fidelity(psi)
        rows.append((J, lc, phi, tau, FG))
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    rows.sort(key=lambda r: r[-1], reverse=True)
    return rows[:topk]

# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(description='5Q bell‑chip sequential GHZ5 sweep')
    sub = parser.add_subparsers(dest='cmd', required=False)

    p_run = sub.add_parser('run', help='single run')
    p_run.add_argument('--J', type=float, default=0.28)
    p_run.add_argument('--lambda_c', type=float, default=0.88)
    p_run.add_argument('--phi', type=float, default=0.0)
    p_run.add_argument('--tau', type=float, default=5.4)

    p_sweep = sub.add_parser('sweep', help='parameter sweep -> CSV')
    p_sweep.add_argument('--Jmin', type=float, default=0.27)
    p_sweep.add_argument('--Jmax', type=float, default=0.31)
    p_sweep.add_argument('--Jnum', type=int, default=9)
    p_sweep.add_argument('--LCmin', type=float, default=0.85)
    p_sweep.add_argument('--LCmax', type=float, default=0.92)
    p_sweep.add_argument('--LCnum', type=int, default=8)
    p_sweep.add_argument('--PHInum', type=int, default=16)
    p_sweep.add_argument('--TAUmin', type=float, default=3.0)
    p_sweep.add_argument('--TAUmax', type=float, default=8.0)
    p_sweep.add_argument('--TAUnum', type=int, default=11)
    p_sweep.add_argument('--csv', type=str, default='5q_GHZ_sweep.csv')
    p_sweep.add_argument('--topk', type=int, default=20)

    args = parser.parse_args()

    if args.cmd == 'run':
        psi = seq_run_5q(args.J, args.lambda_c, args.phi, args.tau)
        FG = compute_GHZ5_fidelity(psi)
        print(f"F_GHZ5={FG:.6f}  (J={args.J}, lambda_c={args.lambda_c}, phi={args.phi}, tau={args.tau})")
        return

    # default: sweep
    J_vals  = np.linspace(getattr(args, 'Jmin', 0.27), getattr(args, 'Jmax', 0.31), getattr(args, 'Jnum', 9))
    LC_vals = np.linspace(getattr(args, 'LCmin', 0.85), getattr(args, 'LCmax', 0.92), getattr(args, 'LCnum', 8))
    PHI_vals = np.linspace(0.0, 2*np.pi, getattr(args, 'PHInum', 16), endpoint=False)
    TAU_vals = np.linspace(getattr(args, 'TAUmin', 3.0), getattr(args, 'TAUmax', 8.0), getattr(args, 'TAUnum', 11))

    best = sweep_5q(J_vals, LC_vals, PHI_vals, TAU_vals, csv_path=getattr(args,'csv','5q_GHZ_sweep.csv'), topk=getattr(args,'topk',20))
    print("Top results (J, lambda_c, phi, tau, F_GHZ5):")
    for row in best:
        print(row)

if __name__ == '__main__':
    main()
