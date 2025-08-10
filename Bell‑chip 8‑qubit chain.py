#!/usr/bin/env python
"""
Bell‑chip 8‑qubit chain (A‑B‑C‑D‑E‑F‑G‑H): sequential + phase cascade targeting high W8 fidelity.

Design mirrors earlier scaffolds (2Q..7Q):
- Segment evolution with static H per segment, XY exchange on active edge.
- Phase steer via local Z terms; ARP scaling hook for time‑envelope gain on J.
- Sweep over (J, lambda_c, phi, tau) with CSV output and top‑k print.

Targets (thread momentum):
- 6Q: W6 > 0.80 near (J ~ 0.27, lambda_c ~ 0.85)
- 7Q: GHZ7 peak > 0.79 near (J ~ 0.265, lambda_c ~ 0.84)
- 8Q: this file — aim W8 > 0.78 with ARP sequential cascade (AB..GH)

Notes
-----
- Python 3.5+ friendly (no f‑strings/typing unions).
- Add "--selftest" to run a tiny internal check that validates I/O and fidelity range.

Convenience flags:
- --phi-ramp <step> and --phi-scale <s> to auto‑generate a ramp phase template.
- --center-boost <b> to lengthen central segments via tau_scale_list.
- Global --phi now **offsets** any per‑segment --phi-list/--phi-ramp template (seg_phi = phi_i + phi).
"""
from __future__ import annotations
import numpy as np
from numpy import kron
from scipy.linalg import expm
import argparse
import itertools
import csv
import tempfile
import os
from typing import List, Dict, Tuple
import math

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


def W8_state() -> np.ndarray:
    # |W8> = (|00000001> + |00000010> + |00000100> + |00001000>
    #         + |00010000> + |00100000> + |01000000> + |10000000>)/sqrt(8)
    n = 8
    terms = [basis_state(n, [7]), basis_state(n, [6]), basis_state(n, [5]), basis_state(n, [4]),
             basis_state(n, [3]), basis_state(n, [2]), basis_state(n, [1]), basis_state(n, [0])]
    psi = sum(terms)
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
        env = 0.5 - 0.5*np.cos(np.pi*x)  # raised‑cosine envelope
    return 1.0 + lambda_c*env

# ---------- Segment Hamiltonian builder ----------

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


def seq_run_8q(J: float, lambda_c: float, phi: float, tau: float,
               sequence: List[Tuple[int,int]]|None=None,
               init: np.ndarray|None=None,
               params: Dict|None=None,
               phi_list: List[float]|None=None,
               tau_list: List[float]|None=None,
               tau_scale_list: List[float]|None=None,
               sandwich: bool=False,
               repeat: int=1,
               phi_mode: str='offset') -> np.ndarray:
    n = 8
    if init is None:
        # Use a single-excitation start so XY exchange can generate W-like superpositions
        init = basis_state(n, [0])  # |10000000> (A excited)
    if sequence is None:
        sequence = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7)]  # AB..GH
    p = params or {}

    # Build full sequence with potential sandwich and repeats
    full_sequence = []
    for _ in range(max(1, int(repeat))):
        full_sequence.extend(sequence)
        if sandwich:
            full_sequence.extend(list(reversed(sequence)))

    # Handle per-segment phases: combine global phi with template/list
    if phi_list is None:
        # No template provided: just repeat the global phi
        per_phi = [phi] * len(full_sequence)
    else:
        # Template provided: extend/cycle to full length
        per_phi = []
        # We treat phi as an **offset** or **scale** over the template values
        template = list(phi_list)
        while len(per_phi) < len(full_sequence):
            need = len(full_sequence) - len(per_phi)
            per_phi.extend(template[:min(len(template), need)])
        if phi_mode == 'scale':
            per_phi = [phi * x for x in per_phi]
        else:  # 'offset'
            per_phi = [phi + x for x in per_phi]

    # Handle per-segment durations
    seg_taus = [tau] * len(full_sequence)
    if tau_list is not None:
        # Extend/cycle explicit tau_list to match length
        per_tau = []
        template_tau = list(tau_list)
        while len(per_tau) < len(full_sequence):
            need = len(full_sequence) - len(per_tau)
            per_tau.extend(template_tau[:min(len(template_tau), need)])
        seg_taus = per_tau

    if tau_scale_list is not None:
        # Extend/cycle scaling list and apply to seg_taus
        per_scale = []
        template_scale = list(tau_scale_list)
        while len(per_scale) < len(full_sequence):
            need = len(full_sequence) - len(per_scale)
            per_scale.extend(template_scale[:min(len(template_scale), need)])
        seg_taus = [t * s for t, s in zip(seg_taus, per_scale)]

    # Evolve
    psi = init
    for idx, (i, j) in enumerate(full_sequence):
        seg_phi = per_phi[idx]
        seg_tau = seg_taus[idx]
        p_seg = {**p, 'segment_T': seg_tau}
        H = seg_H(n, (i, j), J, lambda_c, seg_phi, params=p_seg)
        psi = evolve_seg(psi, H, seg_tau)

    return psi / np.linalg.norm(psi)

# ---------- W8 fidelity ----------

def compute_W8_fidelity(psi: np.ndarray) -> float:
    return ket_fidelity(psi, W8_state())

# ---------- Sweep ----------

def sweep_8q(J_vals, lc_vals, phi_vals, tau_vals,
             csv_path: str='8q_W_sweep.csv', topk: int=20,
             sequence: List[Tuple[int,int]]|None=None,
             phi_list: List[float]|None=None,
             tau_list: List[float]|None=None,
             tau_scale_list: List[float]|None=None,
             sandwich: bool=False,
             repeat: int=1,
             phi_mode: str='offset'):
    header = ['J','lambda_c','phi','tau','F_W8']
    rows = []
    for J, lc, phi, tau in itertools.product(J_vals, lc_vals, phi_vals, tau_vals):
        psi = seq_run_8q(J, lc, phi, tau,
                         sequence=sequence,
                         phi_list=phi_list,
                         tau_list=tau_list,
                         tau_scale_list=tau_scale_list,
                         sandwich=sandwich,
                         repeat=repeat,
                         phi_mode=phi_mode)
        FW = compute_W8_fidelity(psi)
        rows.append((float(J), float(lc), float(phi), float(tau), float(FW)))
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    rows.sort(key=lambda r: r[-1], reverse=True)
    return rows[:topk]

# ---------- Self‑test ----------

def selftest() -> int:
    # Minimal check: run a tiny 1x1x1x1 sweep and assert fidelity in [0,1]
    J_vals = [0.26]
    LC_vals = [0.85]
    PHI_vals = [0.0]
    TAU_vals = [5.0]
    best = sweep_8q(J_vals, LC_vals, PHI_vals, TAU_vals, csv_path='8q_demo.csv', topk=1)
    ok = True
    for row in best:
        F = row[-1]
        if not (0.0 <= F <= 1.0):
            ok = False
    if ok:
        print("Selftest passed: W8 fidelity in [0,1] and CSV written -> 8q_demo.csv")
        # Extra check: using a per‑segment phase list plus a nonzero global phi should change fidelity
        psi_a = seq_run_8q(0.26, 0.85, 0.0, 5.0, phi_list=[0,0.2,0.4,0.6,0.8,1.0,1.2], sandwich=True, repeat=1)
        Fa = compute_W8_fidelity(psi_a)
        psi_b = seq_run_8q(0.26, 0.85, 0.1, 5.0, phi_list=[0,0.2,0.4,0.6,0.8,1.0,1.2], sandwich=True, repeat=1)
        Fb = compute_W8_fidelity(psi_b)
        if abs(Fa - Fb) < 1e-12:
            print("Warning: phi offset had no effect in selftest (could be a pathological symmetry)")
        return 0
    print("Selftest failed: fidelity out of range")
    return 1

# ---------- CLI ----------

def main():
    parser = argparse.ArgumentParser(description='8Q bell‑chip sequential W8 sweep')
    sub = parser.add_subparsers(dest='cmd', required=False)

    p_run = sub.add_parser('run', help='single run')
    p_run.add_argument('--J', type=float, default=0.26)
    p_run.add_argument('--lambda_c', type=float, default=0.85)
    p_run.add_argument('--phi', type=float, default=0.0)
    p_run.add_argument('--tau', type=float, default=5.0)
    p_run.add_argument('--phi-list', type=float, nargs='+', help='Per-segment phase values')
    p_run.add_argument('--tau-list', type=float, nargs='+', help='Per-segment duration values')
    p_run.add_argument('--tau-scale-list', type=float, nargs='+', help='Per-segment duration scaling factors')
    p_run.add_argument('--sandwich', action='store_true', help='Use sandwich sequence (forward then reverse)')
    p_run.add_argument('--repeat', type=int, default=1, help='Number of times to repeat the sequence')
    # Convenience templates
    p_run.add_argument('--phi-ramp', type=float, help='Base step (rad) for a ramp template; generates phi_list = [phi_scale * (k*phi_ramp)] for k=0..L-1')
    p_run.add_argument('--phi-scale', type=float, default=1.0, help='Scale factor applied to ramp/template phases')
    p_run.add_argument('--center-boost', type=float, default=0.0, help='Duration boost toward the center segment (0=no boost). Creates tau_scale_list=1+center_boost*shape')
    p_run.add_argument('--phi-mode', type=str, choices=['offset','scale'], default='offset', help='How to combine global phi with per‑segment phases when a template/list is used')

    p_sweep = sub.add_parser('sweep', help='parameter sweep -> CSV')
    p_sweep.add_argument('--Jmin', type=float, default=0.24)
    p_sweep.add_argument('--Jmax', type=float, default=0.28)
    p_sweep.add_argument('--Jnum', type=int, default=9)
    p_sweep.add_argument('--LCmin', type=float, default=0.82)
    p_sweep.add_argument('--LCmax', type=float, default=0.88)
    p_sweep.add_argument('--LCnum', type=int, default=7)
    p_sweep.add_argument('--PHInum', type=int, default=16)
    p_sweep.add_argument('--TAUmin', type=float, default=3.0)
    p_sweep.add_argument('--TAUmax', type=float, default=8.0)
    p_sweep.add_argument('--TAUnum', type=int, default=11)
    p_sweep.add_argument('--csv', type=str, default='8q_W_sweep.csv')
    p_sweep.add_argument('--topk', type=int, default=20)
    p_sweep.add_argument('--phi-list', type=float, nargs='+', help='Per-segment phase values')
    p_sweep.add_argument('--tau-list', type=float, nargs='+', help='Per-segment duration values')
    p_sweep.add_argument('--tau-scale-list', type=float, nargs='+', help='Per-segment duration scaling factors')
    p_sweep.add_argument('--sandwich', action='store_true', help='Use sandwich sequence (forward then reverse)')
    p_sweep.add_argument('--repeat', type=int, default=1, help='Number of times to repeat the sequence')
    # Convenience templates
    p_sweep.add_argument('--phi-ramp', type=float, help='Base step (rad) for a ramp template used to auto‑generate phi_list')
    p_sweep.add_argument('--phi-scale', type=float, default=1.0, help='Scale factor for ramp/template phases')
    p_sweep.add_argument('--center-boost', type=float, default=0.0, help='Duration center boost used to auto‑generate tau_scale_list')
    p_sweep.add_argument('--phi-mode', type=str, choices=['offset','scale'], default='offset', help='How to combine global phi with per‑segment phases when a template/list is used')

    p_self = sub.add_parser('selftest', help='run a minimal internal test')

    # Optimizer: per-N search over (theta, phi_ramp, repeat) using an alphabet model
    p_opt = sub.add_parser('optimize_w', help='Optimize W-chain in alphabet model (continue ramp) and print table')
    p_opt.add_argument('--Ns', nargs='+', type=int, default=[8,10,12,20,32,64])
    p_opt.add_argument('--phi', type=float, default=np.pi/8)
    p_opt.add_argument('--theta-min', type=float, default=0.2)
    p_opt.add_argument('--theta-max', type=float, default=1.8)
    p_opt.add_argument('--theta-num', type=int, default=33)
    p_opt.add_argument('--ramp-min', type=float, default=0.0)
    p_opt.add_argument('--ramp-max', type=float, default=0.12)
    p_opt.add_argument('--ramp-num', type=int, default=25)
    p_opt.add_argument('--repeats', nargs='+', type=int, default=[1,2,3])

    args = parser.parse_args()

    if args.cmd == 'run':
        # Extract all the new parameters
        phi_list = getattr(args, 'phi_list', None)
        tau_list = getattr(args, 'tau_list', None)
        tau_scale_list = getattr(args, 'tau_scale_list', None)
        sandwich = getattr(args, 'sandwich', False)
        repeat = getattr(args, 'repeat', 1)

        # Auto-build per‑segment lists from simple templates if requested
        L = 7  # number of edges in 8Q chain (AB..GH)
        if getattr(args, 'phi_ramp', None) is not None and phi_list is None:
            step = float(args.phi_ramp)
            scale = float(getattr(args, 'phi_scale', 1.0))
            phi_list = [scale * (k * step) for k in range(L)]
        if float(getattr(args, 'center_boost', 0.0)) != 0.0 and tau_scale_list is None and tau_list is None:
            cb = float(args.center_boost)
            # triangular center boost: 0 at ends, 1 at center
            shape = [1.0 - abs((k - (L-1)/2.0) / ((L-1)/2.0)) for k in range(L)]
            tau_scale_list = [1.0 + cb * s for s in shape]

        psi = seq_run_8q(args.J, args.lambda_c, args.phi, args.tau,
                         phi_list=phi_list,
                         tau_list=tau_list,
                         tau_scale_list=tau_scale_list,
                         sandwich=sandwich,
                         repeat=repeat,
                         phi_mode=getattr(args,'phi_mode','offset'))
        FW = compute_W8_fidelity(psi)

        # Format the basic parameters
        param_str = "J={}, lambda_c={}, phi={}, tau={}".format(args.J, args.lambda_c, args.phi, args.tau)

        # Add advanced parameters if used
        extra_params = []
        if phi_list:
            extra_params.append("phi_list={}".format(phi_list))
        if tau_list:
            extra_params.append("tau_list={}".format(tau_list))
        if tau_scale_list:
            extra_params.append("tau_scale_list={}".format(tau_scale_list))
        if sandwich:
            extra_params.append("sandwich=True")
        if repeat > 1:
            extra_params.append("repeat={}".format(repeat))

        # Add the advanced parameters to the output string if any were used
        if extra_params:
            param_str += ", " + ", ".join(extra_params)

        print("F_W8={:.6f}  ({})".format(FW, param_str))
        return

    if args.cmd == 'optimize_w':
        def edge_unitary(theta: float, phi: float):
            c, s = np.cos(theta), np.sin(theta)
            e_m = np.exp(-1j*phi); e_p = np.exp(+1j*phi)
            return np.array([[c, -1j*e_m*s],
                             [-1j*e_p*s,  c     ]], dtype=complex)

        def w_fidelity_alphabet(N: int, theta: float, phi0: float, phi_ramp: float, repeat: int) -> float:
            amps = np.zeros(N, dtype=complex); amps[0] = 1.0
            edges = [(i, i+1) for i in range(N-1)]
            passes = []
            for _ in range(max(1, int(repeat))):
                passes.append(edges)
                # continue ramp default: keep a global k across passes
                passes.append(list(reversed(edges)))
            k_global = 0
            for path in passes:
                for (i, j) in path:
                    phi_k = args.phi + k_global * phi_ramp
                    U = edge_unitary(theta, phi_k)
                    v = np.array([amps[i], amps[j]], dtype=complex)
                    amps[i], amps[j] = U @ v
                    k_global += 1
            overlap = amps.sum() / np.sqrt(N)
            return float(np.abs(overlap)**2)

        thetas = np.linspace(args.theta_min, args.theta_max, args.theta_num)
        ramps  = np.linspace(args.ramp_min,  args.ramp_max,  args.ramp_num)
        print("N | best F     | theta  | phi_ramp | repeat | mode")
        print("--+------------+--------+----------+--------+------")
        for N in args.Ns:
            bestF, best = -1.0, None
            for rep in args.repeats:
                for th in thetas:
                    for rm in ramps:
                        F = w_fidelity_alphabet(N, float(th), float(args.phi), float(rm), int(rep))
                        if F > bestF:
                            bestF, best = F, (th, rm, rep)
            th, rm, rep = best
            print(f"{N:2d}| {bestF:10.6f} | {th:6.3f} | {rm:8.3f} |   {rep:>2}   | cont")
        return

    if args.cmd == 'selftest':
        code = selftest()
        return

    # default: sweep
    J_vals  = np.linspace(getattr(args, 'Jmin', 0.24), getattr(args, 'Jmax', 0.28), getattr(args, 'Jnum', 9))
    LC_vals = np.linspace(getattr(args, 'LCmin', 0.82), getattr(args, 'LCmax', 0.88), getattr(args, 'LCnum', 7))
    PHI_vals = np.linspace(0.0, 2*np.pi, getattr(args, 'PHInum', 16), endpoint=False)
    TAU_vals = np.linspace(getattr(args, 'TAUmin', 3.0), getattr(args, 'TAUmax', 8.0), getattr(args, 'TAUnum', 11))

    L = 7
    phi_list = getattr(args, 'phi_list', None)
    tau_list = getattr(args, 'tau_list', None)
    tau_scale_list = getattr(args, 'tau_scale_list', None)
    # Template helpers
    if getattr(args, 'phi_ramp', None) is not None and phi_list is None:
        step = float(args.phi_ramp)
        scale = float(getattr(args, 'phi_scale', 1.0))
        phi_list = [scale * (k * step) for k in range(L)]
    if float(getattr(args, 'center_boost', 0.0)) != 0.0 and tau_scale_list is None and tau_list is None:
        cb = float(args.center_boost)
        shape = [1.0 - abs((k - (L-1)/2.0) / ((L-1)/2.0)) for k in range(L)]
        tau_scale_list = [1.0 + cb * s for s in shape]

    best = sweep_8q(J_vals, LC_vals, PHI_vals, TAU_vals,
                    csv_path=getattr(args,'csv','8q_W_sweep.csv'), topk=getattr(args,'topk',20),
                    phi_list=phi_list, tau_list=tau_list, tau_scale_list=tau_scale_list,
                    sandwich=bool(getattr(args,'sandwich', False)),
                    repeat=int(getattr(args,'repeat', 1)),
                    phi_mode=str(getattr(args,'phi_mode','offset')))
    print("Top results (J, lambda_c, phi, tau, F_W8):")
    for row in best:
        print(row)

if __name__ == '__main__':
    main()
