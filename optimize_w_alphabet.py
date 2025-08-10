# file: optimize_w_alphabet.py
import numpy as np
from typing import Tuple

def edge_unitary(theta: float, phi: float):
    c, s = np.cos(theta), np.sin(theta)
    e_m = np.exp(-1j*phi); e_p = np.exp(+1j*phi)
    return np.array([[c, -1j*e_m*s],
                     [-1j*e_p*s,  c     ]], dtype=complex)

def w_fidelity_alphabet(
    N: int,
    theta: float,
    phi0: float,
    phi_ramp: float,
    repeat: int = 1,
    sandwich: bool = True,
    ramp_mode: str = "mirror",   # "mirror" or "continue"
) -> float:
    """
    Single-excitation alphabet model with phase ramp options:
      - mirror: reverse pass uses reversed segment index (our previous behavior)
      - continue: reverse pass keeps increasing the ramp index (more dispersion)
    """
    amps = np.zeros(N, dtype=complex); amps[0] = 1.0
    edges = [(i, i+1) for i in range(N-1)]
    passes = []
    for _ in range(max(1, int(repeat))):
        passes.append(edges)
        if sandwich:
            passes.append(list(reversed(edges)))

    k_global = 0
    for path_i, path in enumerate(passes):
        for k_local, (i, j) in enumerate(path):
            if ramp_mode == "mirror":
                k = k_local
            else:  # "continue"
                k = k_global
            phi_k = phi0 + k * phi_ramp
            U = edge_unitary(theta, phi_k)
            v = np.array([amps[i], amps[j]], dtype=complex)
            amps[i], amps[j] = U @ v
            k_global += 1

    overlap = amps.sum() / np.sqrt(N)  # <W|ψ>
    return float(np.abs(overlap)**2)

def optimize_w(
    N: int,
    phi0: float,
    theta_grid: np.ndarray,
    ramp_grid: np.ndarray,
    repeats=(1,2,3),
    ramp_modes=("mirror","continue"),
) -> Tuple[float, dict]:
    bestF, best = -1.0, {}
    for mode in ramp_modes:
        for rep in repeats:
            for theta in theta_grid:
                for rmp in ramp_grid:
                    F = w_fidelity_alphabet(N, theta, phi0, rmp, repeat=rep, sandwich=True, ramp_mode=mode)
                    if F > bestF:
                        bestF = F
                        best = dict(N=N, theta=float(theta), phi0=float(phi0),
                                    phi_ramp=float(rmp), repeat=int(rep),
                                    ramp_mode=mode, F=float(F))
    return bestF, best

if __name__ == "__main__":
    # Example sweep ranges (fast): tune if you want deeper search.
    Ns = [8, 10, 12, 20, 32, 64]
    phi0 = np.pi/8
    theta_grid = np.linspace(0.2, 1.8, 33)          # effective rotation per segment
    ramp_grid  = np.linspace(0.0, 0.12, 25)         # radians per segment

    print("N | best F     | theta  | phi_ramp | repeat | mode")
    print("--+------------+--------+----------+--------+--------")
    for N in Ns:
        F, cfg = optimize_w(N, phi0, theta_grid, ramp_grid)
        print(f"{N:2d}| {F:10.6f} | {cfg['theta']:6.3f} | {cfg['phi_ramp']:8.3f} |"
              f"   {cfg['repeat']}   | {cfg['ramp_mode']}")
