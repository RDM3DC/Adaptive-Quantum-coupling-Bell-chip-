#!/usr/bin/env python3
"""
fit_bell_chip_units.py
----------------------
Units-aware wrapper for fit_bell_chip.py.

- Accepts Δ, Ω, Γ2 units (GHz cycles/s, rad/s, or normalized by Ω).
- Normalizes everything to the model's dimensionless frame with Ω_model = 1.
- Calls `fit_bell_chip.py` on the normalized CSVs.

Example:
  python fit_bell_chip_units.py \
    --p1 P1_curves.csv --p2 P2_T2_summary.csv \
    --p3-no P3_noARP.csv --p3-arp P3_ARP.csv \
    --delta-unit GHz --omega-ref 0.8 --gamma2-unit us_inv \
    --out ./fit_results

Supported unit flags:
  --delta-unit    : GHz | rad_s | norm   (norm means already Δ/Ω)
  --gamma2-unit   : GHz | rad_s | us_inv | ns_inv | norm  (norm means already Γ2/Ω)
  --omega-ref     : reference Ω in GHz (cycles/s). If provided and delta-unit/gamma2-unit in GHz or rad_s,
                    we normalize by Ω_ref (converted to rad/s) so that Ω_model=1.
"""
import argparse, os, csv, numpy as np, pandas as pd, subprocess, sys, tempfile, math

def to_angular_from_ghz(x_ghz):
    # GHz cycles/s -> rad/s
    return x_ghz * 2.0 * math.pi * 1e9

def inv_time_to_rad_over_s(x_inv, unit):
    # Convert inverse time units to rad/s if possible
    if unit == "us_inv":
        return x_inv * 1e6 * 2.0 * math.pi   # 1/us -> 1e6 /s cycles -> *2π rad/s
    if unit == "ns_inv":
        return x_inv * 1e9 * 2.0 * math.pi
    return None

def normalize_delta_array(delta, unit, omega_ref_ghz=None):
    if unit == "norm":
        return delta  # already Δ/Ω
    if omega_ref_ghz is None:
        raise ValueError("omega-ref is required when delta-unit is GHz or rad_s")
    omega_ref_rad = to_angular_from_ghz(omega_ref_ghz)
    if unit == "GHz":
        # Δ(Ω) in GHz cycles/s -> rad/s -> divide by Ω_ref(rad/s)
        delta_rad = to_angular_from_ghz(delta)
        return delta_rad / omega_ref_rad
    if unit == "rad_s":
        return delta / omega_ref_rad
    raise ValueError("Unsupported delta-unit")

def normalize_gamma2_array(g2, unit, omega_ref_ghz=None):
    if unit == "norm":
        return g2
    if unit == "GHz":
        if omega_ref_ghz is None:
            raise ValueError("omega-ref is required when gamma2-unit is GHz")
        omega_ref_rad = to_angular_from_ghz(omega_ref_ghz)
        g2_rad = to_angular_from_ghz(g2)
        return g2_rad / omega_ref_rad
    if unit == "rad_s":
        if omega_ref_ghz is None:
            raise ValueError("omega-ref is required when gamma2-unit is rad_s")
        omega_ref_rad = to_angular_from_ghz(omega_ref_ghz)
        return g2 / omega_ref_rad
    if unit in ("us_inv", "ns_inv"):
        if omega_ref_ghz is None:
            raise ValueError("omega-ref is required when gamma2-unit is inverse time")
        omega_ref_rad = to_angular_from_ghz(omega_ref_ghz)
        g2_rad = inv_time_to_rad_over_s(g2, unit)
        return g2_rad / omega_ref_rad
    raise ValueError("Unsupported gamma2-unit")

def normalize_p1_csv(path, delta_unit, omega_ref_ghz, out_path):
    df = pd.read_csv(path)
    df["Delta"] = normalize_delta_array(df["Delta"].to_numpy(), delta_unit, omega_ref_ghz)
    df.to_csv(out_path, index=False)

def normalize_p2_csv(path, gamma2_unit, omega_ref_ghz, out_path):
    df = pd.read_csv(path)
    # If T2 present, compute T2_inv first (in 1/time), then normalize
    if "T2_inv" in df.columns:
        g2_raw = df["T2_inv"].to_numpy()
    elif "T2" in df.columns:
        # units of T2 must be provided via gamma2-unit (us_inv or ns_inv)
        g2_raw = 1.0 / df["T2"].to_numpy()
    else:
        raise ValueError("P2 CSV must contain T2 or T2_inv")
    g2_norm = normalize_gamma2_array(g2_raw, gamma2_unit, omega_ref_ghz)
    df["T2_inv"] = g2_norm
    df.to_csv(out_path, index=False)

def normalize_p3_csv(path, delta_unit, omega_ref_ghz, out_path):
    df = pd.read_csv(path)
    df["Delta"] = normalize_delta_array(df["Delta"].to_numpy(), delta_unit, omega_ref_ghz)
    df.to_csv(out_path, index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p1", help="P1 CSV (Delta,z_up,z_dn)")
    ap.add_argument("--p2", help="P2 CSV (y2|y|A and T2|T2_inv)")
    ap.add_argument("--p3-no", dest="p3_no", help="P3 (no ARP) CSV (Delta,R)")
    ap.add_argument("--p3-arp", dest="p3_arp", help="P3 (with ARP) CSV (Delta,R)")
    ap.add_argument("--delta-unit", choices=["GHz","rad_s","norm"], default="norm")
    ap.add_argument("--gamma2-unit", choices=["GHz","rad_s","us_inv","ns_inv","norm"], default="norm")
    ap.add_argument("--omega-ref", type=float, help="Reference Ω in GHz (cycles/s) for normalization")
    ap.add_argument("--out", default="./fit_results")
    ap.add_argument("--alpha", type=float, default=0.6)
    ap.add_argument("--mu", type=float, default=0.2)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    tmp = os.path.join(args.out, "_normalized"); os.makedirs(tmp, exist_ok=True)

    # Normalize provided CSVs
    p1n = p2n = p3n = p3an = None
    if args.p1:
        p1n = os.path.join(tmp, "P1_curves_normalized.csv")
        normalize_p1_csv(args.p1, args.delta_unit, args.omega_ref, p1n)
    if args.p2:
        p2n = os.path.join(tmp, "P2_T2_summary_normalized.csv")
        normalize_p2_csv(args.p2, args.gamma2_unit, args.omega_ref, p2n)
    if args.p3_no:
        p3n = os.path.join(tmp, "P3_noARP_normalized.csv")
        normalize_p3_csv(args.p3_no, args.delta_unit, args.omega_ref, p3n)
    if args.p3_arp:
        p3an = os.path.join(tmp, "P3_ARP_normalized.csv")
        normalize_p3_csv(args.p3_arp, args.delta_unit, args.omega_ref, p3an)

    # Call base fitter (assumes fit_bell_chip.py present in PATH or same folder)
    fitter = os.path.join(os.path.dirname(__file__), "fit_bell_chip.py")
    if not os.path.exists(fitter):
        print("WARNING: fit_bell_chip.py not found alongside this script. Ensure it's accessible.")
        fitter = "fit_bell_chip.py"

    cmd = ["python", fitter, "--out", args.out, "--alpha", str(args.alpha), "--mu", str(args.mu)]
    if p1n: cmd += ["--p1", p1n]
    if p2n: cmd += ["--p2", p2n]
    if p3n and p3an: cmd += ["--p3-no", p3n, "--p3-arp", p3an]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
