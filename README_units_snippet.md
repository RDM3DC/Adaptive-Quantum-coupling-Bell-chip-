
## Units & Conventions

- Labs often report detuning Δ and drive Ω in **GHz (cycles/s)**. Our model integrates in a dimensionless frame with **Ω=1**.
- To map data → model, we normalize by the reference drive Ω:
  $$\Delta_{\text{norm}} = \frac{2\pi\,\Delta_{\text{GHz}}}{2\pi\,\Omega_{\text{GHz}}} = \frac{\Delta}{\Omega}, \quad
    \Gamma_{2,\text{norm}} = \frac{\Gamma_{2,\text{angular}}}{2\pi\,\Omega_{\text{GHz}}}.$$
- Use the wrapper below to convert CSVs automatically.

### Units-aware fitting
```bash
python fit_bell_chip_units.py   --p1 P1_curves.csv   --p2 P2_T2_summary.csv   --p3-no P3_noARP.csv   --p3-arp P3_ARP.csv   --delta-unit GHz   --gamma2-unit us_inv   --omega-ref 0.8 \  # Ω_ref in GHz (cycles/s)
  --out ./fit_results
```
- Supported `--delta-unit`: `GHz`, `rad_s`, `norm` (already Δ/Ω)
- Supported `--gamma2-unit`: `GHz`, `rad_s`, `us_inv`, `ns_inv`, `norm`

Templates are in `templates/`:
- `P1_curves_template.csv` (Δ, z_up, z_dn)
- `P2_T2_summary_template.csv` (A, T2)
- `P3_noARP_template.csv`, `P3_ARP_template.csv` (Δ, R)
