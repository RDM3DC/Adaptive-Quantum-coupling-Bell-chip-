
# Hardware Pathways for the Bell‑chip ARP Model
*A practical appendix for implementing \(g_x, g_y\) and running Predictions P1–P3 on real platforms.*

## 0) Scope & Notation
Working model (qubit frame):
\[
H = \frac{\hbar}{2}\left(\Delta\,\sigma_z + \Omega\,\sigma_x\right) \;+\; \hbar\lambda\,(g_x\,\sigma_x + g_y\,\sigma_y),
\quad
\dot g_x = \tfrac{\alpha}{2}x - \mu g_x,\;\;
\dot g_y = \tfrac{\alpha}{2}y - \mu g_y.
\]
- \(x,y,z\): Bloch components of \(\rho\).
- Control knobs: detuning \(\Delta\), Rabi \(\Omega\), coupling \(\lambda\), write‑in \(\alpha\), forgetting \(\mu\).
- Effective coupling \(\Lambda = \lambda \alpha / \mu\) governs hysteresis and skew strength.

Implementation principle: **compute \(g_{x,y}(t)\) in real time** (fast DSP/FPGA or analog IIR), then **apply as I/Q quadratures** that add to the native drive.

---

## 1) Superconducting qubits (transmons)
**Mapping to hardware**
- \(\sigma_x\): in‑phase (I) drive; \(\sigma_y\): quadrature (Q) drive with \(+90^\circ\) phase.
- Realize \(g_x, g_y\) by modulating I/Q envelopes: \(I(t) \gets I_0(t) + \lambda\,g_x(t)\), \(Q(t) \gets Q_0(t) + \lambda\,g_y(t)\).

**Memory loop (FPGA in the loop)**
Discrete‑time update (sample period \(\Delta t\)):
\[
g_x[n{+}1]=(1{-}\mu\Delta t)\,g_x[n]+\frac{\alpha\Delta t}{2}\,x[n],\quad
g_y[n{+}1]=(1{-}\mu\Delta t)\,g_y[n]+\frac{\alpha\Delta t}{2}\,y[n].
\]
- Obtain \(x[n],y[n]\) from a streaming estimator (e.g., demodulated \(\langle\sigma_{x,y}\rangle\) via weak measurement) **or** from model‑predicted \(x,y\) (closed‑loop sim) for first demonstrations.
- Stability tips: choose \(\mu\Delta t \ll 1\) and \(\alpha\Delta t \lesssim 1\). Start with \(\mu \approx 0.1\,\alpha\) if you target \(\Lambda\sim 1\).

**Prediction protocols**
- **P1 (Hysteresis):** CW drive; sweep \(\Delta\) up then down over a symmetric window. Record steady \(z^\*(\Delta)\). Control: \(\lambda=0\) path to verify no hysteresis without ARP.
- **P2 (Decoherence floor):** Build phase with constant drive (\(\Omega\)) to raise \(y\), then set \(\Omega\to0\) and measure FID on \(y(t)\). Fit \(T_2^{-1}=\Gamma_2+\eta g_y^2\). Since \(g_y\approx(\alpha/2\mu)\,y\) (fast memory), \(T_2^{-1}\propto y^2\).
- **P3 (Line‑shape skew):** Drive \(\Omega(t)=A\cos\omega t\), lock‑in detect \(z(t)\) at \(\omega\). Compare \(\lambda=0\) vs \(\lambda>0\) to reveal odd‑in‑\(\Delta\) skew.

**Diagnostics & pitfalls**
- IQ imbalance / mixer leakage can fake \(g_y\) → run \(\lambda=0\) control and calibrate quadrature balance.
- LO drift alters effective \(\Delta\) during sweeps → log Δ in the control frame and/or interleave up/down points.
- Ensure latency (estimator→I/Q) \(\ll\) Rabi period; otherwise use the **slow‑memory** model in analysis.

---

## 2) Trapped ions
**Mapping to hardware**
- \(\sigma_x,\sigma_y\) via Raman beatnote phase; in‑phase/quadrature control uses DDS phase advances.
- Implement \(g_{x,y}\) by updating amplitude/phase on the DDS each sample using the same IIR rules above.

**Practical notes**
- Use detuning around the carrier; avoid motional sideband crosstalk for P1–P3 by staying in the Lamb‑Dicke regime.
- Latency is typically small; still validate that the controller update rate \(\gtrsim 10\times\) the drive frequency.

---

## 3) Integrated photonics (dual‑rail / MZI)
**Mapping to hardware**
- Dual‑rail qubit: \(\sigma_x\) is balanced mixing (MZI), \(\sigma_y\) adds a \(\pi/2\) phase between rails.
- Implement \(g_x\) via differential phase shifters; \(g_y\) via an added phase arm (or EO phase shifter at \(\pi/2\)).
- Memory loop realized in control electronics that drive the phase shifters (EO preferred over thermal for speed).

**Readout**
- Convert rail intensities to Bloch \(z\); retrieve \(x,y\) using calibrated phase dithers (homodyne‑style).

---

## 4) Parameter targets & units
- Normalize by your reference \(\Omega\): report \(\Delta/\Omega\), \(\Gamma_2/\Omega\), and \(\Lambda=\lambda\alpha/\mu\).
- A good starting regime for visibility: \(\Lambda \in [1,3]\), sweep center \(\Delta_0 \approx 0\).

---

## 5) Data & file formats
- **P1:** `Delta,z_up,z_dn` (same grid).  
- **P2:** `y,T2` *or* `y2,T2` (units noted).  
- **P3:** `Delta,R` for both \(\lambda=0\) and \(\lambda>0\) runs.  
Use the `templates/` CSVs and the units‑aware wrapper to normalize.

---

## 6) Repro pipeline (what to run)
1. Generate predictions: `bell_chip_predictions.py --which ALL --out ./outputs`  
2. Phase map (choose region): `parameter_map_hysteresis.py`  
3. Lab playbook quick check: `bell_chip_lab_playbook.py`  
4. Fit data (units aware): `fit_bell_chip_units.py ...` → `fit_results/fit_params.json`

---

## 7) Controls & ablations
- \(\lambda=0\) (disable ARP) — removes hysteresis & skew; baseline for all figures.  
- \(\alpha=0\) (no write‑in) — memory frozen; probes imprint mechanism.  
- Large \(\mu\) (fast forgetting) — suppresses effects as \(g\to 0\).

---

## 8) Risk log (things that mimic ARP)
- AM‑to‑PM conversion in mixers → fake \(\sigma_y\) terms.  
- Slow thermal drift in photonics → apparent hysteresis on \(\Delta\).  
- Nonlinear readout saturation → skewed line‑shapes.  
Mitigation: interleave controls, calibrate IQ, restrict drive range, use \(\lambda=0\) baselines.
