# Core Equations & Fits (ARP, Noise, ECC)

## ARP update with noise
\[
\dot{G}_{ij} = \alpha |I_{ij}| - \mu G_{ij} - \eta_{\mathrm{noise}}\,\xi(t)
\]
- \(\alpha\): drive gain; \(\mu\): leakage; \(\eta_{\mathrm{noise}}\xi(t)\): effective stochastic noise.

## Noise → segment probabilities (per qubit)
- Dephasing (phase flip): \(p_\phi = 1 - e^{-\Gamma_2\, \tau}\)
- Amplitude damping: \(p_\mathrm{amp} = 1 - e^{-\eta\, \tau}\)

Implemented via Kraus channels after each segment (per-qubit):
- Dephasing: \(\mathcal{E}(\rho) = (1-p_\phi)\,\rho + p_\phi\, Z\rho Z\)
- Amplitude damping: with Kraus \(K_0 = \begin{bmatrix}1&0\\0&\sqrt{1-p_\mathrm{amp}}\end{bmatrix}\),
  \(K_1 = \begin{bmatrix}0&\sqrt{p_\mathrm{amp}}\\0&0\end{bmatrix}\)

## Linear interaction-time fit (continue‑ramp)
\[
\tau(N) \approx 0.88\,N + 11.5
\]
confirmed from 65,536Q → 2,097,152Q (and beyond), with DD + filtered \(\phi\)-ramp + ECC on.

## Fidelity scaling (representative)
For reported threads:
- With baseline noise (\(\eta=0.005, \Gamma_2=0.008\)) and ECC on, fidelity remains \(F > 10^{-50}\) up to \(N \sim 3.2\times 10^6\).
- Without ECC, exponential decay \(F \approx A e^{-\beta N}\) describes the floor; \(\beta\) dominated by \(\Gamma_2\).

## Targets
- W/GHZ chains under sequential + sandwich, continue‑ramp across passes.
- Alphabet model O(N) cost; MPS backend feasible for very large N.

