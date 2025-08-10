That script simulates two qubits with:
H=\tfrac12(\Delta_A\sigma_z^A+\Omega_A\sigma_x^A)+\tfrac12(\Delta_B\sigma_z^B+\Omega_B\sigma_x^B)
+J\,\sigma_x^A\!\otimes\!\sigma_x^B
+\lambda(\,g_x^A\sigma_x^A+g_y^A\sigma_y^A+g_x^B\sigma_x^B+g_y^B\sigma_y^B\,)
+\lambda_c\,g_c\,\sigma_x^A\!\otimes\!\sigma_x^B,
with memory dynamics
\dot g_x^{A,B}=\tfrac{\alpha}{2}\langle\sigma_x^{A,B}\rangle-\mu g_x^{A,B},\quad
\dot g_y^{A,B}=\tfrac{\alpha}{2}\langle\sigma_y^{A,B}\rangle-\mu g_y^{A,B},\quad
\dot g_c=\alpha_c\langle\sigma_x^A\!\otimes\!\sigma_x^B\rangle-\mu_c g_c.

It computes concurrence C(t) (and has fidelity to |\Phi^+\rangle under the hood) and compares with ARP vs no ARP.

Suggested parameter set (normalized units, \Omega=1)
	•	Local: \alpha=0.6,\ \mu=0.2,\ \lambda=0.8  → \Lambda=\lambda\alpha/\mu\approx2.4
	•	Shared: \alpha_c=0.4,\ \mu_c=0.15,\ \lambda_c=0.6
	•	Coupling: J=0.15 (weak XX)
	•	Drive/detuning: \Omega_A=\Omega_B=1,\ \Delta_A=\Delta_B=0
	•	Init state: |00\rangle

To map IBM-style GHz (cycles/s) to the model: multiply Δ, Ω, Γ₂ by 2\pi to get angular, then normalize by your chosen \Omega so \Omega_{\text{model}}=1.

What this shows

The included plot compares concurrence with and without ARP shared memory. You’ll see memory reshapes the entangling dynamics, creating history-dependent entanglement plateaus/oscillations that aren’t there with \lambda=\lambda_c=0.
