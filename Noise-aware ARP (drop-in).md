# per-segment noise: phase-flip (Γ₂) + amplitude damping (η)
def apply_noise(rho, n, p_phi, p_amp):
    I2 = np.eye(2); Z = np.array([[1,0],[0,-1]], complex)
    K0_phi, K1_phi = np.sqrt(1-p_phi)*I2, np.sqrt(p_phi)*Z
    K0_amp = np.array([[1,0],[0,np.sqrt(1-p_amp)]], complex)
    K1_amp = np.array([[0,np.sqrt(p_amp)],[0,0]], complex)
    def lift(K,q):
        ops=[I2]*n; ops[q]=K
        M=ops[0]
        for k in ops[1:]: M=np.kron(M,k)
        return M
    for q in range(n):  # dephasing on each qubit
        L0,L1 = lift(K0_phi,q), lift(K1_phi,q)
        rho = L0@rho@L0.conj().T + L1@rho@L1.conj().T
    for q in range(n):  # amplitude damping on each qubit
        L0,L1 = lift(K0_amp,q), lift(K1_amp,q)
        rho = L0@rho@L0.conj().T + L1@rho@L1.conj().T
    return rho

# usage inside each segment (duration = τ_seg)
# p_phi = 1 - exp(-Γ₂ * τ_seg)      # dephasing
# p_amp = 1 - exp(-η *  τ_seg)      # amplitude damping
