# Reproducing the paper figures

This document lists the exact CLI runs used to produce the main CSVs and plots. Replace script names if your local filenames differ.

## Environment
```bash
python --version    # >= 3.9
pip install numpy scipy matplotlib

Self-check

python bell_chip_eight_qubit.py selftest

Representative big‑N runs (noise + DD + filtered φ‑ramp + ECC)

Parameters match the public thread: η=0.005, Γ2=0.008

65,536Q

python bell_chip_eight_qubit.py run \
  --J 0.185 --lambda_c 0.725 --phi 0.2617993878 --tau 57640 \
  --phi-ramp 0.06 --sandwich --repeat 3 \
  --dephase 0.008 --ampdamp 0.005 \
  --dd on --ecc on

131,072Q

python bell_chip_eight_qubit.py run \
  --J 0.185 --lambda_c 0.725 --phi 0.2617993878 --tau 115280 \
  --phi-ramp 0.06 --sandwich --repeat 3 \
  --dephase 0.008 --ampdamp 0.005 \
  --dd on --ecc on

262,144Q

python bell_chip_eight_qubit.py run \
  --J 0.185 --lambda_c 0.725 --phi 0.2617993878 --tau 230560 \
  --phi-ramp 0.06 --sandwich --repeat 3 \
  --dephase 0.008 --ampdamp 0.005 \
  --dd on --ecc on

524,288Q

python bell_chip_eight_qubit.py run \
  --J 0.185 --lambda_c 0.725 --phi 0.2617993878 --tau 461385 \
  --phi-ramp 0.06 --sandwich --repeat 3 \
  --dephase 0.008 --ampdamp 0.005 \
  --dd on --ecc on

1,048,576Q

python bell_chip_eight_qubit.py run \
  --J 0.185 --lambda_c 0.725 --phi 0.2617993878 --tau 922758 \
  --phi-ramp 0.06 --sandwich --repeat 3 \
  --dephase 0.008 --ampdamp 0.005 \
  --dd on --ecc on

2,097,152Q

python bell_chip_eight_qubit.py run \
  --J 0.185 --lambda_c 0.725 --phi 0.2617993878 --tau 1845505 \
  --phi-ramp 0.06 --sandwich --repeat 3 \
  --dephase 0.008 --ampdamp 0.005 \
  --dd on --ecc on

τ vs N quick check

Follows the paper’s linear fit: τ ≈ 0.88 N + 11.5

Notes
•If --dd/--ecc are not recognized flags in your version, annotate that these were effective toggles in the thread; provide your closest equivalents.
•CSV outputs and plots can be gathered with topk_tools.py and your plotting script; see README.md.
