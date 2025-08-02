> **ℹ️ Note:** This implementation has been **officially [merged](https://github.com/Qiskit/qiskit/pull/13975) into [Qiskit](https://github.com/Qiskit/qiskit)**.
>
> **Quick usage:**
```python
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import HalfAdderGate

qc = QuantumCircuit(7)
qc.append(HalfAdderGate(3), list(range(7)))
transpile(qc, basis_gates=['u', 'cx']).count_ops()["cx"]
```

# Ancilla-Free Quantum Adder with Sublinear Depth

## Overview
Implementation of an **ancilla-free quantum adder** based on [Ancilla-free Quantum Adder with Sublinear Depth](http://arxiv.org/abs/2501.16802)
that achieves **sublinear depth** while using only classical reversible logic (Toffoli, CNOT, and X gates). The
implementation is based on recent techniques in optimising CNOT and Toffoli ladder circuits, allowing for efficient
ripple-carry addition without requiring ancilla qubits. In particular, this implementation uses the
approach of **conditionally clean ancillae** [2] to optimise the depth of the adder circuit.

- Optimised CNOT Ladder: Reduces an n-CNOT ladder to O(log n) depth while maintaining O(n) gate count.
- Optimised Toffoli Ladder: Substitutes an n-Toffoli ladder with a circuit of O(log² n) depth using O(n log n) gates.
- Efficient Ripple-Carry Addition: Performs n-bit addition with an overall O(log² n) depth.

## Noted Typos in the Paper
The following typos were present in the paper as of Jan 30, 2025 when I first began implementing this
work. These were verified with the author via email correspondence:

- **Algorithm 2:** Line 6 should read $X' \leftarrow X_{\alpha_0}$.
- **Algorithm 2:** Line 9 should read $C_L \leftarrow \text{MCX}(X_{\alpha_{k-3}}, \cdots, X_{\alpha_{k-2}} )$
- **Algorithm 2:** Line 17 should read $\bm{\alpha}' \leftarrow \bm{\alpha}'::\alpha_{k-3} - \alpha_0 - k/2 + 2$

## Installation
```bash
git clone https://github.com/patelvyom/QAdderSubLinearDepth.git
cd QAdderSubLinearDepth
pip install -r requirements.txt
```

## Usage
```python
from adder import AdderSublinearDepth
n = 5
adder_gate = AdderSublinearDepth(n)
adder_gate.definition.draw()
```

## References
1. [Ancilla-free Quantum Adder with Sublinear Depth](http://arxiv.org/abs/2501.16802)
2. [Rise of conditionally clean ancillae for optimizing quantum circuits](https://arxiv.org/abs/2407.17966)
