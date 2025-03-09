from qiskit.circuit import Gate
from qiskit import QuantumCircuit, QuantumRegister
from typing import List, Tuple


def get_linear_depth_ladder_ops(qreg: List[int]) -> Tuple[QuantumCircuit, List[int]]:
    """
    Args:
        qreg: List of qubit indices to apply the ladder operations on. qreg[0] is assumed to be ancilla.
    Returns:
        QuantumCircuit: Linear-depth ladder operations.
        int: Index of control qubit to apply the final CCX gate.
    """
    n = len(qreg)
    assert n > 3, "n = n_ctrls + 1 => n_ctrls >= 3 to use MCX ladder. Otherwise, use CCX"
    qc = QuantumCircuit(n)

    # up-ladder
    for i in range(2, n - 2, 2):
        qc.ccx(qreg[i + 1], qreg[i + 2], qreg[i])
        qc.x(qreg[i])

    # down-ladder
    if n % 2 != 0:
        x, y, t = n - 3, n - 5, n - 6
    else:
        x, y, t = n - 1, n - 4, n - 5

    if t > 0:
        qc.ccx(qreg[x], qreg[y], qreg[t])
        qc.x(qreg[t])

    for i in range(t, 2, -2):
        qc.ccx(qreg[i], qreg[i - 1], qreg[i - 2])
        qc.x(qreg[i - 2])

    mid_second_ctrl = 1 + max(0, 6 - n)
    final_ctrl = qreg[mid_second_ctrl] - 1
    return qc, final_ctrl


class MCXLinearDepth(Gate):
    def __init__(self, num_controls, clean=True):
        self.n_ctrls = num_controls
        self.n_anc = 1
        self.n_qubits = num_controls + self.n_anc + 1  # control + ancilla + target
        self.clean = clean
        super().__init__(f"linear_mcx_{num_controls}_{self.n_anc}", self.n_qubits, [])

    def _define(self):
        ctrl, targ = QuantumRegister(self.n_ctrls, "ctrl"), QuantumRegister(1, "targ")
        qc = QuantumCircuit(ctrl, targ)

        if self.n_ctrls <= 2:
            qc.mcx(ctrl, targ)
        else:
            anc = QuantumRegister(self.n_anc, "anc")
            qc.add_register(anc)
            ladder_ops, final_ctrl = get_linear_depth_ladder_ops(list(range(self.n_ctrls + self.n_anc)))
            qc.ccx(ctrl[0], ctrl[1], anc)  #                                      # create conditionally clean ancilla
            qc.compose(ladder_ops, anc[:] + ctrl[:], inplace=True)  #             # up-ladder
            qc.ccx(anc, ctrl[final_ctrl], targ)  #                                # target
            qc.compose(ladder_ops.inverse(), anc[:] + ctrl[:], inplace=True)  #   # down-ladder
            qc.ccx(ctrl[0], ctrl[1], anc)

            if not self.clean:
                # toggle-detection if dirty ancilla
                qc.compose(ladder_ops, anc[:] + ctrl[:], inplace=True)
                qc.ccx(anc, ctrl[final_ctrl], targ)
                qc.compose(ladder_ops.inverse(), anc[:] + ctrl[:], inplace=True)

        self.definition = qc


class CCXN(Gate):
    def __init__(self, n):
        self.n = n
        self.n_qubits = 3 * n
        super().__init__(f"ccxn_{n}", self.n_qubits, [])

    def _define(self):
        x, y, t = (
            QuantumRegister(self.n, "x"),
            QuantumRegister(self.n, "y"),
            QuantumRegister(self.n, "t"),
        )
        qc = QuantumCircuit(x, y, t)
        for x, y, t in zip(x, y, t):
            qc.x(t)
            qc.ccx(x, y, t)

        self.definition = qc


def build_logn_depth_ccx_ladder(
    alloc_anc: int, ctrls: List[int], skip_cond_clean=False
) -> Tuple[QuantumCircuit, List[int]]:
    qc = QuantumCircuit(len(ctrls) + 1)
    anc = [alloc_anc]
    final_ctrls = []

    while len(ctrls) > 1:
        next_batch_len = min(len(anc) + 1, len(ctrls))
        ctrls, nxt_batch = ctrls[next_batch_len:], ctrls[:next_batch_len]
        new_anc = []
        while len(nxt_batch) > 1:
            ccx_n = len(nxt_batch) // 2
            st = int(len(nxt_batch) % 2)
            ccx_x, ccx_y, ccx_t = (
                nxt_batch[st : st + ccx_n],
                nxt_batch[st + ccx_n :],
                anc[-ccx_n:],
            )
            assert len(ccx_x) == len(ccx_y) == len(ccx_t) == ccx_n >= 1
            if ccx_t != [alloc_anc]:
                qc.compose(CCXN(ccx_n).definition, ccx_x + ccx_y + ccx_t, inplace=True)
            else:
                if not skip_cond_clean:
                    qc.ccx(ccx_x[0], ccx_y[0], ccx_t[0])  #   # create conditionally clean ancilla
            new_anc += nxt_batch[st:]  #                      # newly created conditionally clean ancilla
            nxt_batch = ccx_t + nxt_batch[:st]
            anc = anc[:-ccx_n]

        anc = sorted(anc + new_anc)
        final_ctrls += nxt_batch

    final_ctrls += ctrls
    final_ctrls = sorted(final_ctrls)
    return qc, final_ctrls[:-1]  #                            # exclude ancilla


class MCXLogDepth(Gate):
    def __init__(self, num_controls, clean=True):
        self.n_ctrl = num_controls
        self.n_anc = 2
        self.num_qubits = num_controls + self.n_anc + 1  #                   # control + ancilla + target
        self.clean = clean
        super().__init__(f"log_mcx_{num_controls}_{self.n_anc}", self.num_qubits, [])

    def _define(self):
        ctrl, targ = QuantumRegister(self.n_ctrl, "ctrl"), QuantumRegister(1, "targ")
        qc = QuantumCircuit(ctrl, targ)

        if self.n_ctrl <= 2:
            qc.mcx(ctrl, targ)
        else:
            anc = QuantumRegister(self.n_anc, "anc")
            qc.add_register(anc)
            ladder_ops, final_ctrls = build_logn_depth_ccx_ladder(self.n_ctrl, list(range(self.n_ctrl)))
            qc.compose(ladder_ops, ctrl[:] + [anc[0]], inplace=True)
            if len(final_ctrls) == 1:  #                                     # Already a toffoli
                qc.ccx(anc[0], ctrl[final_ctrls[0]], targ)
            else:
                mid_mcx = MCXLinearDepth(len(final_ctrls) + 1, clean=True)
                qc.compose(
                    mid_mcx.definition,
                    [anc[0]] + ctrl[final_ctrls] + targ[:] + [anc[1]], #     # ctrls, targ, anc
                    inplace=True,
                )
            qc.compose(ladder_ops.inverse(), ctrl[:] + [anc[0]], inplace=True)

            if not self.clean:
                # toggle-detection if not clean
                ladder_ops_new, final_ctrls = build_logn_depth_ccx_ladder(
                    self.n_ctrl, list(range(self.n_ctrl)), skip_cond_clean=True
                )
                qc.compose(ladder_ops_new, ctrl[:] + [anc[0]], inplace=True)
                if len(final_ctrls) == 1:
                    qc.ccx(anc[0], ctrl[final_ctrls[0]], targ)
                else:
                    qc.compose(
                        mid_mcx.definition,
                        [anc[0]] + ctrl[final_ctrls] + targ[:] + [anc[1]],
                        inplace=True,
                    )
                qc.compose(ladder_ops_new.inverse(), ctrl[:] + [anc[0]], inplace=True)

        self.definition = qc


def mcx_ladder(N, alpha):
    """
    Create a log-depth MCX ladder circuit.

    Args:
        N (int): Number of MCX gates in the ladder.
        alpha (int): Number of controls per MCX gate.
    """

    def helper(X, alphas):
        k = len(alphas) + 1
        if k == 1:
            return []
        if k == 2:
            return [X[: alphas[0] + 1]]

        X_bar = [X[alphas[0]]]
        alpha_bar = []
        right_pairs = [X[0 : alphas[0] + 1]]
        left_pairs = [X[alphas[k - 3] : alphas[-1] + 1]]

        for i in range(1, ceil(k / 2) - 1):
            left_pairs += [X[alphas[2 * i - 2] : alphas[2 * i - 1] + 1]]
            right_pairs += [X[alphas[2 * i - 1] : alphas[2 * i] + 1]]
            X_bar += X[alphas[2 * i - 2] + 1 : alphas[2 * i - 1]] + X[alphas[2 * i - 1] + 1 : alphas[2 * i] + 1]
            alpha_bar.append(alphas[2 * i] - alphas[0] - i)

        if k % 2 == 0:
            X_bar += X[alphas[k - 4] + 1 : alphas[k - 3] + 1]
            alpha_bar.append(alphas[k - 3] - alphas[0] - int(k / 2) + 2)

        return left_pairs + helper(X_bar, alpha_bar) + right_pairs

    n = N * alpha + 1
    qc = QuantumCircuit(n)
    X, alphas = list(range(n)), list(range(alpha, n, alpha))
    mcxs = helper(X, alphas)
    for mcx in mcxs:
        if len(mcx) <= 3:  # already a Toffoli
            qc.mcx(mcx[:-1], mcx[-1])
        else:
            # for each mcx with n_ctrls > 2, use 2 qubits above the first ctrl index as ancillae
            ancilla_idx = [mcx[0] - 2, mcx[0] - 1]
            gate = MCXLogDepth(len(mcx) - 1).definition
            qc.compose(
                gate,
                # ctrls, targ, anc
                mcx[:-1] + [mcx[-1]] + ancilla_idx,
                inplace=True,
            )

    return qc


class AdderSublinearDepth(Gate):
    def __init__(self, n_qubits: int, name: str = "AncillaFreeSublinearDepthAdder") -> None:
        if n_qubits < 1:
            raise ValueError("The number of qubits must be at least 1.")
        super().__init__(name=name, num_qubits=2 * n_qubits + 1, params=[])
        self.n_bits = n_qubits

    def _define(self):
        n = self.n_bits
        qr_a = QuantumRegister(n, "a")
        qr_b = QuantumRegister(n, "b")
        z = QuantumRegister(1, "z")
        qc = QuantumCircuit(qr_a, qr_b, z)

        ab_interleaved = [q for pair in zip(qr_a, qr_b) for q in pair]

        qc.cx(qr_a[1:], qr_b[1:])
        qc.compose(mcx_ladder(n - 1, 1), qr_a[1:] + z[:], inplace=True)
        qc.compose(mcx_ladder(n, 2).inverse(), ab_interleaved + z[:], inplace=True)
        qc.cx(qr_a[1:], qr_b[1:])
        qc.x(qr_b[1 : n - 1])
        qc.compose(mcx_ladder(n - 1, 2), ab_interleaved[:-1], inplace=True)
        qc.x(qr_b[1 : n - 1])
        qc.compose(mcx_ladder(n - 2, 1).inverse(), qr_a[1:], inplace=True)
        qc.cx(qr_a, qr_b)

        self.definition = qc