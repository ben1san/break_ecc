import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit.library import QFT

class ModularArithmetic:
    def __init__(self, N, n_qubits):
        self.N = N
        self.n = n_qubits

    def _get_qft(self, n, inverse=False):
        return QFT(n, do_swaps=True, inverse=inverse).to_gate()

    def _get_phase_add(self, n, val):
        qc = QuantumCircuit(n, name=f"PhiAdd({val:.2f})")
        for i in range(n):
            angle = 2 * np.pi * val / (2**(n - i))
            if abs(angle) > 1e-9:
                qc.p(angle, i)
        return qc.to_gate()

    def cc_phase_add(self, circuit, ctrl_list, target_reg, val):
        """
        修正版: 内部で % N を行わず、渡された val をそのまま位相回転させる。
        これにより、inverse=True 時に負の val を渡すことで正確な引き算が可能。
        """
        gate = self._get_phase_add(len(target_reg), val)
        if len(ctrl_list) > 0:
            gate = gate.control(len(ctrl_list))
            circuit.append(gate, list(ctrl_list) + list(target_reg))
        else:
            circuit.append(gate, list(target_reg))

    def modular_square(self, circuit, src_reg, out_reg, inverse=False):
        n_in = len(src_reg)
        n_out = len(out_reg)
        
        if not inverse:
            circuit.append(self._get_qft(n_out), out_reg)
            sign = 1
        else:
            circuit.append(self._get_qft(n_out), out_reg)
            sign = -1

        for i in range(n_in):
            val = (2**(2*i)) % self.N
            self.cc_phase_add(circuit, [src_reg[i]], out_reg, sign * val)
            for j in range(i + 1, n_in):
                val = (2 * (2**(i+j))) % self.N
                self.cc_phase_add(circuit, [src_reg[i], src_reg[j]], out_reg, sign * val)

        circuit.append(self._get_qft(n_out, inverse=True), out_reg)

    def modular_general_multiply(self, circuit, reg_a, reg_b, out_reg, inverse=False):
        n_a = len(reg_a)
        n_b = len(reg_b)
        n_out = len(out_reg)
        
        circuit.append(self._get_qft(n_out), out_reg)
        sign = 1 if not inverse else -1
        
        for i in range(n_a):
            for j in range(n_b):
                val = (2**(i+j)) % self.N
                self.cc_phase_add(circuit, [reg_a[i], reg_b[j]], out_reg, sign * val)
                
        circuit.append(self._get_qft(n_out, inverse=True), out_reg)

    def modular_scalar_mult(self, circuit, src_reg, out_reg, scalar, inverse=False):
        n_in = len(src_reg)
        n_out = len(out_reg)
        
        circuit.append(self._get_qft(n_out), out_reg)
        sign = 1 if not inverse else -1
        
        for i in range(n_in):
            val = (scalar * (2**i)) % self.N
            self.cc_phase_add(circuit, [src_reg[i]], out_reg, sign * val)
            
        circuit.append(self._get_qft(n_out, inverse=True), out_reg)

    def modular_sub(self, circuit, src_reg, target_reg, inverse=False):
        n_src = len(src_reg)
        n_target = len(target_reg)
        
        if not inverse:
            circuit.append(self._get_qft(n_target), target_reg)
            sign = 1
        else:
            circuit.append(self._get_qft(n_target), target_reg)
            sign = -1

        for i in range(n_src):
            val_base = (self.N - (2**i)) % self.N
            self.cc_phase_add(circuit, [src_reg[i]], target_reg, sign * val_base)
            
        circuit.append(self._get_qft(n_target, inverse=True), target_reg)