from qiskit import QuantumRegister

class QuantumECC:
    def __init__(self, arithmetic):
        self.arith = arithmetic
        self.p = arithmetic.N

    def calculate_H_R(self, circuit, P_regs, const_point, ancilla_regs):
        X1, Y1, Z1 = P_regs
        x2, y2 = const_point
        # ancilla_regs: [T1(Z3), T2(H), T3(H^2), T4(R), ...]
        T1, T2, T3, T4 = ancilla_regs[:4]
        
        # 1. T1 = Z^2, 2. T2 = U2, 3. T3 = Z^3, 4. T4 = S2
        self.arith.modular_square(circuit, Z1, T1)
        self.arith.modular_scalar_mult(circuit, T1, T2, x2)
        self.arith.modular_general_multiply(circuit, Z1, T1, T3)
        self.arith.modular_scalar_mult(circuit, T3, T4, y2)
        
        # 5. H = U2 - X1, 6. R = S2 - Y1
        self.arith.modular_sub(circuit, X1, T2) # T2 is now H
        self.arith.modular_sub(circuit, Y1, T4) # T4 is now R

    def calculate_Z3_and_cleanup(self, circuit, P_regs, ancilla_regs):
        X1, Y1, Z1 = P_regs
        T1, T2, T3, T4 = ancilla_regs[:4]
        # Uncompute T3(Z^3), T1(Z^2) using inverse ops
        self.arith.modular_general_multiply(circuit, Z1, T1, T3, inverse=True)
        self.arith.modular_square(circuit, Z1, T1, inverse=True)
        # Calculate Z3 into T1
        self.arith.modular_general_multiply(circuit, Z1, T2, T1)

    def calculate_X3_Y3_and_final_cleanup(self, circuit, P_regs, const_point, ancilla_regs):
        X1, Y1, Z1 = P_regs
        x2, y2 = const_point
        
        Z3_reg = ancilla_regs[0] # T1
        H_reg  = ancilla_regs[1] # T2
        H2_reg = ancilla_regs[2] # T3
        R_reg  = ancilla_regs[3] # T4
        V_reg  = ancilla_regs[4] # T5
        X3_reg = ancilla_regs[5] # T6
        Y3_reg = ancilla_regs[6] # T7
        tmp_reg = ancilla_regs[7] # T8

        # --- Step 1: Compute H^2, V ---
        self.arith.modular_square(circuit, H_reg, H2_reg)
        self.arith.modular_general_multiply(circuit, X1, H2_reg, V_reg)

        # --- Step 2: Compute X3 ---
        self.arith.modular_square(circuit, R_reg, X3_reg)
        self.arith.modular_general_multiply(circuit, H_reg, H2_reg, X3_reg, inverse=True) # -H^3
        self.arith.modular_sub(circuit, V_reg, X3_reg) # -V
        self.arith.modular_sub(circuit, V_reg, X3_reg) # -V

        # --- Step 3: Compute Y3 ---
        self.arith.modular_sub(circuit, X3_reg, V_reg) # V = V - X3
        self.arith.modular_general_multiply(circuit, R_reg, V_reg, Y3_reg)
        
        # Y3 -= Y1 * H^3 (via tmp)
        self.arith.modular_general_multiply(circuit, Y1, H_reg, tmp_reg) # tmp = Y1*H
        self.arith.modular_general_multiply(circuit, tmp_reg, H2_reg, Y3_reg, inverse=True) # -tmp*H^2
        self.arith.modular_general_multiply(circuit, Y1, H_reg, tmp_reg, inverse=True) # Uncompute tmp

        # --- Step 4: Global Cleanup ---
        self.arith.modular_sub(circuit, X3_reg, V_reg, inverse=True) # Restore V
        self.arith.modular_general_multiply(circuit, X1, H2_reg, V_reg, inverse=True)
        self.arith.modular_square(circuit, H_reg, H2_reg, inverse=True)
        
        self.arith.modular_sub(circuit, Y1, R_reg, inverse=True) # R -> S2
        self.arith.modular_sub(circuit, X1, H_reg, inverse=True) # H -> U2
        
        # Recompute Z^2, Z^3 to cleanup S2, U2
        self.arith.modular_square(circuit, Z1, tmp_reg)
        self.arith.modular_general_multiply(circuit, Z1, tmp_reg, V_reg)
        
        self.arith.modular_scalar_mult(circuit, V_reg, R_reg, y2, inverse=True)
        self.arith.modular_scalar_mult(circuit, tmp_reg, H_reg, x2, inverse=True)
        
        self.arith.modular_general_multiply(circuit, Z1, tmp_reg, V_reg, inverse=True)
        self.arith.modular_square(circuit, Z1, tmp_reg, inverse=True)