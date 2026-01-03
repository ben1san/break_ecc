from qiskit import QuantumRegister, QuantumCircuit

class QuantumECC:
    def __init__(self, arithmetic):
        self.arith = arithmetic
        self.p = arithmetic.N

    def calculate_H_R(self, circuit, P_regs, const_point, ancilla_regs):
        X1, Y1, Z1 = P_regs
        x2, y2 = const_point
        
        # ビット幅 n を取得
        n = len(X1) 

        # 【修正】ancilla_regs を n ビットごとの塊（レジスタ）として定義
        T1 = ancilla_regs[0*n : 1*n]
        T2 = ancilla_regs[1*n : 2*n]
        T3 = ancilla_regs[2*n : 3*n]
        T4 = ancilla_regs[3*n : 4*n]
        
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
        
        # ビット幅 n を取得
        n = len(X1)
        
        # 【修正】スライスを使ってレジスタ定義
        T1 = ancilla_regs[0*n : 1*n] # Z3 will be stored here
        T2 = ancilla_regs[1*n : 2*n] # H is here
        T3 = ancilla_regs[2*n : 3*n]
        T4 = ancilla_regs[3*n : 4*n]

        # Uncompute T3(Z^3), T1(Z^2) using inverse ops
        self.arith.modular_general_multiply(circuit, Z1, T1, T3, inverse=True)
        self.arith.modular_square(circuit, Z1, T1, inverse=True)
        # Calculate Z3 into T1 (Z3 = Z1 * H)
        self.arith.modular_general_multiply(circuit, Z1, T2, T1)

    def calculate_X3_Y3_and_final_cleanup(self, circuit, P_regs, const_point, ancilla_regs):
        X1, Y1, Z1 = P_regs
        x2, y2 = const_point
        
        # ビット幅 n を取得
        n = len(X1)
        
        # 【修正】ここも全てスライスに変更
        Z3_reg = ancilla_regs[0*n : 1*n] # T1
        H_reg  = ancilla_regs[1*n : 2*n] # T2
        H2_reg = ancilla_regs[2*n : 3*n] # T3
        R_reg  = ancilla_regs[3*n : 4*n] # T4
        V_reg  = ancilla_regs[4*n : 5*n] # T5
        X3_reg = ancilla_regs[5*n : 6*n] # T6
        Y3_reg = ancilla_regs[6*n : 7*n] # T7
        tmp_reg = ancilla_regs[7*n : 8*n] # T8

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


class ScalarMultiplication:
    def __init__(self, quantum_ecc, arithmetic, a=0):
        self.ecc = quantum_ecc
        self.arith = arithmetic
        self.a = a

    def create_controlled_add_gate(self, const_point, reg_specs):
        """指定された定数点 const_point を加算する制御付きゲートを作成"""
        n = reg_specs['p_size']
        q_p = [QuantumRegister(n, name=reg) for reg in ['X1', 'Y1', 'Z1']]
        q_anc = QuantumRegister(reg_specs['anc_size'], name='ancilla')
        
        qc_add = QuantumCircuit(*q_p, q_anc)
        P_regs = [q_p[0], q_p[1], q_p[2]]
        ancilla_regs = list(q_anc)

        # ECC加算ロジック
        self.ecc.calculate_H_R(qc_add, P_regs, const_point, ancilla_regs)
        self.ecc.calculate_Z3_and_cleanup(qc_add, P_regs, ancilla_regs)
        self.ecc.calculate_X3_Y3_and_final_cleanup(qc_add, P_regs, const_point, ancilla_regs)

        # 結果をメインレジスタへSWAP
        Z3_reg = ancilla_regs[0*n : 1*n]
        X3_reg = ancilla_regs[5*n : 6*n]
        Y3_reg = ancilla_regs[6*n : 7*n]

        for i in range(n):
            qc_add.swap(P_regs[0][i], X3_reg[i])
            qc_add.swap(P_regs[1][i], Y3_reg[i])
            qc_add.swap(P_regs[2][i], Z3_reg[i])

        gate_label = f"Add({const_point[0]},{const_point[1]})"
        return qc_add.to_gate(label=gate_label).control(1)

    def _classical_point_doubling(self, point, p):
        """古典的な点2倍算 (次の加算点を計算するため)"""
        x, y = point
        if x is None or y == 0: return (None, None)
        numerator = (3 * x**2 + self.a) % p
        denominator = (2 * y) % p
        inv = pow(denominator, -1, p)
        lam = (numerator * inv) % p
        x3 = (lam**2 - 2*x) % p
        y3 = (lam * (x - x3) - y) % p
        return (x3, y3)

    def build_scalar_mult_circuit(self, circuit, k_regs, P_regs, ancilla_regs, base_point_P, p_mod):
        """k_regsの各ビットに応じて、2^i * P を累積加算するループを構築"""
        reg_specs = {
            'p_size': len(P_regs[0]),
            'anc_size': len(ancilla_regs)
        }
        current_Q = base_point_P
        
        # kの最下位ビットから順に処理
        for i in range(len(k_regs)):
            # 1. 制御付き加算ゲート作成
            ctrl_gate = self.create_controlled_add_gate(current_Q, reg_specs)
            
            # 2. 回路に適用
            qubits = [k_regs[i]] + list(P_regs[0]) + list(P_regs[1]) + list(P_regs[2]) + list(ancilla_regs)
            circuit.append(ctrl_gate, qubits)
            
            # 3. 定数点を2倍に更新
            current_Q = self._classical_point_doubling(current_Q, p_mod)
            
    # 検証用ヘルパー
    def _classical_scalar_mult(self, k, point, p):
        current_P = point
        result = (None, None)
        for i in range(k.bit_length()):
            if (k >> i) & 1:
                result = self._classical_add(result, current_P, p)
            current_P = self._classical_point_doubling(current_P, p)
        return result

    def _classical_add(self, p1, p2, N):
        if p1 == (None, None): return p2
        if p2 == (None, None): return p1
        x1, y1 = p1
        x2, y2 = p2
        if x1 == x2 and y1 == y2:
            return self._classical_point_doubling(p1, N)
        num = (y2 - y1) % N
        den = (x2 - x1) % N
        inv = pow(den, -1, N)
        lam = (num * inv) % N
        x3 = (lam**2 - x1 - x2) % N
        y3 = (lam * (x1 - x3) - y1) % N
        return (x3, y3)