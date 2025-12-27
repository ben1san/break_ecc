import time
import operator
import warnings
# Qiskitのバージョンによる警告を抑制
warnings.filterwarnings('ignore', category=DeprecationWarning)

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import QFTGate  # 修正: QFTクラスではなくQFTGateを使用
from qiskit_aer import AerSimulator

from arithmetic import ModularArithmetic
from ecc import QuantumECC, ScalarMultiplication

# --- 1. ショアのアルゴリズム統合クラス ---
class ShorECDLP:
    def __init__(self, scalar_mult_system, p_mod):
        self.sm_system = scalar_mult_system
        self.p_mod = p_mod

    def construct_circuit(self, num_ctrl_qubits, point_P, point_Q):
        n = self.sm_system.arith.n
        reg_a = QuantumRegister(num_ctrl_qubits, name='reg_a')
        reg_b = QuantumRegister(num_ctrl_qubits, name='reg_b')
        reg_X = QuantumRegister(n, name='X')
        reg_Y = QuantumRegister(n, name='Y')
        reg_Z = QuantumRegister(n, name='Z')
        reg_ancilla = QuantumRegister(8 * n, name='ancilla')
        cr_a = ClassicalRegister(num_ctrl_qubits, name='c_a')
        cr_b = ClassicalRegister(num_ctrl_qubits, name='c_b')

        qc = QuantumCircuit(reg_a, reg_b, reg_X, reg_Y, reg_Z, reg_ancilla, cr_a, cr_b)

        # 初期化
        qc.h(reg_a)
        qc.h(reg_b)
        qc.x(reg_Z[0]) 

        P_regs = [reg_X, reg_Y, reg_Z]
        ancilla_regs = list(reg_ancilla)

        # ダブルスカラー倍算
        self.sm_system.build_scalar_mult_circuit(qc, reg_a, P_regs, ancilla_regs, point_P, self.p_mod)
        self.sm_system.build_scalar_mult_circuit(qc, reg_b, P_regs, ancilla_regs, point_Q, self.p_mod)

        # 逆QFT (QFTGateを使用)
        iqft = QFTGate(num_ctrl_qubits).inverse()
        qc.append(iqft, reg_a)
        qc.append(iqft, reg_b)

        qc.measure(reg_a, cr_a)
        qc.measure(reg_b, cr_b)
        return qc

# --- 2. 古典後処理クラス ---
class ShorPostProcessor:
    def __init__(self, curve_order, p_mod, a, b):
        self.r = curve_order
        self.p = p_mod
        self.a = a
        self.b = b

    def _classical_point_mult(self, k, point):
        if k == 0: return (None, None)
        if k == 1: return point
        result = (None, None)
        add_point = point
        for i in range(k.bit_length()):
            if (k >> i) & 1:
                result = self._point_add(result, add_point)
            add_point = self._point_double(add_point)
        return result

    def _point_add(self, p1, p2):
        if p1 == (None, None): return p2
        if p2 == (None, None): return p1
        x1, y1 = p1
        x2, y2 = p2
        if x1 == x2 and y1 != y2: return (None, None)
        if x1 == x2 and y1 == y2: return self._point_double(p1)
        num = (y2 - y1) % self.p
        den = (x2 - x1) % self.p
        try:
            inv = pow(den, -1, self.p)
        except ValueError: return (None, None) # ゼロ除算などは無視
        lam = (num * inv) % self.p
        x3 = (lam**2 - x1 - x2) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return (x3, y3)

    def _point_double(self, p1):
        if p1 == (None, None): return (None, None)
        x1, y1 = p1
        if y1 == 0: return (None, None)
        num = (3 * x1**2 + self.a) % self.p
        den = (2 * y1) % self.p
        try:
            inv = pow(den, -1, self.p)
        except ValueError: return (None, None)
        lam = (num * inv) % self.p
        x3 = (lam**2 - 2*x1) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return (x3, y3)

    def solve_k(self, counts, point_P, point_Q):
        sorted_counts = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
        print(f"  Analyzing {len(counts)} patterns...")
        
        # 候補ごとの得票数を集計
        k_votes = {}

        for bitstring, count in sorted_counts:
            clean_bits = bitstring.replace(" ", "")
            mid = len(clean_bits) // 2
            s_b, s_a = clean_bits[:mid], clean_bits[mid:]
            
            val_a = int(s_a, 2)
            val_b = int(s_b, 2)
            
            # 関係式 k = b/a mod r を仮定 (aP + bQ = 0 => bQ = -aP => Q = (-a/b)P ? いや単純に比率を見る)
            # 実際には Q = kP で位相推定を行うと、測定値 u, v に対して v/u ≒ k などの関係が出る
            # ここではシンプルに k = b * a^-1 mod r を試行
            
            if val_a % self.r == 0: continue
            
            try:
                inv_a = pow(val_a, -1, self.r)
                k_cand = (val_b * inv_a) % self.r
                k_votes[k_cand] = k_votes.get(k_cand, 0) + count
            except ValueError: continue

        # 得票数の多い順に検証
        sorted_candidates = sorted(k_votes.items(), key=operator.itemgetter(1), reverse=True)
        
        for cand_k, votes in sorted_candidates:
            print(f"  Checking candidate k={cand_k} (Votes: {votes})...", end="")
            calc_Q = self._classical_point_mult(cand_k, point_P)
            if calc_Q == point_Q:
                print(" ✅ VALID")
                return cand_k
            else:
                print(" ❌ Invalid")
                
        return None

# --- 3. メインルーチン ---
def verify_qday():
    test_cases = [
        {
            "name": "Bit size 4",
            "p": 13, "order": 7, "n_qubits": 4, 
            "G": (11, 5), "Q": (11, 8), "true_d": 6,
            "ctrl_bits": 5  # 修正: 3 -> 5 に増加
        }
    ]

    for case in test_cases:
        print(f"\n=== Testing {case['name']} (p={case['p']}) ===")
        print(f"Target Private Key d: {case['true_d']}")
        
        arith = ModularArithmetic(case['p'], case['n_qubits'])
        ecc = QuantumECC(arith)
        sm = ScalarMultiplication(ecc, arith, a=0) 
        shor = ShorECDLP(sm, case['p'])
        
        print(f"Constructing circuit (Control Bits: {case['ctrl_bits']})...")
        qc = shor.construct_circuit(case['ctrl_bits'], case['G'], case['Q'])
        
        simulator = AerSimulator(method='matrix_product_state')
        print("Transpiling...")
        t_qc = transpile(qc, simulator)
        
        print(f"Running simulation (Depth: {t_qc.depth()})...")
        # 修正: ショット数を128に増加
        result = simulator.run(t_qc, shots=128).result()
        counts = result.get_counts()
        print("Raw Counts (Top 5):", sorted(counts.items(), key=lambda x: x[1], reverse=True)[:5])
        
        processor = ShorPostProcessor(case['order'], case['p'], a=0, b=7)
        found_k = processor.solve_k(counts, case['G'], case['Q'])
        
        if found_k == case['true_d']:
            print(f"SUCCESS! 🎯 Found private key: {found_k}")
        else:
            print(f"FAILED. Found: {found_k}, Expected: {case['true_d']}")

if __name__ == "__main__":
    verify_qday()