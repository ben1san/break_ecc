import operator
import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit.circuit.library import QFTGate
from qiskit_aer import AerSimulator

class QDayOptimizedSolver:
    def __init__(self, p, a, b):
        self.p = p
        self.a = a
        self.b = b

    def _point_add(self, p1, p2):
        if p1 == (None, None): return p2
        if p2 == (None, None): return p1
        x1, y1 = p1
        x2, y2 = p2
        if x1 == x2 and y1 != y2: return (None, None)
        if x1 == x2 and y1 == y2: return self._point_double(p1)
        
        num = (y2 - y1) % self.p
        den = (x2 - x1) % self.p
        if den == 0: return (None, None)
        
        lam = (num * pow(den, -1, self.p)) % self.p
        x3 = (lam**2 - x1 - x2) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return (x3, y3)

    def _point_double(self, p1):
        if p1 == (None, None): return (None, None)
        x1, y1 = p1
        if y1 == 0: return (None, None)
        
        num = (3 * x1**2 + self.a) % self.p
        den = (2 * y1) % self.p
        if den == 0: return (None, None)
        
        lam = (num * pow(den, -1, self.p)) % self.p
        x3 = (lam**2 - 2*x1) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return (x3, y3)

    def _scalar_mult(self, k, point):
        res = (None, None)
        temp = point
        for i in range(k.bit_length()):
            if (k >> i) & 1:
                res = self._point_add(res, temp)
            temp = self._point_double(temp)
        return res

    def create_oracle_circuit(self, point_P, point_Q, n_ctrl):
        """
        |a>|b>|0> -> |a>|b>|aP + bQ> を実現する最適化回路を作成
        """
        # ターゲットレジスタのビット数（x, y座標それぞれ）
        n_target = self.p.bit_length()
        
        # レジスタ定義
        qa = QuantumRegister(n_ctrl, 'a')
        qb = QuantumRegister(n_ctrl, 'b')
        # ターゲットは (x, y) 座標 + 無限遠点フラグ(inf)
        # 簡略化のため、(0,0)を無限遠点として扱う場合はフラグ不要だが、
        # 正確を期して (x, y) レジスタを用意
        qx = QuantumRegister(n_target, 'x')
        qy = QuantumRegister(n_target, 'y')
        
        qc = QuantumCircuit(qa, qb, qx, qy)
        
        # --- Lookup Tableの構築 ---
        # すべての a, b の組み合わせについて計算し、制御ゲートでターゲットをセットする
        print(f"Compiling Oracle for {2**(2*n_ctrl)} states...")
        
        for a_val in range(2**n_ctrl):
            for b_val in range(2**n_ctrl):
                # 1. 古典的に計算
                P_part = self._scalar_mult(a_val, point_P)
                Q_part = self._scalar_mult(b_val, point_Q)
                Res = self._point_add(P_part, Q_part)
                
                if Res == (None, None):
                    # 無限遠点の場合は何もしない（|00...0>のまま）
                    # または特定のコード（例: 全ビット1）を割り当てる
                    continue
                
                rx, ry = Res
                
                # 2. 制御状態の認識 (Multi-Controlled X)
                # 制御ビットが現在の (a_val, b_val) に一致するときだけ発火させる
                # 一致させるために、0の部分はXで反転させてから制御し、終わったら戻す
                
                # a_valのビットパターンに合わせてXゲート
                for i in range(n_ctrl):
                    if not ((a_val >> i) & 1):
                        qc.x(qa[i])
                
                # b_valのビットパターンに合わせてXゲート
                for i in range(n_ctrl):
                    if not ((b_val >> i) & 1):
                        qc.x(qb[i])
                
                # 3. ターゲットへの書き込み (Toffoli / MCX)
                # 制御ビットすべてが1のとき、ターゲットを rx, ry にする
                ctrl_qubits = list(qa) + list(qb)
                
                # X座標の書き込み
                for i in range(n_target):
                    if (rx >> i) & 1:
                        qc.mcx(ctrl_qubits, qx[i])
                
                # Y座標の書き込み
                for i in range(n_target):
                    if (ry >> i) & 1:
                        qc.mcx(ctrl_qubits, qy[i])
                
                # 4. 制御状態の復元 (Xゲートを戻す)
                for i in range(n_ctrl):
                    if not ((b_val >> i) & 1):
                        qc.x(qb[i])
                for i in range(n_ctrl):
                    if not ((a_val >> i) & 1):
                        qc.x(qa[i])
                        
        return qc

    def build_shor_circuit(self, point_P, point_Q, n_ctrl):
        # メイン回路
        n_target = self.p.bit_length()
        qa = QuantumRegister(n_ctrl, 'a')
        qb = QuantumRegister(n_ctrl, 'b')
        qx = QuantumRegister(n_target, 'x')
        qy = QuantumRegister(n_target, 'y')
        ca = ClassicalRegister(n_ctrl, 'read_a')
        cb = ClassicalRegister(n_ctrl, 'read_b')
        
        qc = QuantumCircuit(qa, qb, qx, qy, ca, cb)
        
        # 1. 重ね合わせ
        qc.h(qa)
        qc.h(qb)
        
        # 2. オラクル (aP + bQ)
        oracle = self.create_oracle_circuit(point_P, point_Q, n_ctrl)
        # 回路を結合
        qc.compose(oracle, qubits=list(qa)+list(qb)+list(qx)+list(qy), inplace=True)
        
        # 3. 逆QFT
        iqft = QFTGate(n_ctrl).inverse()
        qc.append(iqft, qa)
        qc.append(iqft, qb)
        
        # 4. 測定
        qc.measure(qa, ca)
        qc.measure(qb, cb)
        
        return qc

# --- 実行パート ---
def run_qday_optimized():
    # Bit size 4 (p=13) の設定
    p = 13
    a_curve = 0
    b_curve = 7 # y^2 = x^3 + 7
    P = (11, 5)
    Q = (11, 8) # Q = 6P
    
    # 制御ビット数 (精度のため少し多めに)
    # ただし、オラクルのサイズは 2^(2*n_ctrl) に比例して大きくなるため、
    # 3〜4程度が現実的です。
    n_ctrl = 6
    
    solver = QDayOptimizedSolver(p, a_curve, b_curve)
    
    print(f"Building Optimized Circuit for p={p}, n_ctrl={n_ctrl}...")
    qc = solver.build_shor_circuit(P, Q, n_ctrl)
    
    # シミュレータで実行
    sim = AerSimulator()
    t_qc = transpile(qc, sim)
    print(f"Circuit Depth: {t_qc.depth()}")
    
    result = sim.run(t_qc, shots=1024).result()
    counts = result.get_counts()
    
    # 結果表示 (Top 10)
    print("\nSimulation Results:")
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:10]
    for k, v in sorted_counts:
        print(f"State {k}: {v} shots")
        
    # 解説
    print("\n--- Verification ---")
    print(f"Correct d = 6")
    # 最頻値の解析
    top_bitstring = sorted_counts[0][0]
    # bitstringは "read_b read_a" の順 (Qiskitの仕様による)
    parts = top_bitstring.split()
    if len(parts) == 2:
        sb, sa = parts
    else:
        # splitできない場合(空白なし)の対策
        sb = top_bitstring[:n_ctrl]
        sa = top_bitstring[n_ctrl:]
        
    val_b = int(sb, 2)
    val_a = int(sa, 2)
    print(f"Top result: b={val_b}, a={val_a}")
    
    # ショアのアルゴリズムの関係式: b + d*a = 0 (mod r) または b/a = -d ?
    # ECDLPの場合の位相推定結果の解釈:
    # Q = dP とすると、状態は |a>|b>|aP+bQ> = |a>|b>|(a+bd)P>
    # 周期性から、測定値 a, b は a + b*d = 0 (mod r) を満たす確率が高い
    # つまり d = -a * b^-1 (mod r)
    
    r = 7 # 位数
    if val_b % r != 0:
        inv_b = pow(val_b, -1, r)
        calc_d = (-val_a * inv_b) % r
        print(f"Calculated d = -a/b mod r = -{val_a}/{val_b} = {calc_d}")
    else:
        print("b is multiple of r, cannot divide.")

if __name__ == "__main__":
    run_qday_optimized()