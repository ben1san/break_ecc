from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.library import QFT

class ShorECDLP:
    def __init__(self, scalar_mult_system, p_mod):
        """
        Args:
            scalar_mult_system: ecc.ScalarMultiplication のインスタンス
            p_mod: 楕円曲線の定義体 (素数 p)
        """
        self.sm_system = scalar_mult_system
        self.p_mod = p_mod

    def construct_circuit(self, num_ctrl_qubits, point_P, point_Q, initial_point_Z1=1):
        """
        ショアのアルゴリズム全体回路を構築します。
        
        Args:
            num_ctrl_qubits (int): 制御レジスタ(|a>, |b>)それぞれのビット数
            point_P (tuple): ベースポイント P (x, y)
            point_Q (tuple): ベースポイント Q (x, y)
            initial_point_Z1 (int): ターゲットレジスタの初期Z値（通常1）
        
        Returns:
            QuantumCircuit: 構築された回路
        """
        # --- 1. レジスタの定義 ---
        # 制御用レジスタ |a> と |b>
        reg_a = QuantumRegister(num_ctrl_qubits, name='reg_a')
        reg_b = QuantumRegister(num_ctrl_qubits, name='reg_b')
        
        # 楕円曲線点計算用レジスタ (X, Y, Z)
        # ビット幅は arithmetic.n から取得
        n = self.sm_system.arith.n
        reg_X = QuantumRegister(n, name='X')
        reg_Y = QuantumRegister(n, name='Y')
        reg_Z = QuantumRegister(n, name='Z')
        
        # アンシラレジスタ (ecc.py の構成に合わせて確保)
        # calculate_X3_Y3... で最大8ブロック(T1~T8)使用するため 8*n 必要
        num_ancilla = 8 * n
        reg_ancilla = QuantumRegister(num_ancilla, name='ancilla')
        
        # 測定用古典レジスタ
        cr_a = ClassicalRegister(num_ctrl_qubits, name='c_a')
        cr_b = ClassicalRegister(num_ctrl_qubits, name='c_b')

        qc = QuantumCircuit(reg_a, reg_b, reg_X, reg_Y, reg_Z, reg_ancilla, cr_a, cr_b)

        # --- 2. 状態初期化 (QPE 前半) ---
        # 制御レジスタにアダマールゲートを適用して重ね合わせ状態を作成
        qc.h(reg_a)
        qc.h(reg_b)

        # ターゲットレジスタの初期化
        # Projective座標 (X:Y:Z) の初期状態を設定
        # 注意: 通常は無限遠点ですが、算術制約上 (0,0,0) を避ける場合は
        # ここで特定のダミー点や Z=1 をセットします。(要件に応じて調整)
        # ここでは Z=1 のみをセットする例とします (X=0, Y=0, Z=1 => 原点(0,0)相当)
        # 実際の運用では point_P や point_Q とは異なるオフセット点を用います
        if initial_point_Z1 == 1:
            qc.x(reg_Z[0]) # Z=1 (最下位ビット反転)

        P_regs = [reg_X, reg_Y, reg_Z]
        ancilla_regs = list(reg_ancilla)

        # --- 3. ダブルスカラー倍算 (Double Scalar Multiplication) ---
        # 回路ロジック: |Target> = (a*P + b*Q) + Initial_Point
        
        print("Building Scalar Mult for P...")
        # Step 3a: a * P を加算
        # ScalarMultiplicationクラスのメソッドを利用
        self.sm_system.build_scalar_mult_circuit(
            circuit=qc,
            k_regs=reg_a,           # 制御レジスタ a
            P_regs=P_regs,          # ターゲット
            ancilla_regs=ancilla_regs,
            base_point_P=point_P,   # 定数点 P
            p_mod=self.p_mod
        )

        print("Building Scalar Mult for Q...")
        # Step 3b: b * Q を加算
        # 同じターゲットレジスタに対して続けて加算を行うことで aP + bQ を実現
        self.sm_system.build_scalar_mult_circuit(
            circuit=qc,
            k_regs=reg_b,           # 制御レジスタ b
            P_regs=P_regs,          # ターゲット
            ancilla_regs=ancilla_regs,
            base_point_P=point_Q,   # 定数点 Q
            p_mod=self.p_mod
        )

        # --- 4. 逆QFT (Inverse QFT) ---
        # Qiskitのライブラリ関数を使用 (制御レジスタごと)
        print("Appending Inverse QFT...")
        qc.append(QFT(num_ctrl_qubits, inverse=True).to_gate(), reg_a)
        qc.append(QFT(num_ctrl_qubits, inverse=True).to_gate(), reg_b)

        # --- 5. 測定 ---
        qc.measure(reg_a, cr_a)
        qc.measure(reg_b, cr_b)

        return qc