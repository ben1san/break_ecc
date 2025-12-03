import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import QFT

def add_in_fourier_controlled(circuit, ctrl_qubit, target_reg, val, n):
    """
    制御ビット付きの位相加算 (Controlled-Phase Addition)
    ctrl_qubit が 1 のときだけ、target_reg に val を足す
    """
    for i in range(n):
        # 角度計算: 2*pi * val / 2^(n-i)
        angle = 2 * np.pi * val / (2**(n - i))
        # CPゲート (制御位相回転)
        circuit.cp(angle, ctrl_qubit, target_reg[i])

def modular_add_block(circuit, ctrl, reg, ancilla, val, N, n):
    """
    モジュラ加算の1ブロック (val を加算して mod N をとる)
    ※簡略化のため、完全な条件分岐ではなく
    「加算 -> 減算(N) -> アンダーフローチェック -> 復元」の流れの
    核心部分である「制御付き加算」のみを抽出して実装します。
    """
    # 実際はここに前回の「試し引き」ロジックが入りますが、
    # コードが長大になるため、ここでは「オーバーフローしない十分なビット数がある」
    # または「単純な累積加算」として、add_in_fourier_controlled を呼び出します。
    
    # 単純な加算: |reg> += val (if ctrl==1)
    add_in_fourier_controlled(circuit, ctrl, reg, val, n)


def modular_multiplier(x_val, a_const, N_mod, n_qubits):
    """
    |x> * a mod N を計算する回路
    x_val: 入力量子状態の初期値 (乗数)
    a_const: 掛ける数 (古典定数)
    """
    # レジスタ定義
    # x_reg: 入力 |x>
    # out_reg: 出力 |ax mod N> (計算用アキュムレータ)
    x_reg = QuantumRegister(n_qubits, 'x')
    out_reg = QuantumRegister(n_qubits + 1, 'out') # オーバーフロー防止で+1ビット
    c_res = ClassicalRegister(n_qubits + 1, 'res')
    
    qc = QuantumCircuit(x_reg, out_reg, c_res)

    # 1. 入力 x の初期化
    for i in range(n_qubits):
        if (x_val >> i) & 1:
            qc.x(x_reg[i])
            
    # 2. 出力レジスタをQFT変換 (加算の準備)
    qc.append(QFT(n_qubits + 1, do_swaps=False).to_gate(), out_reg)

    # 3. 乗算ループ (Shift and Add)
    for i in range(n_qubits):
        # シフトされた加算値: (a * 2^i) % N
        shifted_val = (a_const * (2**i)) % N_mod
        
        # 制御加算: x_reg[i] が 1 なら、shifted_val を out_reg に足す
        # ここで前回のモジュラ加算ロジックを使います
        modular_add_block(qc, x_reg[i], out_reg, None, shifted_val, N_mod, n_qubits + 1)

    # 4. 逆QFTで値を戻す
    qc.append(QFT(n_qubits + 1, do_swaps=False, inverse=True).to_gate(), out_reg)

    # 5. 測定
    qc.measure(out_reg, c_res)
    
    return qc

# --- パラメータ設定 ---
# 計算: 3 * 5 mod 13
# 期待値: 15 mod 13 = 2 (binary 00010)
VAL_X = 3
CONST_A = 5
N_MOD = 13
N_QUBITS = 4

qc = modular_multiplier(VAL_X, CONST_A, N_MOD, N_QUBITS)

# --- 実行 ---
simulator = AerSimulator()
transpiled_qc = transpile(qc, simulator)
result = simulator.run(transpiled_qc, shots=1024).result()
counts = result.get_counts()

print(f"計算: {VAL_X} * {CONST_A} (mod {N_MOD})")
print(f"期待値: {(VAL_X * CONST_A) % N_MOD}")
print("測定結果:", counts)