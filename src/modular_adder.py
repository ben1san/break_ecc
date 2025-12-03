import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import QFT

def add_in_fourier(circuit, reg, val, n, inverse=False):
    """
    QFT空間での定数加算/減算を行う関数
    inverse=True で減算になる
    """
    sign = -1 if inverse else 1
    for i in range(n):
        angle = sign * 2 * np.pi * val / (2**(n - i))
        circuit.p(angle, reg[i])

def controlled_add_in_fourier(circuit, ctrl_bit, reg, val, n):
    """
    制御ビット付きの定数加算 (もしctrlが1なら val を足す)
    """
    for i in range(n):
        angle = 2 * np.pi * val / (2**(n - i))
        circuit.cp(angle, ctrl_bit, reg[i])

def modular_adder(a_val, b_val, N_mod, n_qubits):
    """
    |a> + |b> mod N を計算する回路
    a_val: 足す数（古典値）
    b_val: 足される数（量子状態の初期値）
    N_mod: 法（モジュラス）
    """
    # レジスタ: 計算用(n+1ビット) + 判定用アンシラ(1ビット) + 出力
    # n+1ビットにする理由: オーバーフロー(符号)を確認するため
    q_reg = QuantumRegister(n_qubits + 1, 'q') 
    ancilla = QuantumRegister(1, 'ancilla')
    c_out = ClassicalRegister(n_qubits + 1, 'out')
    
    qc = QuantumCircuit(q_reg, ancilla, c_out)

    # 1. 初期化 (|b> をセット)
    # ここでは b_val をセットします
    for i in range(n_qubits):
        if (b_val >> i) & 1:
            qc.x(q_reg[i])
            
    # アンシラは |0> スタート、q_regの最上位ビットも |0>

    # 2. QFT変換 (加算の準備)
    qc.append(QFT(n_qubits + 1, do_swaps=False).to_gate(), q_reg)

    # 3. 加算: |b> + a -> |a+b>
    add_in_fourier(qc, q_reg, a_val, n_qubits + 1)

    # 4. 試し引き: |a+b> - N
    # もし結果が負になれば、最上位ビット(符号ビット)が1になる(2の補数表現的に)
    add_in_fourier(qc, q_reg, N_mod, n_qubits + 1, inverse=True)

    # 5. アンダーフロー判定のために一度IQFTで戻す
    # QFT空間ではビットの値を確認できないため
    qc.append(QFT(n_qubits + 1, do_swaps=False, inverse=True).to_gate(), q_reg)

    # 6. 修復判定
    # 最上位ビット(q_reg[n])が 1 なら「負になった（引きすぎた）」ことを意味する
    # この情報をアンシラにコピー (CNOT)
    qc.cx(q_reg[n_qubits], ancilla[0])

    # 7. 再びQFT空間へ (修復加算のため)
    qc.append(QFT(n_qubits + 1, do_swaps=False).to_gate(), q_reg)

    # 8. 条件付き修復加算: もしアンシラが1なら、Nを足し戻す
    controlled_add_in_fourier(qc, ancilla[0], q_reg, N_mod, n_qubits + 1)

    # 9. 最終的なIQFT
    qc.append(QFT(n_qubits + 1, do_swaps=False, inverse=True).to_gate(), q_reg)

    # --- 重要: アンシラのアンコンピュテーション（今回は省略） ---
    # 本来はここでアンシラを|0>に戻す処理が必要ですが、
    # 測定して終わりの場合は省略可能です。

    # 10. 測定
    qc.measure(q_reg, c_out)
    
    return qc

# --- パラメータ設定 ---
# 5 + 6 = 11.  11 mod 7 = 4.
# 期待される出力: 4 (binary 00100) ※最上位ビットは符号用なので0
N_MOD = 7
VAL_A = 5
VAL_B = 6
N_QUBITS = 4 # 数値を表現するのに十分なビット数

qc = modular_adder(VAL_A, VAL_B, N_MOD, N_QUBITS)

# --- 実行 ---
simulator = AerSimulator()
transpiled_qc = transpile(qc, simulator)
result = simulator.run(transpiled_qc, shots=1024).result()
counts = result.get_counts()

print(f"計算: {VAL_A} + {VAL_B} (mod {N_MOD})")
print(f"期待値: {(VAL_A + VAL_B) % N_MOD}")
print("測定結果:", counts)