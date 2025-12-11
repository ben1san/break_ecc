import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import QFT

def add_in_fourier(circuit, reg, val, n, inverse=False):
    """
    QFT空間での定数加算/減算を行う関数
    """
    sign = -1 if inverse else 1
    for i in range(n):
        # q[0]がLSB、q[n-1]がMSBである前提の角度計算
        angle = sign * 2 * np.pi * val / (2**(n - i))
        circuit.p(angle, reg[i])

def controlled_add_in_fourier(circuit, ctrl_bit, reg, val, n):
    """
    制御ビット付きの定数加算
    """
    for i in range(n):
        angle = 2 * np.pi * val / (2**(n - i))
        circuit.cp(angle, ctrl_bit, reg[i])

def modular_adder(a_val, b_val, N_mod, n_qubits):
    q_reg = QuantumRegister(n_qubits + 1, 'q') 
    ancilla = QuantumRegister(1, 'ancilla')
    c_out = ClassicalRegister(n_qubits + 1, 'out')
    
    qc = QuantumCircuit(q_reg, ancilla, c_out)

    # 1. 初期化
    for i in range(n_qubits):
        if (b_val >> i) & 1:
            qc.x(q_reg[i])

    # 2. QFT変換 (do_swaps=True に変更)
    # これによりビット順序が直感的(q[0]=LSB)になります
    qc.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 3. 加算: |b> + a
    add_in_fourier(qc, q_reg, a_val, n_qubits + 1)

    # 4. 試し引き: - N
    add_in_fourier(qc, q_reg, N_mod, n_qubits + 1, inverse=True)

    # 5. アンダーフロー判定のためにIQFT (do_swaps=True)
    qc.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    # 6. 修復判定 (最上位ビットをアンシラへ)
    # 最上位ビットはインデックス n_qubits です
    qc.cx(q_reg[n_qubits], ancilla[0])

    # 7. 再びQFT空間へ (do_swaps=True)
    qc.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 8. 条件付き修復加算: + N
    controlled_add_in_fourier(qc, ancilla[0], q_reg, N_mod, n_qubits + 1)

    # 9. 最終的なIQFT (do_swaps=True)
    qc.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    # 10. 測定
    qc.measure(q_reg, c_out)
    
    return qc

# --- パラメータ設定 ---
N_MOD = 7
VAL_A = 5
VAL_B = 6
N_QUBITS = 4 

qc = modular_adder(VAL_A, VAL_B, N_MOD, N_QUBITS)

# --- 実行 ---
simulator = AerSimulator()
transpiled_qc = transpile(qc, simulator)
result = simulator.run(transpiled_qc, shots=1024).result()
counts = result.get_counts()

print(f"計算: {VAL_A} + {VAL_B} (mod {N_MOD})")
print(f"期待値: {(VAL_A + VAL_B) % N_MOD}")
# 見やすいようにソートして表示
sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
print("測定結果(Top 5):", sorted_counts[:5])