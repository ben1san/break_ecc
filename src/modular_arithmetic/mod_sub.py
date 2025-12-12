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

def modular_subtractor(subtrahend_val, minuend_val, N_mod, n_qubits):
    """
    |b> - |a> mod N を計算する回路
    subtrahend_val: 引く数 a (古典値)
    minuend_val: 引かれる数 b (量子状態の初期値)
    """
    q_reg = QuantumRegister(n_qubits + 1, 'q') 
    ancilla = QuantumRegister(1, 'ancilla')
    c_out = ClassicalRegister(n_qubits + 1, 'out')
    
    qc = QuantumCircuit(q_reg, ancilla, c_out)

    # 1. 初期化 (|b> をセット)
    for i in range(n_qubits):
        if (minuend_val >> i) & 1:
            qc.x(q_reg[i])

    # 2. QFT変換 (do_swaps=True)
    qc.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 3. 減算: |b> - a -> |b-a>
    # inverse=True で減算を行います
    add_in_fourier(qc, q_reg, subtrahend_val, n_qubits + 1, inverse=True)

    # --- ここからロジック変更 ---
    # 加算器にあった「試し引き」は不要です。
    # すでに b - a を計算したので、これが負かどうかチェックするだけです。

    # 4. アンダーフロー判定のためにIQFT
    qc.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    # 5. 修復判定
    # 最上位ビット(符号ビット)が 1 なら「負になった」ことを意味する
    qc.cx(q_reg[n_qubits], ancilla[0])

    # 6. 再びQFT空間へ (修復加算のため)
    qc.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 7. 条件付き修復: もし負(アンシラが1)なら、Nを足す
    # 負の数に N を足せば、正の剰余に戻ります (-2 + 7 = 5)
    controlled_add_in_fourier(qc, ancilla[0], q_reg, N_mod, n_qubits + 1)

    # 8. 最終的なIQFT
    qc.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    qc.measure(q_reg, c_out)
    
    return qc