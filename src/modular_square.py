import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import QFT

def cc_phase_add(circuit, ctrl1, ctrl2, target_reg, val, n_target):
    """
    2重制御付き位相加算 (Controlled-Controlled-Phase Addition)
    ctrl1 と ctrl2 が両方 1 のときだけ、target_reg に val を足す
    """
    for k in range(n_target):
        # 加える角度: 2*pi * val / 2^(n-k)
        angle = 2 * np.pi * val / (2**(n_target - k))
        
        # CPゲート(制御位相)にさらにもう一つ制御を加えたもの -> MCPhase (Multi-Controlled Phase)
        # Qiskitでは cp(theta, c, t) ですが、2つ制御がある場合は mcp または cpをラップする
        if abs(angle) > 1e-5: # 最適化: 角度がほぼ0ならゲートを置かない
            circuit.mcphase(angle, [ctrl1, ctrl2], target_reg[k])

def c_phase_add(circuit, ctrl, target_reg, val, n_target):
    """
    単一制御付き位相加算 (対角成分用)
    """
    for k in range(n_target):
        angle = 2 * np.pi * val / (2**(n_target - k))
        if abs(angle) > 1e-5:
            circuit.cp(angle, ctrl, target_reg[k])

def modular_square(val_x, N_mod, n_input):
    """
    |x> -> |x^2 mod N> を計算する回路
    val_x: 入力の初期値（整数）
    """
    # ビット幅の計算: x^2 は最大 2*n ビットになるため、出力レジスタは大きめに取る
    n_out = 2 * n_input 
    if n_out > N_mod.bit_length() + 2: # Nより十分大きければNのビット数+余裕でOK
        n_out = N_mod.bit_length() + 2

    x_reg = QuantumRegister(n_input, 'x')
    out_reg = QuantumRegister(n_out, 'out') # アキュムレータ
    c_res = ClassicalRegister(n_out, 'res')
    
    qc = QuantumCircuit(x_reg, out_reg, c_res)

    # 1. 入力 x の初期化
    for i in range(n_input):
        if (val_x >> i) & 1:
            qc.x(x_reg[i])

    # 2. アキュムレータをQFT変換
    qc.append(QFT(n_out, do_swaps=False).to_gate(), out_reg)

    # 3. 二乗計算ループ: sum_{i,j} (x_i * x_j) * 2^(i+j)
    for i in range(n_input):
        for j in range(i, n_input): # i <= j の範囲で回す（対称性を利用）
            
            # 加算する値: 2^(i+j) mod N
            term_val = (2**(i + j)) % N_mod
            
            if i == j:
                # 対角成分: x_i * x_i = x_i
                # 制御ビットは x[i] ひとつだけ
                c_phase_add(qc, x_reg[i], out_reg, term_val, n_out)
            else:
                # 非対角成分: x_i * x_j
                # 制御ビットは x[i] と x[j] のふたつ
                # i != j なので、本来は term_val を2回足す必要がある (i,j と j,i)
                # したがって 2 * term_val を加算する
                doubled_term = (2 * term_val) % N_mod
                cc_phase_add(qc, x_reg[i], x_reg[j], out_reg, doubled_term, n_out)

    # 4. 逆QFT
    qc.append(QFT(n_out, do_swaps=False, inverse=True).to_gate(), out_reg)

    # 5. 測定
    qc.measure(out_reg, c_res)
    
    return qc

# --- パラメータ設定 ---
# 計算: 3^2 mod 13
# 期待値: 9
VAL_X = 3
N_MOD = 13
N_INPUT_QUBITS = 3 # 3を表現するのに十分なビット数

qc = modular_square(VAL_X, N_MOD, N_INPUT_QUBITS)

# --- 実行 ---
simulator = AerSimulator()
transpiled_qc = transpile(qc, simulator)
result = simulator.run(transpiled_qc, shots=1024).result()
counts = result.get_counts()

# --- 結果表示 (ビット列を整数に変換) ---
print(f"計算: {VAL_X}^2 (mod {N_MOD})")
# Qiskitの出力キーは 'bitstring' (例: '001001')
# これをパースして最も確率が高いものを表示
sorted_counts = sorted(counts.items(), key=lambda item: item[1], reverse=True)
most_frequent_bitstring = sorted_counts[0][0]
measured_int = int(most_frequent_bitstring, 2)

print(f"測定結果(Binary): {most_frequent_bitstring}")
print(f"測定結果(Int): {measured_int}")
print(f"正解: {(VAL_X**2) % N_MOD}")