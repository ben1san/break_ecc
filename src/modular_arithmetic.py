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
    c_phase_add(circuit, ctrl, reg, val, n)

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
    c_phase_add(qc, ancilla[0], q_reg, N_mod, n_qubits + 1)

    # 9. 最終的なIQFT (do_swaps=True)
    qc.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    # 10. 測定
    qc.measure(q_reg, c_out)
    
    return qc

def modular_substractor(subtrahend_val, minuend_val, N_mod, n_qubits):
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
    c_phase_add(qc, ancilla[0], q_reg, N_mod, n_qubits + 1)

    # 8. 最終的なIQFT
    qc.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    qc.measure(q_reg, c_out)
    
    return qc


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

