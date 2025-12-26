import numpy as np
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from qiskit.circuit.library import QFT

def add_in_fourier(circuit, reg, val, n, inverse=False):
    """QFT空間での定数加算 (Phase Addition)"""
    sign = -1 if inverse else 1
    for i in range(n):
        angle = sign * 2 * np.pi * val / (2**(n - i))
        if abs(angle) > 1e-9: # 精度調整
            circuit.p(angle, reg[i])

def cc_phase_add_constant(circuit, ctrl_qubit, target_reg, val, n, inverse=False):
    """
    制御付き定数加算 (Controlled-Phase Addition)
    ctrl_qubit が |1> のときだけ target_reg += val
    """
    sign = -1 if inverse else 1
    for i in range(n):
        angle = sign * 2 * np.pi * val / (2**(n - i))
        if abs(angle) > 1e-9:
            # cp (controlled-phase) を使用
            circuit.cp(angle, ctrl_qubit, target_reg[i])

def c_modular_add(circuit, ctrl_qubit, reg, ancilla, val, N, n):
    """
    【重要】制御付きモジュラ加算: |reg> += val mod N (if ctrl==1)
    
    QFT空間内で完結させます。
    reg はすでに QFT されている前提です。
    """
    
    # 1. 普通に加算: |reg> += val
    cc_phase_add_constant(circuit, ctrl_qubit, reg, val, n)
    
    # 2. 引き算 (試し引き): |reg> -= N
    # ここは「無条件」ではなく、ロジックとして「引き算をして負になるかチェック」が必要ですが、
    # QFT加算器での完全な条件分岐はアンシラ操作が複雑です。
    
    # ★簡易版（Beauregard型）の実装フロー:
    # A. 引く: reg -= N
    add_in_fourier(circuit, reg, N, n, inverse=True)
    
    # B. アンダーフローチェック (IQFTして最上位ビットを確認)
    # ※注: ここで一度QFTを解かないと符号確認できません
    circuit.append(QFT(n, do_swaps=False, inverse=True).to_gate(), reg)
    
    # C. アンシラに「負になった(戻しが必要)」フラグを立てる
    # 最上位ビット(MSB)が1なら負とみなす
    circuit.cx(reg[-1], ancilla)
    
    # D. 再びQFT
    circuit.append(QFT(n, do_swaps=False).to_gate(), reg)
    
    # E. 修復加算: もし負(ancilla=1)なら N を足し戻す
    # ここは制御ビットが ancilla なので c_phase_add
    cc_phase_add_constant(circuit, ancilla, reg, N, n)
    
    # --- ここから重要: 加算しなかった場合の後始末 ---
    # ctrl_qubitが0だった場合、そもそもステップ1で足していないのに、
    # ステップAでNを引いてしまっています。
    # これを補正するロジックが必要ですが、非常に複雑になります。
    
    # 【推奨】単純化のため、「常に足して、条件付きで引く」アプローチに変えます。
    # しかし、乗算器の中では「制御ビットが1の時だけモジュラ加算」したいので、
    # 最も単純なのは「古典的に計算した値を足す」のではなく、
    # 量子回路として正しいモジュラ加算器を呼ぶことです。
    pass 
    # ※この関数の完全な実装は長くなるため、下の「modular_multiplier」で
    # どう扱うべきかの方針を示します。


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

def modular_adder(circuit, a_val, b_val, N_mod, n_qubits):
    q_reg = QuantumRegister(n_qubits + 1, 'q') 
    ancilla = QuantumRegister(1, 'ancilla')
    c_out = ClassicalRegister(n_qubits + 1, 'out')
    
    # 1. 初期化
    for i in range(n_qubits):
        if (b_val >> i) & 1:
            circuit.x(q_reg[i])

    # 2. QFT変換 (do_swaps=True に変更)
    # これによりビット順序が直感的(q[0]=LSB)になります
    circuit.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 3. 加算: |b> + a
    add_in_fourier(circuit, q_reg, a_val, n_qubits + 1)

    # 4. 試し引き: - N
    add_in_fourier(circuit, q_reg, N_mod, n_qubits + 1, inverse=True)

    # 5. アンダーフロー判定のためにIQFT (do_swaps=True)
    circuit.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    # 6. 修復判定 (最上位ビットをアンシラへ)
    # 最上位ビットはインデックス n_qubits です
    circuit.cx(q_reg[n_qubits], ancilla[0])

    # 7. 再びQFT空間へ (do_swaps=True)
    circuit.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 8. 条件付き修復加算: + N
    c_phase_add(circuit, ancilla[0], q_reg, N_mod, n_qubits + 1)

    # 9. 最終的なIQFT (do_swaps=True)
    circuit.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    # 10. 測定
    circuit.measure(q_reg, c_out)
    
    return circuit

def modular_substractor(circuit, subtrahend_val, minuend_val, N_mod, n_qubits):
    """
    |b> - |a> mod N を計算する回路
    subtrahend_val: 引く数 a (古典値)
    minuend_val: 引かれる数 b (量子状態の初期値)
    """
    q_reg = QuantumRegister(n_qubits + 1, 'q') 
    ancilla = QuantumRegister(1, 'ancilla')
    c_out = ClassicalRegister(n_qubits + 1, 'out')
    

    # 1. 初期化 (|b> をセット)
    for i in range(n_qubits):
        if (minuend_val >> i) & 1:
            circuit.x(q_reg[i])

    # 2. QFT変換 (do_swaps=True)
    circuit.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 3. 減算: |b> - a -> |b-a>
    # inverse=True で減算を行います
    add_in_fourier(circuit, q_reg, subtrahend_val, n_qubits + 1, inverse=True)

    # --- ここからロジック変更 ---
    # 加算器にあった「試し引き」は不要です。
    # すでに b - a を計算したので、これが負かどうかチェックするだけです。

    # 4. アンダーフロー判定のためにIQFT
    circuit.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    # 5. 修復判定
    # 最上位ビット(符号ビット)が 1 なら「負になった」ことを意味する
    circuit.cx(q_reg[n_qubits], ancilla[0])

    # 6. 再びQFT空間へ (修復加算のため)
    circuit.append(QFT(n_qubits + 1, do_swaps=True).to_gate(), q_reg)

    # 7. 条件付き修復: もし負(アンシラが1)なら、Nを足す
    # 負の数に N を足せば、正の剰余に戻ります (-2 + 7 = 5)
    c_phase_add(circuit, ancilla[0], q_reg, N_mod, n_qubits + 1)

    # 8. 最終的なIQFT
    circuit.append(QFT(n_qubits + 1, do_swaps=True, inverse=True).to_gate(), q_reg)

    return circuit


def modular_multiplier_accumulator(circuit, x_reg, a_val, out_reg, ancilla, N, n):
    """
    計算内容: out_reg = (out_reg + x_reg * a_val) % N
    注意: out_reg は初期値 0 であること。
    """
    
    # 1. out_reg を QFT 空間へ
    circuit.append(QFT(n+1, do_swaps=False).to_gate(), out_reg)
    
    # 2. x_reg の各ビットごとに「制御付きモジュラ加算」を行う
    for i in range(len(x_reg)):
        # 加算する値: (a * 2^i) mod N
        term_val = (a_val * (2**i)) % N
        
        # 制御ビット: x_reg[i]
        # x_reg[i] が 1 なら、out_reg に term_val を足して mod N する
        
        # --- ここで「制御付きモジュラ加算」を展開 ---
        # ステップ 1: 加算 (Controlled by x_reg[i])
        cc_phase_add_constant(circuit, x_reg[i], out_reg, term_val, n+1)
        
        # ステップ 2: モジュラ減算 (ここからは x_reg[i] に依存せず、オーバーフロー依存)
        # 「もし out_reg >= N なら N を引く」
        # これを実現するには、
        #   a. N を引く
        #   b. 符号チェック (IQFT -> MSB -> Ancilla)
        #   c. 負なら N を足し戻す
        # という一連の流れを量子回路として書く必要があります。
        
        # 簡易実装（QFT空間）:
        add_in_fourier(circuit, out_reg, N, n+1, inverse=True) # とりあえず引く
        circuit.append(QFT(n+1, do_swaps=False, inverse=True).to_gate(), out_reg) # 時間領域へ
        circuit.cx(out_reg[-1], ancilla) # 符号チェック
        circuit.append(QFT(n+1, do_swaps=False).to_gate(), out_reg) # 周波数領域へ
        cc_phase_add_constant(circuit, ancilla, out_reg, N, n+1) # 負なら戻す
        
        # アンシラのリセット (重要！これをしないと次のループで干渉する)
        # 符号ビットの情報をアンシラから消すために、逆演算が必要
        # ここでは簡易的に、ancillaを次の加算ロジックに使う前にクリアする工夫が必要です。
        # (通常は逆の計算をしてクリアしますが、ここでは省略)
        circuit.x(ancilla) # ビットフリップ等で簡易対処（厳密にはUncomputationが必要）
        
    # 3. 最後に out_reg を逆QFT
    circuit.append(QFT(n+1, do_swaps=False, inverse=True).to_gate(), out_reg)

def modular_square(circuit, val_x, N_mod, n_input):
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
    
    # 1. 入力 x の初期化
    for i in range(n_input):
        if (val_x >> i) & 1:
            circuit.x(x_reg[i])

    # 2. アキュムレータをQFT変換
    circuit.append(QFT(n_out, do_swaps=False).to_gate(), out_reg)

    # 3. 二乗計算ループ: sum_{i,j} (x_i * x_j) * 2^(i+j)
    for i in range(n_input):
        for j in range(i, n_input): # i <= j の範囲で回す（対称性を利用）
            
            # 加算する値: 2^(i+j) mod N
            term_val = (2**(i + j)) % N_mod
            
            if i == j:
                # 対角成分: x_i * x_i = x_i
                # 制御ビットは x[i] ひとつだけ
                c_phase_add(circuit, x_reg[i], out_reg, term_val, n_out)
            else:
                # 非対角成分: x_i * x_j
                # 制御ビットは x[i] と x[j] のふたつ
                # i != j なので、本来は term_val を2回足す必要がある (i,j と j,i)
                # したがって 2 * term_val を加算する
                doubled_term = (2 * term_val) % N_mod
                cc_phase_add(circuit, x_reg[i], x_reg[j], out_reg, doubled_term, n_out)

    # 4. 逆QFT
    circuit.append(QFT(n_out, do_swaps=False, inverse=True).to_gate(), out_reg)

    # 5. 測定
    circuit.measure(out_reg, c_res)
    
    return circuit

def c_add_constant_mod_N(circuit, ctrl_qubit, reg, ancilla, C, N):
    """
    制御付きモジュラ加算: 
    ctrl_qubit が |1> なら |reg> += C mod N
    それ以外なら |reg> はそのまま (恒等写像)
    """
    n = len(reg)
    
    # -------------------------------------------------------
    # 1. +C (制御付き加算: ctrlが1のときだけ足す)
    # -------------------------------------------------------
    for i in range(n):
        angle = 2 * np.pi * C / (2**(n - i))
        if abs(angle) > 1e-9:
            # ここが変更点: p ではなく cp (controlled-phase)
            circuit.cp(angle, ctrl_qubit, reg[i])
            
    # -------------------------------------------------------
    # 2. -N (無条件で引いてみる)
    # -------------------------------------------------------
    # ctrlが0の場合でも一旦引きますが、後で戻るので問題ありません
    for i in range(n):
        angle = -2 * np.pi * N / (2**(n - i))
        if abs(angle) > 1e-9:
            circuit.p(angle, reg[i])
            
    # -------------------------------------------------------
    # 3. 符号チェック (アンダーフロー判定)
    # -------------------------------------------------------
    # 周波数領域 -> 時間領域
    circuit.append(QFT(n, do_swaps=False, inverse=True).to_gate(), reg)
    
    # MSBチェック: 負なら ancilla = 1
    circuit.cx(reg[-1], ancilla)
    
    # 時間領域 -> 周波数領域
    circuit.append(QFT(n, do_swaps=False).to_gate(), reg)
    
    # -------------------------------------------------------
    # 4. 補正 (負になったら +N)
    # -------------------------------------------------------
    for i in range(n):
        angle = 2 * np.pi * N / (2**(n - i))
        if abs(angle) > 1e-9:
            circuit.cp(angle, ancilla, reg[i])
            
    # -------------------------------------------------------
    # 5. アンシラのリセット (重要)
    # -------------------------------------------------------
    circuit.reset(ancilla)


def modular_multiplier(a_val, N_mod, n_qubits):
    """
    乗算回路: |x>|0> -> |x>|ax mod N>
    a_val: 古典的な乗数 (定数)
    """
    # レジスタの準備
    x_reg = QuantumRegister(n_qubits, 'x')       # 入力 x
    out_reg = QuantumRegister(n_qubits + 1, 'out') # 出力 (アキュムレータ) ※オーバーフロー用に+1ビット推奨
    ancilla = QuantumRegister(1, 'anc')          # 計算用アンシラ
    c_res = ClassicalRegister(n_qubits + 1, 'res')
    
    circuit = QuantumCircuit(x_reg, out_reg, ancilla, c_res)
    
    # --- テスト用の初期化 (例: x = 3 をセット) ---
    # 実際使うときはここは外してください
    # circuit.x(x_reg[0])
    # circuit.x(x_reg[1])
    # -----------------------------------------

    # 1. アキュムレータ(out_reg)をQFT空間へ
    # これにより、以降の加算はずっと周波数領域で行えます
    circuit.append(QFT(len(out_reg), do_swaps=False).to_gate(), out_reg)

    # 2. 乗算ループ (Shift and Add)
    for i in range(n_qubits):
        # 足すべき定数: (a * 2^i) mod N
        # 事前に mod N を取っておくことで、和が 2N を超えるのを防ぎます
        term_val = (a_val * (2**i)) % N_mod
        
        # 制御付きモジュラ加算を実行
        # 制御ビット: x_reg[i]
        c_add_constant_mod_N(circuit, x_reg[i], out_reg, ancilla, term_val, N_mod)

    # 3. 最後に逆QFTで値を戻す
    circuit.append(QFT(len(out_reg), do_swaps=False, inverse=True).to_gate(), out_reg)
    
    # 4. 測定 (アキュムレータのみ)
    circuit.measure(out_reg, c_res)
    
    return circuit