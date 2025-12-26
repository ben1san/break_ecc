from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from src.arithmetic import ModularArithmetic
from src.ecc import QuantumECC

def test_full_addition_and_cleanup():
    p = 13
    n_in = 4   # 入力は mod 13 なので 4bit
    n_anc = 12  # 中間計算はオーバーフロー対策で少し余裕を持つ (5bit)

    print(f"--- Project: ECC Mixed Addition (Full Step) ---")
    print(f"Target: Calculate (X3, Y3, Z3) and clean up ancillas.")

    # --- 1. インスタンス生成 ---
    # デフォルトのビット数は n_anc に合わせますが、各レジスタ定義が優先されます
    arith = ModularArithmetic(p, n_anc)
    ecc = QuantumECC(arith)

    # --- 2. レジスタ定義 ---
    # 入力座標 P1 = (X1, Y1, Z1)
    qr_X1 = QuantumRegister(n_in, 'x1')
    qr_Y1 = QuantumRegister(n_in, 'y1')
    qr_Z1 = QuantumRegister(n_in, 'z1')

    # アンシラ (T1~T8)
    # マッピング:
    # T1: Z3 (出力)
    # T2: H -> Clean
    # T3: H^2 -> Clean
    # T4: R -> Clean
    # T5: V -> Clean
    # T6: X3 (出力)
    # T7: Y3 (出力)
    # T8: tmp -> Clean
    
    # 出力用(X3, Y3, Z3)は n_anc(5bit) あると安心ですが、値は mod p なので下位4bitが重要
    qr_T1 = QuantumRegister(n_anc, 't1_z3') 
    qr_T2 = QuantumRegister(n_anc, 't2_h')
    qr_T3 = QuantumRegister(n_anc, 't3_h2')
    qr_T4 = QuantumRegister(n_anc, 't4_r')
    qr_T5 = QuantumRegister(n_anc, 't5_v')
    qr_T6 = QuantumRegister(n_anc, 't6_x3')
    qr_T7 = QuantumRegister(n_anc, 't7_y3')
    qr_T8 = QuantumRegister(n_anc, 't8_tmp')

    # アンシラリストを作成
    ancillas = [qr_T1, qr_T2, qr_T3, qr_T4, qr_T5, qr_T6, qr_T7, qr_T8]

    # 測定用クラシカルレジスタ
    cr_res = ClassicalRegister(n_anc * 3, 'res_coords') # X3, Y3, Z3
    cr_cln = ClassicalRegister(n_anc * 5, 'res_clean')  # T2,3,4,5,8 (0になるべき)

    qc = QuantumCircuit(qr_X1, qr_Y1, qr_Z1, *ancillas, cr_res, cr_cln)

    print(f"Total Qubits: {qc.num_qubits} (MPS simulation required)")

    # --- 3. 初期化 ---
    # P1 = (4, 8, 1)  (Affine: (4, 8))
    # P2 = (11, 5)    (Affine: (11, 5))
    # 期待される結果: P1 + P2
    # 計算過程 (Mod 13):
    #   Affine加算: slope = (5-8)/(11-4) = -3/7 = 10*2 = 20 = 7
    #   x3 = 7^2 - 4 - 11 = 49 - 15 = 34 = 8
    #   y3 = 7(4-8) - 8 = 7(-4) - 8 = -28 - 8 = -36 = 3
    #   => Affine (8, 3)
    #
    # Projective計算 (コード上の期待値):
    #   Z3 = 7
    #   X3 = 2  (Affine変換: 2 / 7^2 = 2/10 = 8. OK)
    #   Y3 = 2  (Affine変換: 2 / 7^3 = 2/5 = 16 = 3. OK)
    
    # P1 セット
    qc.x(qr_X1[2]) # 4 (100)
    qc.x(qr_Y1[3]) # 8 (1000)
    qc.x(qr_Z1[0]) # 1 (001)

    const_P2 = (11, 5)

    # --- 4. 量子回路構築 (全ステップ) ---
    print("Building Circuit...")
    
    # Step 1: H, R
    ecc.calculate_H_R(qc, (qr_X1, qr_Y1, qr_Z1), const_P2, ancillas)
    
    # Step 2: Z3 & Uncompute T3, T1
    ecc.calculate_Z3_and_cleanup(qc, (qr_X1, qr_Y1, qr_Z1), ancillas)
    
    # Step 3: X3, Y3 & Final Cleanup
    # (前回実装した新しいメソッド)
    ecc.calculate_X3_Y3_and_final_cleanup(qc, (qr_X1, qr_Y1, qr_Z1), const_P2, ancillas)

# --- 5. 測定 ---
    # スライス幅を 5 ではなく n_anc (12) に合わせる
    
    # 出力座標 (X3, Y3, Z3) -> cr_res (size: 3 * n_anc)
    # Qiskitの並び順: T7(Y3), T6(X3), T1(Z3)
    qc.measure(qr_T7, cr_res[0:n_anc])             # Y3
    qc.measure(qr_T6, cr_res[n_anc : 2*n_anc])     # X3
    qc.measure(qr_T1, cr_res[2*n_anc : 3*n_anc])   # Z3

    # お掃除チェック -> cr_cln (size: 5 * n_anc)
    qc.measure(qr_T2, cr_cln[0:n_anc])
    qc.measure(qr_T3, cr_cln[n_anc : 2*n_anc])
    qc.measure(qr_T4, cr_cln[2*n_anc : 3*n_anc])
    qc.measure(qr_T5, cr_cln[3*n_anc : 4*n_anc])
    qc.measure(qr_T8, cr_cln[4*n_anc : 5*n_anc])

# --- 6. シミュレーション ---
    print("Starting Simulation (MPS)...")
    simulator = AerSimulator(method='matrix_product_state')
    
    # 【修正点】
    # backend=simulator を渡すとビット数制限(63)に引っかかるため削除します。
    # 代わりに basis_gates を指定して、ゲートを分解させます。
    # MPSがサポートする基本ゲート: u, cx, p, swap, x, id, measure など
    basis_gates = ['u', 'cx', 'p', 'swap', 'x', 'id', 'measure', 'rz', 'sx']
    
    qc_transpiled = transpile(
        qc, 
        basis_gates=basis_gates, 
        optimization_level=0
    )
    
    # 実行
    job = simulator.run(qc_transpiled, shots=100)
    counts = job.result().get_counts()
    
    # 最頻値を取得
    top_meas = max(counts, key=counts.get)
    
    # ビット列のパース (Qiskitは逆順出力: cln res の順)
    # 定義順: cr_res(15), cr_cln(25)
    # 出力文字列: "cr_cln_bits cr_res_bits" (スペース区切り)
    parts = top_meas.split()
    bin_cln = parts[0]
    bin_res = parts[1]

    # 結果の数値化 (スライス幅を n_anc に変更)
    val_z3 = int(bin_res[0 : n_anc], 2)
    val_x3 = int(bin_res[n_anc : 2*n_anc], 2)
    val_y3 = int(bin_res[2*n_anc : 3*n_anc], 2)
    
    val_clean_total = int(bin_cln, 2)

    # --- 7. 検証 ---
    expected_x3 = 2
    expected_y3 = 2
    expected_z3 = 7

    print(f"\n=== RESULTS ===")
    print(f"Coordinates (Mod 13):")
    print(f"  X3: {val_x3 % p} (Expected: {expected_x3}) -> {'OK' if (val_x3 % p == expected_x3) else 'NG'}")
    print(f"  Y3: {val_y3 % p} (Expected: {expected_y3}) -> {'OK' if (val_y3 % p == expected_y3) else 'NG'}")
    print(f"  Z3: {val_z3 % p} (Expected: {expected_z3}) -> {'OK' if (val_z3 % p == expected_z3) else 'NG'}")
    
    print(f"\nCleanup Verification:")
    print(f"  Sum of Garbage Regs: {val_clean_total} (Expected: 0)")
    
    if (val_x3 % p == expected_x3) and \
       (val_y3 % p == expected_y3) and \
       (val_z3 % p == expected_z3) and \
       (val_clean_total == 0):
        print("\n✅ TEST PASSED: Full Mixed Addition Successful")
    else:
        print("\n❌ TEST FAILED")

if __name__ == "__main__":
    test_full_addition_and_cleanup()