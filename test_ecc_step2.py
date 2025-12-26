from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from src.arithmetic import ModularArithmetic
from src.ecc import QuantumECC

def test_Z3_and_cleanup():
    p = 13
    n_in = 4
    
    # --- レジスタサイズ調整 ---
    # T1: Z3用。H(5bit) * Z1(4bit) の結果を受けるため 5bit に拡張
    n_t1 = 5 
    n_t2 = 5 # H用
    n_t3 = 4 # Z^3用 (ここは mod p されるので 4bit で足りる)
    n_t4 = 5 # R用
    
    print(f"--- Step 2: Z3 Calc & Uncomputation (Optimized) ---")

    # arithmeticのビット数はデフォルト値。個別のレジスタサイズが優先される
    arith = ModularArithmetic(p, n_t1)
    ecc = QuantumECC(arith)
    
    qr_X1 = QuantumRegister(n_in, 'X1')
    qr_Y1 = QuantumRegister(n_in, 'Y1')
    qr_Z1 = QuantumRegister(n_in, 'Z1')
    
    qr_T1 = QuantumRegister(n_t1, 'T1_Z3') 
    qr_T2 = QuantumRegister(n_t2, 'T2_H')
    qr_T3 = QuantumRegister(n_t3, 'T3_Clean')
    qr_T4 = QuantumRegister(n_t4, 'T4_R')
    
    cr_Z3 = ClassicalRegister(n_t1, 'c_Z3')
    cr_clean = ClassicalRegister(n_t3, 'c_clean')
    
    qc = QuantumCircuit(qr_X1, qr_Y1, qr_Z1, qr_T1, qr_T2, qr_T3, qr_T4, cr_Z3, cr_clean)
    
    print(f"Total Qubits: {qc.num_qubits} (Limit check: 30+ is OK with MPS)")

    # 初期化: P1 = (4, 8, 1)
    qc.x(qr_X1[2]) # 4
    qc.x(qr_Y1[3]) # 8
    qc.x(qr_Z1[0]) # 1
    
    const_P2 = (11, 5)
    
    # Step 1: H, R 計算
    print("1. Calculating H and R...")
    ecc.calculate_H_R(qc, (qr_X1, qr_Y1, qr_Z1), const_P2, [qr_T1, qr_T2, qr_T3, qr_T4])
    
    # Step 2: Z3 計算 & お掃除
    print("2. Uncomputing temps and calculating Z3...")
    ecc.calculate_Z3_and_cleanup(qc, (qr_X1, qr_Y1, qr_Z1), [qr_T1, qr_T2, qr_T3, qr_T4])
    
    # 測定
    qc.measure(qr_T1, cr_Z3)
    qc.measure(qr_T3, cr_clean)
    
    # シミュレーション
    print("Starting Simulation...")
    simulator = AerSimulator(method='matrix_product_state') 
    
    # coupling_map=None でビット数制限を回避
    qc_transpiled = transpile(qc, simulator, coupling_map=None)
    
    job = simulator.run(qc_transpiled, shots=100)
    counts = job.result().get_counts()
    
    top_measurement = max(counts, key=counts.get)
    # ビット数に合わせて分割位置を調整 (T3=4bit, T1=5bit)
    # 文字列は "c_clean c_Z3" の順 (右側が T1) ではなく、レジスタ定義順による
    # Qiskitの並び: cr_clean(4) cr_Z3(5) の順で出力されるはず (スペース区切り)
    parts = top_measurement.split()
    res_clean_bin = parts[0]
    res_Z3_bin = parts[1]
    
    val_Z3 = int(res_Z3_bin, 2)
    val_clean = int(res_clean_bin, 2)
    
    print(f"\nMeasured Z3: {val_Z3} (Mod 13: {val_Z3 % p}) | Expected: 7")
    print(f"Cleanup Check (T3): {val_clean} (Mod 13: {val_clean % p}) | Expected: 0")
    
    # 判定条件: Mod 13 で一致すればOK
    pass_z3 = (val_Z3 % p == 7)
    pass_clean = (val_clean % p == 0)
    
    if pass_z3 and pass_clean:
        print("RESULT: PASS ✅")
    else:
        print("RESULT: FAIL ❌")

if __name__ == "__main__":
    test_Z3_and_cleanup()