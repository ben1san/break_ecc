from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from src.arithmetic import ModularArithmetic
from src.ecc import QuantumECC

def test_H_R_calculation():
    p = 13
    n_in = 4
    
    # --- ビット数最適化 ---
    # 中間値 (Z^2, Z^3) は mod 13 済みなので 4bit でOK
    # 加減算結果 (H, R) は一時的に 25 (approx 2p) まで行くので 5bit 必要
    n_temp_small = 4
    n_temp_large = 5
    
    print(f"--- Setting up with Optimized Registers (4/5 bits) ---")

    # クラス初期化時の n_qubits はデフォルト値として使うだけなので、
    # メソッド呼び出し側でサイズを変えれば動的に対応します。
    arith = ModularArithmetic(p, n_temp_large)
    ecc = QuantumECC(arith)
    
    # 入力レジスタ (P1)
    qr_X1 = QuantumRegister(n_in, 'X1')
    qr_Y1 = QuantumRegister(n_in, 'Y1')
    qr_Z1 = QuantumRegister(n_in, 'Z1')
    
    # 作業用 (サイズを混在させる)
    qr_T1 = QuantumRegister(n_temp_small, 'T1_Zsq') # Z^2
    qr_T2 = QuantumRegister(n_temp_large, 'T2_H')   # H (要5bit)
    qr_T3 = QuantumRegister(n_temp_small, 'T3_Zcb') # Z^3
    qr_T4 = QuantumRegister(n_temp_large, 'T4_R')   # R (要5bit)
    
    # 測定用
    cr_H = ClassicalRegister(n_temp_large, 'c_H')
    cr_R = ClassicalRegister(n_temp_large, 'c_R')
    
    qc = QuantumCircuit(qr_X1, qr_Y1, qr_Z1, qr_T1, qr_T2, qr_T3, qr_T4, cr_H, cr_R)
    
    print(f"Total Qubits in Circuit: {qc.num_qubits}") # Should be 30

    # 初期化: P1 = (4, 8, 1)
    qc.x(qr_X1[2]) # 4
    qc.x(qr_Y1[3]) # 8
    qc.x(qr_Z1[0]) # 1
    
    const_P2 = (11, 5)
    
    # 計算実行
    print("Building circuit instructions...")
    ecc.calculate_H_R(qc, (qr_X1, qr_Y1, qr_Z1), const_P2, [qr_T1, qr_T2, qr_T3, qr_T4])
    
    # 測定
    qc.measure(qr_T2, cr_H)
    qc.measure(qr_T4, cr_R)
    
    # シミュレーション
    print("Starting Simulation...")
    simulator = AerSimulator(method='matrix_product_state') 
    
    qc_transpiled = transpile(qc, simulator)
    job = simulator.run(qc_transpiled, shots=100)
    counts = job.result().get_counts()
    
    # 結果解析
    top_measurement = max(counts, key=counts.get)
    res_R_bin, res_H_bin = top_measurement.split()
    
    val_H = int(res_H_bin, 2)
    val_R = int(res_R_bin, 2)
    
    print(f"\nMeasured H: {val_H} (Mod 13: {val_H % p}) | Expected: 7")
    print(f"Measured R: {val_R} (Mod 13: {val_R % p}) | Expected: 10 (-3 mod 13)")
    
    if (val_H % p == 7) and (val_R % p == 10):
        print("RESULT: PASS ✅")
    else:
        print("RESULT: FAIL ❌")

if __name__ == "__main__":
    test_H_R_calculation()