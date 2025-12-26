from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile # ← transpileを追加
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
from arithmetic import ModularArithmetic

def test_square_and_uncompute():
    # パラメータ設定 (Bit size 4)
    n_in = 4
    n_out = 8 
    
    # レジスタ準備
    qr_src = QuantumRegister(n_in, 'z')
    qr_out = QuantumRegister(n_out, 'z_sq')
    cr_src = ClassicalRegister(n_in, 'c_z')
    cr_out = ClassicalRegister(n_out, 'c_sq')
    
    # 回路作成
    qc = QuantumCircuit(qr_src, qr_out, cr_src, cr_out)
    arith = ModularArithmetic()

    # 1. 初期化: Z = 3 (0011)
    print("Initializing Z = 3...")
    qc.x(qr_src[0])
    qc.x(qr_src[1])
    qc.barrier()

    # 2. 計算: Z^2 を計算
    print("Computing Z^2...")
    arith.compute_square(qc, qr_src, qr_out)
    qc.barrier()

    # 3. アンコンピュテーション: Z^2 を元に戻す
    print("Uncomputing Z^2...")
    arith.uncompute_square(qc, qr_src, qr_out)
    qc.barrier()

    # 4. 測定
    qc.measure(qr_src, cr_src)
    qc.measure(qr_out, cr_out)

    # 5. シミュレーション実行
    simulator = AerSimulator()
    
    # ★修正箇所: シミュレータ用に回路を変換(トランスパイル)する
    print("Transpiling circuit...")
    qc_transpiled = transpile(qc, simulator) 
    
    # 変換した回路を実行する
    job = simulator.run(qc_transpiled, shots=1000)
    result = job.result()
    counts = result.get_counts()

    print("\nResult Counts (c_sq c_z):")
    # 期待値: c_sq=00000000, c_z=0011 (つまり '00000000 0011')
    print(counts)
    
    # 結果の検証
    expected_key = "00000000 0011"
    if counts.get(expected_key, 0) > 990:
        print("\nSUCCESS: Uncomputation worked perfectly! Output register is clean.")
    else:
        print("\nFAILURE: Output register was not reset.")

if __name__ == "__main__":
    test_square_and_uncompute()