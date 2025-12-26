from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from src.arithmetic import ModularArithmetic

def run_test(test_name, qc, output_reg, expected_mod, N):
    print(f"--- {test_name} ---")
    simulator = AerSimulator()
    qc_transpiled = transpile(qc, simulator)
    job = simulator.run(qc_transpiled, shots=100)
    result = job.result()
    counts = result.get_counts()
    
    # 最頻値を取得
    top_measurement = max(counts, key=counts.get)
    # 結果はビット文字列なので整数に変換 (空白除去)
    measured_int = int(top_measurement.replace(" ", ""), 2)
    
    print(f"Measured (Raw): {measured_int} (bin: {top_measurement})")
    print(f"Measured (Mod {N}): {measured_int % N}")
    print(f"Expected (Mod {N}): {expected_mod}")
    
    if measured_int % N == expected_mod:
        print("RESULT: PASS ✅\n")
    else:
        print("RESULT: FAIL ❌\n")

def test_modular_operations():
    p = 13
    n_in = 4
    n_out = 8 # 加算結果がNを超えても良いように大きめに取る
    
    arith = ModularArithmetic(p, n_out)

    # ==========================================
    # Test 1: Modular Square (3^2 mod 13 = 9)
    # ==========================================
    qr_src = QuantumRegister(n_in, 'src')
    qr_out = QuantumRegister(n_out, 'out')
    cr = ClassicalRegister(n_out, 'c')
    qc = QuantumCircuit(qr_src, qr_out, cr)
    
    # Input: 3 (0011)
    qc.x(qr_src[0])
    qc.x(qr_src[1])
    
    arith.modular_square(qc, qr_src, qr_out)
    qc.measure(qr_out, cr)
    
    run_test("Square: 3^2 -> 9", qc, qr_out, 9, p)

    # ==========================================
    # Test 2: Modular Scalar Mult (3 * 4 mod 13 = 12)
    # ==========================================
    qr_src = QuantumRegister(n_in, 'src')
    qr_out = QuantumRegister(n_out, 'out')
    cr = ClassicalRegister(n_out, 'c')
    qc = QuantumCircuit(qr_src, qr_out, cr)
    
    # Input: 3 (0011)
    qc.x(qr_src[0])
    qc.x(qr_src[1])
    
    # Scalar: 4
    arith.modular_scalar_mult(qc, qr_src, qr_out, 4)
    qc.measure(qr_out, cr)
    
    run_test("Scalar Mult: 3 * 4 -> 12", qc, qr_out, 12, p)
    
    # ==========================================
    # Test 3: Modular Subtraction (9 - 3 = 6)
    # ==========================================
    qr_src = QuantumRegister(n_in, 'src')
    qr_target = QuantumRegister(n_out, 'target') # Target starts with 9
    cr = ClassicalRegister(n_out, 'c')
    qc = QuantumCircuit(qr_src, qr_target, cr)
    
    # Src: 3 (0011)
    qc.x(qr_src[0])
    qc.x(qr_src[1])
    
    # Target: 9 (1001) - 初期値としてセット
    qc.x(qr_target[0])
    qc.x(qr_target[3])
    
    arith.modular_sub(qc, qr_src, qr_target)
    qc.measure(qr_target, cr)
    
    run_test("Sub: 9 - 3 -> 6", qc, qr_target, 6, p)

if __name__ == "__main__":
    test_modular_operations()