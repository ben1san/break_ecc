from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from src.arithmetic import ModularArithmetic
from src.ecc import QuantumECC, ScalarMultiplication

# --- パラメータ設定 ---
N = 13          
n_qubits = 4    
a = 0           
P_base = (11, 5) 

# --- インスタンス化 ---
arith = ModularArithmetic(N, n_qubits)
ecc = QuantumECC(arith)
scalar_mult = ScalarMultiplication(ecc, arith, a=a)

# --- 回路構築 ---
k_reg = QuantumRegister(1, 'k')
x_reg = QuantumRegister(n_qubits, 'x1')
y_reg = QuantumRegister(n_qubits, 'y1')
z_reg = QuantumRegister(n_qubits, 'z1')
anc_reg = QuantumRegister(8 * n_qubits, 'anc')
c_res = ClassicalRegister(2 * n_qubits, 'res')

qc = QuantumCircuit(k_reg, x_reg, y_reg, z_reg, anc_reg, c_res)

# --- 【修正】初期化 (Xゲートを使用) ---
def set_value(circuit, reg, val):
    """整数 val をビットパターンとしてレジスタにセットする"""
    for i in range(len(reg)):
        if (val >> i) & 1:
            circuit.x(reg[i])

# P(11, 5, 1) をセット
set_value(qc, x_reg, 11) # 1011
set_value(qc, y_reg, 5)  # 0101
set_value(qc, z_reg, 1)  # 0001

# 制御ビット k を |1> にする
qc.x(k_reg[0])

# --- 制御付き加算ゲートの適用 ---
reg_specs = {'p_size': n_qubits, 'anc_size': len(anc_reg)}
ctrl_add_gate = scalar_mult.create_controlled_add_gate(P_base, reg_specs)

qubits = [k_reg[0]] + list(x_reg) + list(y_reg) + list(z_reg) + list(anc_reg)
qc.append(ctrl_add_gate, qubits)

# --- 測定 ---
qc.measure(x_reg, c_res[:n_qubits])
qc.measure(y_reg, c_res[n_qubits:])

# --- 実行 ---
print("Simulating (MPS)...")
backend = AerSimulator(method='matrix_product_state')
t_qc = transpile(qc, backend, basis_gates=['u', 'cx', 'p', 'swap', 'x', 'id', 'measure'], optimization_level=1)
result = backend.run(t_qc, shots=1).result()
counts = result.get_counts()

# --- 結果確認 ---
expected_2P = scalar_mult._classical_point_doubling(P_base, N)
print(f"Expected 2P: {expected_2P}")

measured = list(counts.keys())[0]
val_y = int(measured[:n_qubits], 2)
val_x = int(measured[n_qubits:], 2)
print(f"Measured 2P: ({val_x}, {val_y})")

if (val_x, val_y) == expected_2P:
    print("\nSUCCESS!")
else:
    print("\nMISMATCH")