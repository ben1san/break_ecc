from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from src.arithmetic import ModularArithmetic
from src.ecc import QuantumECC, ScalarMultiplication

# --- パラメータ設定 ---
N = 13
# 【修正】中間値の爆発(N^6程度)に耐えるため、ビット数を24に増やす
n_qubits = 24 
a = 0           
P_base = (11, 5) 

# テスト: 2P からスタートし、k=3 (011) を加算
# 期待される動作:
# 1. Start: 2P
# 2. Bit 0 (1): Add 1P -> 3P
# 3. Bit 1 (1): Add 2P -> 5P (Result)
start_point_idx = 2 
k_val = 3           

# --- インスタンス化 ---
arith = ModularArithmetic(N, n_qubits)
ecc = QuantumECC(arith)
scalar_mult = ScalarMultiplication(ecc, arith, a=a)

# --- 古典計算での期待値 ---
start_point = scalar_mult._classical_scalar_mult(start_point_idx, P_base, N)
expected_idx = start_point_idx + k_val
expected_point = scalar_mult._classical_scalar_mult(expected_idx, P_base, N)
print(f"Start Point: {start_point}")
print(f"Target Point: {expected_point}")

# --- 回路構築 ---
k_reg = QuantumRegister(3, 'k')
x_reg = QuantumRegister(n_qubits, 'x1')
y_reg = QuantumRegister(n_qubits, 'y1')
z_reg = QuantumRegister(n_qubits, 'z1')
anc_reg = QuantumRegister(8 * n_qubits, 'anc')
c_res = ClassicalRegister(3 * n_qubits, 'res')

qc = QuantumCircuit(k_reg, x_reg, y_reg, z_reg, anc_reg, c_res)

# --- 初期化 (2P) ---
def set_value(circuit, reg, val):
    for i in range(len(reg)):
        if (val >> i) & 1: circuit.x(reg[i])

set_value(qc, x_reg, start_point[0])
set_value(qc, y_reg, start_point[1])
set_value(qc, z_reg, 1)

# --- スカラー k=3 (011) ---
qc.x(k_reg[0])
qc.x(k_reg[1])

# --- ループ構築 ---
print(f"Building scalar loop with {n_qubits} qubits (this may be large)...")
scalar_mult.build_scalar_mult_circuit(qc, k_reg, [x_reg, y_reg, z_reg], list(anc_reg), P_base, N)

# --- 測定 ---
qc.measure(x_reg, c_res[0:n_qubits])
qc.measure(y_reg, c_res[n_qubits:2*n_qubits])
qc.measure(z_reg, c_res[2*n_qubits:3*n_qubits])

# --- 実行 ---
print("Simulating (MPS)... might take a minute due to large registers.")
backend = AerSimulator(method='matrix_product_state')
# 最適化レベルを1にしてトランスパイル時間を短縮
t_qc = transpile(qc, basis_gates=['u', 'cx', 'p', 'swap', 'x', 'id', 'measure'], optimization_level=1)
result = backend.run(t_qc, shots=1).result()
counts = result.get_counts()

# --- 結果検証 ---
measured_hex = list(counts.keys())[0]
val_z = int(measured_hex[0 : n_qubits], 2)
val_y = int(measured_hex[n_qubits : 2*n_qubits], 2)
val_x = int(measured_hex[2*n_qubits : 3*n_qubits], 2)

X_proj = val_x % N
Y_proj = val_y % N
Z_proj = val_z % N

print(f"Measured Projective: ({X_proj}, {Y_proj}, {Z_proj})")

if Z_proj == 0:
    affine_res = (None, None)
else:
    inv_z = pow(Z_proj, -1, N)
    x_affine = (X_proj * (inv_z**2)) % N
    y_affine = (Y_proj * (inv_z**3)) % N
    affine_res = (x_affine, y_affine)

print(f"Expected: {expected_point}")
print(f"Measured: {affine_res}")

if affine_res == expected_point:
    print("\nSUCCESS! Loop handled multiple additions without overflow.")
else:
    print("\nMISMATCH")