from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from general.arithmetic import ModularArithmetic
from general.ecc import QuantumECC, ScalarMultiplication

# --- パラメータ設定 (軽量テスト用) ---
# N=5, Curve: y^2 = x^3 + x + 1 (mod 5) -> a=1, b=1
N = 5
n_qubits = 12   # 12ビットあれば N=5 の数ステップの爆発には耐えられる
a = 1 
b = 1
P_base = (0, 1) # ベースポイント

# テスト: 2P からスタートし、k=3 (011) を加算
# Start: 2P = (4, 2)
# Add k=3 (binary 011):
#   Bit 0 (1): Add 1P -> 3P = (2, 1)
#   Bit 1 (1): Add 2P -> 5P = (0, 4) ... result
#   Bit 2 (0): Skip
start_point_idx = 2 
k_val = 3           
k_bits = 3

# --- インスタンス化 ---
arith = ModularArithmetic(N, n_qubits)
ecc = QuantumECC(arith)
scalar_mult = ScalarMultiplication(ecc, arith, a=a)

# --- 古典計算で期待値を算出 ---
# Python上で正解を計算しておく
start_point = scalar_mult._classical_scalar_mult(start_point_idx, P_base, N)
expected_idx = start_point_idx + k_val
expected_point = scalar_mult._classical_scalar_mult(expected_idx, P_base, N)

print(f"--- Small N Test (N={N}) ---")
print(f"Start Point (2P): {start_point}")
print(f"Target Point (5P): {expected_point}")

# --- 回路構築 ---
# ビット数を減らした構成
k_reg = QuantumRegister(k_bits, 'k')
x_reg = QuantumRegister(n_qubits, 'x1')
y_reg = QuantumRegister(n_qubits, 'y1')
z_reg = QuantumRegister(n_qubits, 'z1')
# アンシラも n=12 なので 96ビット程度に収まる
anc_reg = QuantumRegister(8 * n_qubits, 'anc') 
c_res = ClassicalRegister(3 * n_qubits, 'res')

qc = QuantumCircuit(k_reg, x_reg, y_reg, z_reg, anc_reg, c_res)

# --- 初期化 ---
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
print(f"Building scalar loop with {n_qubits} qubits...")
scalar_mult.build_scalar_mult_circuit(qc, k_reg, [x_reg, y_reg, z_reg], list(anc_reg), P_base, N)

# --- 測定 ---
qc.measure(x_reg, c_res[0:n_qubits])
qc.measure(y_reg, c_res[n_qubits:2*n_qubits])
qc.measure(z_reg, c_res[2*n_qubits:3*n_qubits])

# --- 実行 ---
print("Simulating (MPS)... Should be fast now.")
backend = AerSimulator(method='matrix_product_state')
# 最適化レベル1で実行
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

print(f"Measured Projective (mod {N}): ({X_proj}, {Y_proj}, {Z_proj})")

if Z_proj == 0:
    affine_res = (None, None)
    print("Result is Point at Infinity")
else:
    inv_z = pow(Z_proj, -1, N)
    x_affine = (X_proj * (inv_z**2)) % N
    y_affine = (Y_proj * (inv_z**3)) % N
    affine_res = (x_affine, y_affine)

print(f"Measured Affine: {affine_res}")
print(f"Expected Affine: {expected_point}")

if affine_res == expected_point:
    print("\nSUCCESS! Loop logic verified with small N.")
else:
    print("\nMISMATCH")