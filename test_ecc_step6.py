from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from src.arithmetic import ModularArithmetic
from src.ecc import QuantumECC, ScalarMultiplication

# --- パラメータ ---
N = 13          
n_qubits = 8    # オーバーフロー防止
a = 0           
P = (11, 5)     

# --- インスタンス ---
arith = ModularArithmetic(N, n_qubits)
ecc = QuantumECC(arith)
scalar_mult = ScalarMultiplication(ecc, arith, a=a)

# --- 回路構築 ---
# スカラーkは2ビット用意します (k=2 としたいところですが、ループ検証のため)
# 今回は手動構成で「P + (k=1のとき2P)」をテストします
k_reg = QuantumRegister(1, 'k_bit1') # 2^1 の位
x_reg = QuantumRegister(n_qubits, 'x')
y_reg = QuantumRegister(n_qubits, 'y')
z_reg = QuantumRegister(n_qubits, 'z')
anc_reg = QuantumRegister(8 * n_qubits, 'anc')
c_res = ClassicalRegister(3 * n_qubits, 'res')

qc = QuantumCircuit(k_reg, x_reg, y_reg, z_reg, anc_reg, c_res)

# --- 初期化 ---
def set_value(circuit, reg, val):
    for i in range(len(reg)):
        if (val >> i) & 1: circuit.x(reg[i])

# 1. アキュムレータを P に初期化 (これは k=1 の寄与分とみなす)
set_value(qc, x_reg, 11)
set_value(qc, y_reg, 5)
set_value(qc, z_reg, 1)

# 2. 制御ビットを |1> にセット (2Pを加算させるため)
qc.x(k_reg[0]) 

# --- ループ処理の模倣 ---
# build_scalar_mult_circuit を使う代わりに、手動で「2ビット目」の加算を行うイメージ
# Bit 0 (1P): 初期化で済ませた
# Bit 1 (2P): ここで加算する
Q_val = scalar_mult._classical_point_doubling(P, N) # これは (7, 5) になるはず
print(f"Preparing to add 2P: {Q_val}")

reg_specs = {'p_size': n_qubits, 'anc_size': len(anc_reg)}
ctrl_add_gate = scalar_mult.create_controlled_add_gate(Q_val, reg_specs)

# 適用: k_reg[0] が制御ビット
qubits = [k_reg[0]] + list(x_reg) + list(y_reg) + list(z_reg) + list(anc_reg)
qc.append(ctrl_add_gate, qubits)

# --- 測定 ---
qc.measure(x_reg, c_res[0:n_qubits])
qc.measure(y_reg, c_res[n_qubits:2*n_qubits])
qc.measure(z_reg, c_res[2*n_qubits:3*n_qubits])

# --- 実行 ---
print("Simulating (MPS) for k=3 (P + 2P)...")
backend = AerSimulator(method='matrix_product_state')
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

# 座標変換
if Z_proj == 0:
    affine_res = (None, None)
else:
    inv_z = pow(Z_proj, -1, N)
    x_affine = (X_proj * (inv_z**2)) % N
    y_affine = (Y_proj * (inv_z**3)) % N
    affine_res = (x_affine, y_affine)

# 期待値 k=3 => 3P
expected = scalar_mult._classical_scalar_mult(3, P, N)
print(f"Expected (3P): {expected}")
print(f"Measured:      {affine_res}")

if affine_res == expected:
    print("\nSUCCESS! Chain addition works.")
else:
    print("\nMISMATCH")