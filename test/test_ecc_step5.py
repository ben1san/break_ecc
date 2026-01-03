from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_aer import AerSimulator
from general.arithmetic import ModularArithmetic
from general.ecc import QuantumECC, ScalarMultiplication

# --- パラメータ設定 ---
N = 13          
n_qubits = 8    
a = 0           
P = (11, 5)     
Q = (7, 5)      

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
# 【修正】結果測定用ビット数を増やす (X, Y, Z の3つ分)
c_res = ClassicalRegister(3 * n_qubits, 'res')

qc = QuantumCircuit(k_reg, x_reg, y_reg, z_reg, anc_reg, c_res)

# --- 初期化 ---
def set_value(circuit, reg, val):
    for i in range(len(reg)):
        if (val >> i) & 1:
            circuit.x(reg[i])

set_value(qc, x_reg, 11)
set_value(qc, y_reg, 5)
set_value(qc, z_reg, 1) # Z初期値は1

qc.x(k_reg[0]) # Enable addition

# --- ゲート適用 ---
reg_specs = {'p_size': n_qubits, 'anc_size': len(anc_reg)}
ctrl_add_gate = scalar_mult.create_controlled_add_gate(Q, reg_specs)

qubits = [k_reg[0]] + list(x_reg) + list(y_reg) + list(z_reg) + list(anc_reg)
qc.append(ctrl_add_gate, qubits)

# --- 測定 (Zも測定する) ---
qc.measure(x_reg, c_res[0 : n_qubits])
qc.measure(y_reg, c_res[n_qubits : 2*n_qubits])
qc.measure(z_reg, c_res[2*n_qubits : 3*n_qubits])

# --- 実行 ---
print(f"Test: Adding P{P} + Q{Q} (Projective Check)...")
print("Simulating (MPS)...")
backend = AerSimulator(method='matrix_product_state')
t_qc = transpile(qc, basis_gates=['u', 'cx', 'p', 'swap', 'x', 'id', 'measure'], optimization_level=1)
result = backend.run(t_qc, shots=1).result()
counts = result.get_counts()

# --- 結果確認 ---
measured_hex = list(counts.keys())[0]
# ビット列をパース (順序は measure の逆順になることに注意: Z, Y, X)
# c_resの定義順: X(low), Y(mid), Z(high) -> Qiskitの出力文字列: Z...Y...X
len_bin = 3 * n_qubits
val_z = int(measured_hex[0 : n_qubits], 2)
val_y = int(measured_hex[n_qubits : 2*n_qubits], 2)
val_x = int(measured_hex[2*n_qubits : 3*n_qubits], 2)

# Modulo N
X_proj = val_x % N
Y_proj = val_y % N
Z_proj = val_z % N

print(f"\nMeasured Projective (mod {N}):")
print(f"  X: {X_proj} (raw {val_x})")
print(f"  Y: {Y_proj} (raw {val_y})")
print(f"  Z: {Z_proj} (raw {val_z})")

# --- 座標変換: Projective -> Affine ---
# x = X / Z^2, y = Y / Z^3
if Z_proj == 0:
    print("Result is Point at Infinity (Z=0)")
    affine_res = (None, None)
else:
    inv_z = pow(Z_proj, -1, N)
    inv_z2 = (inv_z ** 2) % N
    inv_z3 = (inv_z ** 3) % N
    
    x_affine = (X_proj * inv_z2) % N
    y_affine = (Y_proj * inv_z3) % N
    affine_res = (x_affine, y_affine)

print(f"\nConverted to Affine: {affine_res}")

expected_3P = scalar_mult._classical_add(P, Q, N)
print(f"Expected Classical:  {expected_3P}")

if affine_res == expected_3P:
    print("\nSUCCESS! Quantum addition matches classical result.")
else:
    print("\nMISMATCH")