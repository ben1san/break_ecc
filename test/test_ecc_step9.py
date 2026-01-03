import time
from qiskit_aer import AerSimulator
from general.arithmetic import ModularArithmetic
from general.ecc import QuantumECC, ScalarMultiplication
from general.shor_ecdlp import ShorECDLP
from qiskit import transpile  # 追加

# --- パラメータ設定 (小規模検証用) ---
# N=5 (3ビット) 程度でロジック確認を推奨
N = 5  # 素数 p
n_qubits = 3 # 演算ビット幅
p_point = (2, 1) # 曲線上の点 P (例)
q_point = (2, 1) # 曲線上の点 Q (本来は Q = kP だが検証のため Pと同じにする等の調整可)

# --- インスタンス生成 ---
arith = ModularArithmetic(N, n_qubits) #
ecc = QuantumECC(arith) #
# scalar_mult は初期化時にダミーの a=0 を渡しておく(内部メソッド利用のみのため)
sm = ScalarMultiplication(ecc, arith, a=2) # aは曲線パラメータ y^2 = x^3 + ax + b

# ShorECDLP インスタンス
shor_solver = ShorECDLP(sm, N)

# --- 回路構築 ---
ctrl_bits = 2 # 制御レジスタのビット数 (小さめに設定)
print(f"Constructing circuit with control bits={ctrl_bits}...")
start_time = time.time()

qc = shor_solver.construct_circuit(
    num_ctrl_qubits=ctrl_bits,
    point_P=p_point,
    point_Q=q_point,
    initial_point_Z1=1
)

# 【ここが修正ポイント】
# シミュレータが実行可能な形式（基本ゲート）に回路を変換（トランスパイル）します
print("Transpiling circuit...")
simulator = AerSimulator(method='matrix_product_state')
transpiled_qc = transpile(qc, simulator) 

print(f"Transpiled Depth: {transpiled_qc.depth()}") # 実際の深さを確認

# --- シミュレーション ---
print("Starting simulation...")
# トランスパイル済みの回路 (transpiled_qc) を渡します
job = simulator.run(transpiled_qc, shots=16) 
result = job.result()
counts = result.get_counts()

print("Measurement Results (reg_b reg_a):")
print(counts)