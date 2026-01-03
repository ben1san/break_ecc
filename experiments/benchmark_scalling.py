import sys
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator

# srcディレクトリをパスに追加（experimentsから見て一つ上）
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))

# 必要なクラスをインポート
# ※ src内のファイル構成に合わせてimportを変更しています
from general.arithmetic import ModularArithmetic
from general.ecc import QuantumECC, ScalarMultiplication
from general.shor_ecdlp import ShorECDLP 

def get_resource_estimates(n_bits):
    """
    指定されたビット数のダミーECDLP回路を構築し、リソース量を見積もる
    注: 実際に解けるパラメータである必要はなく、回路規模が測定できれば良い
    """
    # 1. パラメータ設定 (ビット数に応じた適当な素数pと曲線パラメータ)
    # nビットで表現できる最大の素数に近いものを使用（回路規模を最大化するため）
    # ここでは簡易的に 2^n - 1 に近い奇数を設定（厳密な素数でなくても回路生成は可能）
    p = 2**n_bits - 1 
    if p % 2 == 0: p -= 1
    
    # ダミーの点P, Q (座標の数値自体は回路構造に大きな影響を与えないため適当でOK)
    # ただし、ModularArithmeticが正常に動く範囲の整数であること
    point_P = (2, 1)
    point_Q = (3, 1)
    
    # 2. クラスの初期化
    # 制御ビット数もビット数に応じてスケーリングさせる (例: n_bitsと同じにする)
    num_ctrl_qubits = n_bits 
    
    arith = ModularArithmetic(p, n_bits)
    ecc = QuantumECC(arith)
    sm = ScalarMultiplication(ecc, arith, a=0)
    shor = ShorECDLP(sm, p)
    
    # 3. 回路構築 (Construct)
    start_build = time.time()
    # construct_circuitの引数仕様に合わせて調整してください
    qc = shor.construct_circuit(num_ctrl_qubits, point_P, point_Q)
    build_time = time.time() - start_build
    
    # 4. リソース計測 (論理ゲート数)
    # トランスパイルなしの純粋な構成ゲート数
    ops = qc.count_ops()
    depth_logical = qc.depth()
    qubits = qc.num_qubits
    
    # 5. トランスパイル (物理ゲート数に近い見積もり)
    # ※ ビット数が大きいと時間がかかるため、Optimization level=0 (最速) で計測
    # ※ 実際のデバイス向けではないため basis_gates は指定しない、または汎用的な ['cx', 'u'] などにする
    print(f"  Transpiling {n_bits}-bit circuit...")
    start_transpile = time.time()
    qc_transpiled = transpile(qc, optimization_level=0, basis_gates=['cx', 'u', 'p'])
    transpile_time = time.time() - start_transpile
    
    ops_transpiled = qc_transpiled.count_ops()
    cx_count = ops_transpiled.get('cx', 0)
    depth_transpiled = qc_transpiled.depth()
    
    return {
        "n_bits": n_bits,
        "p_approx": p,
        "qubits": qubits,
        "logical_depth": depth_logical,
        "logical_gate_count": sum(ops.values()),
        "transpiled_cx_count": cx_count,
        "transpiled_depth": depth_transpiled,
        "build_time": build_time,
        "transpile_time": transpile_time
    }

def run_benchmark(max_bits=5):
    """
    ビット数を変えながらベンチマークを実行
    注意: 汎用回路は指数関数的に深くなるため、max_bitsは慎重に設定する
    """
    results = []
    print(f"Starting Phase 1 Benchmark (Max bits: {max_bits})...")
    
    for n in range(2, max_bits + 1):
        print(f"Benchmarking n={n}...")
        try:
            res = get_resource_estimates(n)
            results.append(res)
            print(f"  -> Qubits: {res['qubits']}, CX: {res['transpiled_cx_count']}, Depth: {res['transpiled_depth']}")
        except Exception as e:
            print(f"  -> Error at n={n}: {e}")
            break
            
    return pd.DataFrame(results)

def plot_results(df):
    """結果をプロットして保存"""
    if df.empty:
        print("No data to plot.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: CX Count vs Bit Size
    axes[0].plot(df['n_bits'], df['transpiled_cx_count'], 'o-', label='CX Count (Phase 1)')
    axes[0].set_title('Scaling of CNOT Gates')
    axes[0].set_xlabel('Bit Size (n)')
    axes[0].set_ylabel('Number of CX Gates')
    axes[0].set_yscale('log') # 対数グラフにする
    axes[0].grid(True, which="both", ls="--")
    axes[0].legend()
    
    # Plot 2: Depth vs Bit Size
    axes[1].plot(df['n_bits'], df['transpiled_depth'], 's-', color='orange', label='Circuit Depth (Phase 1)')
    axes[1].set_title('Scaling of Circuit Depth')
    axes[1].set_xlabel('Bit Size (n)')
    axes[1].set_ylabel('Depth')
    axes[1].set_yscale('log') # 対数グラフにする
    axes[1].grid(True, which="both", ls="--")
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig('phase1_scaling_benchmark.png')
    print("Plot saved to experiments/phase1_scaling_benchmark.png")

if __name__ == "__main__":
    # n=6くらいでかなり重くなる可能性があります。まずは小さく試してください。
    df = run_benchmark(max_bits=6) 
    df.to_csv('phase1_benchmark_results.csv', index=False)
    print("Results saved to experiments/phase1_benchmark_results.csv")
    plot_results(df)