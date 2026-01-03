import operator

class ShorPostProcessor:
    def __init__(self, curve_order, p_mod, a, b):
        self.r = curve_order
        self.p = p_mod
        self.a = a
        self.b = b

    def _classical_point_mult(self, k, point):
        """検証用の古典スカラー倍算 (k * P)"""
        if k == 0: return (None, None)
        if k == 1: return point
        
        # Double-and-Add
        result = (None, None)
        add_point = point
        
        for i in range(k.bit_length()):
            if (k >> i) & 1:
                result = self._point_add(result, add_point)
            add_point = self._point_double(add_point)
        return result

    def _point_add(self, p1, p2):
        if p1 == (None, None): return p2
        if p2 == (None, None): return p1
        x1, y1 = p1
        x2, y2 = p2
        
        if x1 == x2 and y1 != y2: return (None, None)
        if x1 == x2 and y1 == y2: return self._point_double(p1)
        
        num = (y2 - y1) % self.p
        den = (x2 - x1) % self.p
        lam = (num * pow(den, -1, self.p)) % self.p
        x3 = (lam**2 - x1 - x2) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return (x3, y3)

    def _point_double(self, p1):
        if p1 == (None, None): return (None, None)
        x1, y1 = p1
        if y1 == 0: return (None, None)
        
        num = (3 * x1**2 + self.a) % self.p
        den = (2 * y1) % self.p
        lam = (num * pow(den, -1, self.p)) % self.p
        x3 = (lam**2 - 2*x1) % self.p
        y3 = (lam * (x1 - x3) - y1) % self.p
        return (x3, y3)

    def solve_k(self, counts, point_P, point_Q):
        """
        測定結果(counts)からkを推定する
        """
        # 1. カウントの多い順にソート
        sorted_counts = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
        
        candidates = {}
        
        print(f"\n--- Post-Processing Results (Order r={self.r}) ---")
        print(f"Target: Find k such that Q = k*P")
        print(f"P={point_P}, Q={point_Q}")

        for bitstring, count in sorted_counts:
            # ビット列 '10 01' を分割 (reg_b, reg_a)
            # 注意: Qiskitの出力順序や空白に依存するため柔軟にパース
            parts = bitstring.split()
            if len(parts) == 2:
                s_b, s_a = parts
            else:
                # 空白がない場合、半分で分割（ビット数が同じと仮定）
                mid = len(bitstring) // 2
                s_b, s_a = bitstring[:mid], bitstring[mid:]
            
            val_a = int(s_a, 2) # u
            val_b = int(s_b, 2) # v
            
            # 関係式: v ≡ k * u (mod r) を解く
            # k = v * u^(-1) (mod r)
            
            # uがrの倍数(0など)だと逆元がないのでスキップ
            if val_a % self.r == 0:
                continue
                
            try:
                inv_a = pow(val_a, -1, self.r)
                k_candidate = (val_b * inv_a) % self.r
                
                # 候補の出現頻度を集計
                candidates[k_candidate] = candidates.get(k_candidate, 0) + count
                
                # 検証
                calculated_Q = self._classical_point_mult(k_candidate, point_P)
                is_correct = (calculated_Q == point_Q)
                
                mark = "✅ CORRECT" if is_correct else "❌"
                print(f"Meas(b={val_b}, a={val_a}) -> counts={count} -> Candidate k={k_candidate} : {mark}")
                
                if is_correct and candidates[k_candidate] >= 2:
                     # 複数回観測されたら早期終了しても良い
                     pass

            except ValueError:
                continue

        print("\n--- Summary ---")
        if not candidates:
            print("No valid candidates found.")
            return None

        best_k = max(candidates.items(), key=operator.itemgetter(1))[0]
        print(f"Most likely k: {best_k} (supported by {candidates[best_k]} shots)")
        
        # 最終確認
        final_Q = self._classical_point_mult(best_k, point_P)
        if final_Q == point_Q:
            print(f"SUCCESS: Discrete Logarithm Found! k = {best_k}")
            return best_k
        else:
            print("FAILED: Verification failed.")
            return None

# --- 実行 ---
# 設定値 (N=5, P=(2,1) の場合)
N = 5
curve_a = 2
curve_b = 4 # y^2 = x^3 + 2x + 4
order_r = 7 # 位数

# 先ほどのシミュレーション結果
# (実際のコードでは result.get_counts() をそのまま渡してください)
counts = {'11 11': 1, '00 11': 2, '01 11': 1, '10 11': 1, '10 10': 2, '10 01': 3, '11 10': 1, '11 00': 1, '01 00': 1, '11 01': 1, '01 01': 2}
p_point = (2, 1)
q_point = (2, 1)

processor = ShorPostProcessor(order_r, N, curve_a, curve_b)
k = processor.solve_k(counts, p_point, q_point)