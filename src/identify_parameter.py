def find_curve_parameters(p, x, y):
    print(f"Searching parameters for point ({x}, {y}) on y^2 = x^3 + ax + b mod {p}")
    valid_params = []
    rhs_target = (y**2) % p
    x3 = (x**3) % p
    
    for a in range(p):
        for b in range(p):
            # 方程式を満たすか確認
            if (x3 + a*x + b) % p == rhs_target:
                # 特異点チェック (判別式 Delta = -16(4a^3 + 27b^2) != 0)
                delta = (4 * a**3 + 27 * b**2) % p
                if delta != 0:
                    valid_params.append((a, b))
                    
    return valid_params

# 実行
params = find_curve_parameters(13, 11, 5)
print(f"Valid (a, b) pairs: {params}")