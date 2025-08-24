import mpmath

# --- 全局精度设置 ---
# 设置 mpmath 的二进制位精度
mpmath.mp.prec = 1000


def StatDist(P, Q):
    """
    计算两个概率分布 P 和 Q 之间的统计距离，并返回其以 2 为底的对数。
    """
    if len(P) != len(Q):
        raise ValueError("输入向量 P 和 Q 的长度必须相同。")

    n = len(P)
    res = mpmath.mpf(0)
    for i in range(n):
        tmp = mpmath.fabs(P[i] - Q[i])
        res += tmp
    res /= 2

    # --- 修正部分 ---
    if res <= 0:
        return mpmath.mpf('-inf')

    # 使用 mpmath.log(x, 2) 计算以 2 为底的对数
    return mpmath.log(res, 2)


def RenyiDiv(P, Q, a):
    """
    计算两个概率分布 P 和 Q 之间 alpha 阶的 Rényi 散度。
    """
    if len(P) != len(Q):
        raise ValueError("输入向量 P 和 Q 的长度必须相同。")
    if a == 1:
        return mpmath.mpf('nan')

    n = len(P)
    res_sum = mpmath.mpf(0)
    a_mpf = mpmath.mpf(a)

    for i in range(n):
        if Q[i] == 0 and P[i] != 0:
            return mpmath.mpf('inf')
        if P[i] == 0:
            continue

        ratio = P[i] / Q[i]
        tmp = mpmath.power(ratio, a_mpf - 1) * P[i]
        res_sum += tmp

    exponent = 1 / (a_mpf - 1)
    res_multiplier = mpmath.power(res_sum, exponent)

    if res_multiplier <= 0:
        return mpmath.mpf('-inf')

    # --- 修正部分 ---
    # 使用 mpmath.log(x, 2)
    return mpmath.log(res_multiplier-1, 2)


def RenyiDivHalf(P, Q, a):
    """
    计算两个"半"分布 P 和 Q 之间的 Rényi 散度。
    """
    if len(P) != len(Q):
        raise ValueError("输入向量 P 和 Q 的长度必须相同。")
    if a == 1:
        return mpmath.mpf('nan')

    n = len(P)
    a_mpf = mpmath.mpf(a)

    res_sum = mpmath.mpf(0)

    if Q[0] == 0 and P[0] != 0: return mpmath.mpf('inf')
    if P[0] != 0:
        res_sum += mpmath.power(P[0] / Q[0], a_mpf - 1) * P[0]

    for i in range(1, n):
        if Q[i] == 0 and P[i] != 0: return mpmath.mpf('inf')
        if P[i] == 0: continue

        ratio = P[i] / Q[i]
        tmp = mpmath.power(ratio, a_mpf - 1) * P[i]
        res_sum += (2 * tmp)

    exponent = 1 / (a_mpf - 1)
    res_multiplier = mpmath.power(res_sum, exponent)

    if res_multiplier <= 0:
        return mpmath.mpf('-inf')

    # --- 修正部分 ---
    # 使用 mpmath.log(x, 2)
    return mpmath.log(res_multiplier-1, 2)


# # --- 使用示例 ---
# if __name__ == '__main__':
#     # 设置 mpmath 打印的小数位数
#     mpmath.mp.dps = 50
#
#     # 创建两个示例概率分布 P 和 Q
#     P_dist = [mpmath.mpf(p) for p in [0.1, 0.2, 0.7]]
#     Q_dist = [mpmath.mpf(q) for q in [0.7, 0.2, 0.1]]
#
#     alpha = 2
#
#     print("-" * 40)
#     print(f"输入分布 P: {[float(p) for p in P_dist]}")
#     print(f"输入分布 Q: {[float(q) for q in Q_dist]}")
#     print(f"Rényi 阶数 alpha: {alpha}")
#     print("-" * 40)
#
#     # 1. 计算统计距离的对数
#     log_stat_dist = StatDist(P_dist, Q_dist)
#     print(f"log2(统计距离): {log_stat_dist}")
#
#     # 2. 计算 Rényi 散度
#     renyi_divergence = RenyiDiv(P_dist, Q_dist, alpha)
#     print(f"Rényi 散度 (比特): {renyi_divergence}")
#     print("-" * 40)
#
#     # 3. 验证当 P=Q 时的情况
#     print("当 P 和 Q 相同时:")
#     log_stat_dist_same = StatDist(P_dist, P_dist)
#     renyi_div_same = RenyiDiv(P_dist, P_dist, alpha)
#     print(f"log2(统计距离): {log_stat_dist_same}")
#     print(f"Rényi 散度 (比特): {renyi_div_same}")