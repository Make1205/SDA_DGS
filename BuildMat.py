import numpy as np
import mpmath
import SDSVP
from SDSVP import SDSVP


def build_mat_3(v, eps):

    # 设置 mpmath 的浮点数精度 (dps: decimal places)
    # C++ NTL 的 SetPrecision(1000) 设置的是二进制位精度，
    # 350 个十进制位的精度足以覆盖它。
    mpmath.mp.dps = 350

    n = len(v)

    # 初始化一个 (n+1)x(n+1) 的高精度浮点数矩阵
    mat_mp = np.zeros((n + 1, n + 1), dtype=object)

    # 填充矩阵的初始值
    for i in range(n):
        mat_mp[i, i] = mpmath.mpf(1)
        mat_mp[i, n] = -v[i]

    # 计算缩放因子 tmp = eps^(n+1)
    tmp = mpmath.power(eps, n + 1)

    # 应用缩放因子
    for i in range(n):
        mat_mp[i, n] /= tmp
        mat_mp[i, i] /= tmp

    # 创建最终的整数矩阵
    mat_zz = np.zeros((n + 1, n + 1), dtype=object)

    # 将高精度浮点数截断为整数
    for i in range(n):
        # Python 的 int() 函数执行截断操作 (truncate)，与 NTL 的 to_ZZ() 行为一致
        mat_zz[i, i] = int(mat_mp[i, i])
        mat_zz[i, n] = int(mat_mp[i, n])

    # 设置右下角的元素为 1
    mat_zz[n, n] = 1

    return mat_zz






