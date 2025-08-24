import numpy as np
import mpmath
import BuildMat
from SDSVP import SDSVP


# --- 1. 初始化参数 ---
# 设置 mpmath 的浮点数精度。
# NTL 的 SetPrecision(1000) 设置的是二进制位精度。
# 1000 bits 约等于 301 个十进制位，这里设置为 350 以确保足够。
mpmath.mp.prec = 100

# 定义与 C++ 代码中相同的变量
n = 1000
s = 18
a = 513
k = 72
full = 1000  # 注意：这个值很大，计算 probfull 会非常耗时
sigma = 1.8205

# 计算 2^k，Python 的 int 支持任意精度，等同于 NTL::ZZ
tk = 2 ** k

# 初始化两个列表来存储概率分布
prob = []
probfull = []

print("参数初始化完成:")
print(f"  s = {s}")
print(f"  full = {full}")
print(f"  sigma = {sigma}")
print(f"  精度 (二进制位): {mpmath.mp.prec}")
print("-" * 30)

# --- 2. 计算第一个概率分布 (prob) ---

print(f"正在计算 'prob' 列表 (循环 {s + 1} 次)...")
# 使用 mpmath.mpf 初始化高精度浮点数
current_sum = mpmath.mpf(0)

# 循环计算未归一化的概率值
for i in range(s + 1):
    # 计算标准正态分布在点 x = i/sigma 处的概率密度值
    # mpmath.npdf(x, mu, sigma) 默认 mu=0, sigma=1
    p = mpmath.npdf(i / sigma)
    prob.append(p)
    current_sum += p

# 归一化概率分布
for i in range(s + 1):
    prob[i] /= current_sum

current_sum_full = mpmath.mpf(0)
# 循环计算未归一化的概率值
for i in range(full + 1):
    p = mpmath.npdf(i / sigma)
    probfull.append(p)
    current_sum_full += p

    # 归一化概率分布
for i in range(full + 1):
    probfull[i] /= current_sum_full

# --- 准备输入数据 ---
# 1. 设置 mpmath 的工作精度
mpmath.mp.dps = 50  # 在主程序中可以根据需要设置精度
v_in = probfull[0:s]
epsilon = mpmath.mpf(pow(2, -k / s))

# print("输入向量 v:")
# for val in v_in:
#     print(f"  {val}")
# print(f"\n精度参数 eps: {epsilon}")
# print("-" * 40)

# --- 调用函数 ---
lattice_basis = BuildMat.build_mat_3(v_in, epsilon)
print(lattice_basis)
print(SDSVP(lattice_basis))
# --- 打印结果 ---
# print("构建的格基矩阵 (BuildMat3 的输出):")
# # 设置 numpy 的打印选项，以便完整显示大整数
# np.set_printoptions(linewidth=200, threshold=np.inf)
# print(lattice_basis)
#
# # 验证矩阵的数据类型
# print(f"\n矩阵元素的类型: {type(lattice_basis[0, 0])}")