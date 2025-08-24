from contextlib import suppress
import numpy as np
import mpmath
from decorator import append
from sympy.codegen.cnodes import sizeof

import BuildMat
import SDSVP
import DistTest
from DistTest import StatDist,RenyiDivHalf,RenyiDiv

# --- 1. 初始化参数 ---

# 设置 mpmath 的浮点数工作精度（dps = decimal places）。
# 这对于模拟 NTL::RR 的高精度至关重要。
mpmath.mp.dps = 500

# 从 C++ 代码中转换过来的变量
t = 2
s = [13, 11, 7]
sigma = [2.8, 2.3, 1.4]
a = [200, 500, 1000]  # 注意：变量 a 在 C++ 片段中已定义，但未使用
k = 15

# 根据 t 的值选择参数
selected_s = s[t]
selected_sigma = sigma[t]

# C++: int full=s[t];
# 在 Python 中，我们将 s[t] 的值赋给 full
full = 1000

# C++: ZZ tk; power2(tk,k);
# Python 的 int 类型原生支持大整数，可以直接进行幂运算
tk = 2 ** k

print("参数初始化完成:")
print(f"  t = {t}")
print(f"  k = {k}")
print(f"  s[{t}] = {selected_s}")
print(f"  sigma[{t}] = {selected_sigma}")
print(f"  full = {full}")
print(f"  2^k = {tk}")
print(f"  mpmath 精度 (dps): {mpmath.mp.dps}")
print("-" * 30)

# --- 2. 计算概率分布 (probfull) ---

# C++: Vec<RR> probfull; RR sum = RR(0.5); probfull.append(RR(0.5));
# 初始化列表和总和变量。
# 使用 mpmath.mpf() 来创建高精度浮点数，确保精度。
# 完全遵循 C++ 代码的逻辑，将第一个元素和 sum 的初始值都设为 0.5。
probfull = []
current_sum = mpmath.mpf('0')

print(f"正在计算 'probfull' (循环 {full} 次)...")
# C++: for(int i=1;i<=full;i++)

p = mpmath.npdf(0 / selected_sigma)
probfull.append(p)
# C++: sum+=probfull[i];
current_sum += p

for i in range(1, full + 1):
    p = mpmath.npdf(i / selected_sigma) + mpmath.npdf(-i / selected_sigma)
    probfull.append(p)
    # C++: sum+=probfull[i];
    current_sum += p

print("计算完成，正在进行归一化...")
# C++: for(int i=0;i<=full;i++) { probfull[i] /= sum; }
# 归一化列表中的每一个元素
for i in range(len(probfull)):
    probfull[i] /= current_sum



probin = probfull[0:s[t]]

# for i in probin:
#     print(i * tk)
# 或者使用切片: probin = probfull[:]
# print(probin[:])
np.set_printoptions(linewidth=200, threshold=np.inf, suppress=True)
res=SDSVP.SDSVP(BuildMat.build_mat_3(probin,pow(2,-k/(s[t]))))
print(res)

# print(sizeof(probin))

probin = probfull[0:s[t]]
# print(probin)

# probFP=[1035./5966,1883./5966,1418./5966,884./5966,456./5966,195./5966,69./5966,20./5966,5./5966,1./5966,0.]
# probFP=[ 298./1024.,  462./1024. , 215./1024. ,  60./1024. ,  10./1024.  ,  1./1024. ,   0./1024.]
# print(probFP)
# print(probin)


# for i in probin:
    # probFP.append(round(i*tk)/tk)

# print(probFP)


print(StatDist(probin[0:s[t]],res[0:s[t]]/res[s[t]]))
print(RenyiDivHalf(res[0:s[t]]/res[s[t]],probin[0:s[t]],a[t]))

# print(StatDist(probfull[0:s[t]],probFP[0:s[t]]))
# print(RenyiDivHalf(probFP[0:s[t]],probfull[0:s[t]],a[t]))

