import numpy as np
import sys


def SDSVP(H, radius=None):
    """
    将 MATLAB 函数 SDSVP 精确翻译为 Python 版本。

    该函数使用球形解码算法来解决最短向量问题（SVP）。

    参数:
        H (array_like): 格基矩阵，矩阵的列构成了格的基向量。
        radius (float, optional): 搜索半径。注意：为了忠实于原始 MATLAB 代码，
                                  算法内部会直接将距离的平方与此半径值进行比较。
                                  如果为 None，则默认为浮点数的最大值，
                                  以进行无约束的搜索。

    返回:
        numpy.ndarray: 在给定半径内找到的最短格向量的整数系数向量 'x'。
                       如果没有找到非零向量，则返回一个零向量。
    """

    # ------------------ 辅助状态和递归函数 ------------------
    # 使用一个简单的类来封装状态，避免使用全局变量
    class State:
        def __init__(self):
            self.SPHDEC_RADIUS = 0.0
            self.RETVAL = None
            self.x = None
            self.NUMX = 0

    def sphdec_core(state, z, R, layer, dist):
        """
        球形解码的核心递归函数。
        """
        n = R.shape[1]

        # 1. 计算 zi 以找到搜索中心 c
        # MATLAB: zi = z - R(:,layer+1:end)*x(layer+1:end);
        if layer == n - 1:
            zi = z.copy()
        else:
            zi = z - R[:, layer + 1:] @ state.x[layer + 1:]

        # 2. 找到当前层的中心整数 c
        # MATLAB: c = round(zi(layer)/R(layer,layer));
        center = np.round(zi[layer] / R[layer, layer])
        state.x[layer] = center

        # 3. 计算以 center 为中心的点的平方距离
        # 距离计算必须基于原始目标向量 z (原点)
        # MATLAB: d = abs(z(layer) - R(layer,layer:end)*x(layer:end))^2 + dist;
        term = R[layer, layer:] @ state.x[layer:]
        d = term ** 2 + dist

        if d <= state.SPHDEC_RADIUS:
            if layer == 0:  # 到达最底层
                if d > 1e-15:  # 确保找到的不是零向量
                    state.RETVAL = state.x.copy()
                    state.SPHDEC_RADIUS = d
                    state.NUMX += 1
            else:  # 继续向下一层递归
                sphdec_core(state, z, R, layer - 1, d)

        # 4. Schnorr-Euchner 枚举策略
        delta = 0
        d_for_while_loop = d
        while d_for_while_loop <= state.SPHDEC_RADIUS:
            delta += 1
            # 检查两个新点: center - delta 和 center + delta
            for k in range(1, 3):  # k=1, 2
                ci = center + delta * ((-1) ** k)
                state.x[layer] = ci

                # 距离计算同样基于 z (原点)
                term_new = R[layer, layer:] @ state.x[layer:]
                d_new = term_new ** 2 + dist

                if d_new <= state.SPHDEC_RADIUS:
                    if layer == 0:
                        if d_new > 1e-15:
                            state.RETVAL = state.x.copy()
                            state.SPHDEC_RADIUS = d_new
                            state.NUMX += 1
                    else:
                        sphdec_core(state, z, R, layer - 1, d_new)

                # while 循环的条件依赖于 k=2 时计算出的 d_new
                if k == 2:
                    d_for_while_loop = d_new

    # ------------------ SDSVP 主函数逻辑 ------------------

    # 1. 初始化和参数检查
    if radius is None:
        radius = sys.float_info.max

    H_matrix = np.asarray(H, dtype=float)
    if H_matrix.ndim != 2:
        raise ValueError("输入 H 必须是一个二维矩阵。")

    # 如果行数小于列数，则用零填充以使其成为方阵
    m, n = H_matrix.shape
    if m < n:
        padding = np.zeros((n - m, n))
        H_matrix = np.vstack([H_matrix, padding])
        m = n  # 更新行数

    # 2. QR分解
    # --- 关键修正 ---
    # 移除 mode='economic' 参数以确保与旧版 NumPy 的最大兼容性。
    # 对于 m >= n 的矩阵，默认行为与 'economic' 相同。
    qr_result = np.linalg.qr(H_matrix)

    # 仍然保留检查，以防万一
    if not (isinstance(qr_result, tuple) and len(qr_result) == 2):
        raise TypeError(f"np.linalg.qr 返回了预料之外的类型。请升级您的 NumPy 库。")

    _, R = qr_result

    # 3. 初始化状态和目标向量
    # MATLAB: z=zeros(min(m,n),1);
    z = np.zeros(n)

    state = State()
    state.SPHDEC_RADIUS = radius
    state.RETVAL = np.zeros(n)
    state.x = np.zeros(n)
    state.NUMX = 0

    # 4. 调用核心递归函数 (初始层级为 n-1，因为Python是0索引)
    sphdec_core(state, z, R, n - 1, 0)

    # 5. 返回结果
    if state.NUMX > 0:
        return state.RETVAL
    else:
        return np.zeros(n)


# 示例用法
# if __name__ == '__main__':
#     # 示例：创建一个格基矩阵 H
#     H = np.array([
#         [999999,-367880],
#         [0, 1]
#     ])
#
#     # 调用函数寻找最短向量的系数（不设半径限制）
#     x_coeffs = SDSVP(H)
#
#     print("找到的整数系数向量 x:")
#     print(x_coeffs)
#
#     # 计算对应的格向量 v = Hx
#     if np.any(x_coeffs):
#         v = H @ x_coeffs
#         print("\n对应的格向量 v = Hx:")
#         print(v)
#         print("\n向量 v 的平方长度 ||v||^2:")
#         print(np.dot(v, v))