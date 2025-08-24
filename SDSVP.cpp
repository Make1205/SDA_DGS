#include <iostream>
#include <vector>
#include <limits>

// NTL库的头文件
#include <NTL/mat_RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h> // <--- 修正 #2: 添加缺失的头文件

// 使用NTL命名空间
NTL_CLIENT

// --- 修正 #1: 手动实现 Cholesky 分解 ---
// 计算对称正定矩阵 G 的 Cholesky 分解 G = R^T * R
// 其中 R 是一个上三角矩阵。
void Cholesky(mat_RR& R, const mat_RR& G) {
    long n = G.NumRows();
    if (G.NumCols() != n) {
        throw std::runtime_error("Cholesky: Matrix must be square.");
    }

    R.SetDims(n, n);
    clear(R); // 初始化为零矩阵

    for (long j = 0; j < n; ++j) {
        RR sum_sq(0);
        for (long k = 0; k < j; ++k) {
            sum_sq += sqr(R[k][j]);
        }

        RR diag_val = G[j][j] - sum_sq;
        if (diag_val <= 0) {
            // 如果对角线元素非正，矩阵不是正定的
            throw std::runtime_error("Cholesky: Matrix is not positive definite.");
        }
        R[j][j] = SqrRoot(diag_val);

        for (long i = j + 1; i < n; ++i) {
            RR sum_prod(0);
            for (long k = 0; k < j; ++k) {
                sum_prod += R[k][i] * R[k][j];
            }
            R[j][i] = (G[j][i] - sum_prod) / R[j][j];
        }
    }
    // 我们需要的是上三角矩阵 R，但上面的计算方式生成的是下三角 L 的转置
    // 所以我们需要将结果转置以得到标准的上三角 R
    R = transpose(R);
}


// 辅助类，用于执行球形译码的核心递归逻辑 (无改动)
class SphereDecoder {
public:
    SphereDecoder(const mat_RR& R, RR& radius_sq, vec_ZZ& best_x, long& num_solutions)
        : R_(R), n_(R.NumCols()), radius_sq_(radius_sq), best_x_(best_x), num_solutions_(num_solutions) {
        x_.SetLength(n_);
        dist_.SetLength(n_);
    }

    void solve() {
        recursive_search(n_ - 1);
    }

private:
    void recursive_search(long k) {
        RR zi_k(0);
        for (long j = k + 1; j < n_; ++j) {
            zi_k += R_[k][j] * to_RR(x_[j]);
        }
        zi_k = -zi_k;

        RR center_prime = zi_k / R_[k][k];
        ZZ center = RoundToZZ(center_prime);

        for (long delta = 0; ; ++delta) {
            bool branch_alive = false;
            for (int sign_val = -1; sign_val <= 1; sign_val += 2) {
                ZZ u;
                if (delta == 0 && sign_val == 1) continue;
                
                if (delta == 0) {
                    u = center;
                } else {
                    u = center + sign_val * delta;
                }

                x_[k] = u;
                RR term = (to_RR(u) - center_prime) * R_[k][k];
                RR current_dist_sq = (k + 1 < n_) ? dist_[k + 1] : RR(0);
                dist_[k] = current_dist_sq + sqr(term);

                if (dist_[k] < radius_sq_) {
                    branch_alive = true;
                    if (k == 0) {
                        if (!IsZero(x_)) {
                            radius_sq_ = dist_[0];
                            best_x_ = x_;
                            num_solutions_++;
                        }
                    } else {
                        recursive_search(k - 1);
                    }
                }
            }
            if (!branch_alive && delta > 0) {
                break;
            }
        }
    }

    const mat_RR& R_;
    long n_;
    vec_ZZ x_;
    vec_RR dist_;
    RR& radius_sq_;
    vec_ZZ& best_x_;
    long& num_solutions_;
};


// 主函数，功能等同于 MATLAB 的 SDSVP
void SDSVP(vec_ZZ& r, const mat_RR& H_in, const RR& radius) {
    mat_RR H = H_in;
    long m = H.NumRows();
    long n = H.NumCols();

    if (m < n) {
        H.SetDims(n, n);
        for (long i = m; i < n; ++i) {
            for (long j = 0; j < n; ++j) {
                H[i][j] = 0;
            }
        }
    }

    mat_RR G;
    G.SetDims(n, n);
    G = transpose(H) * H;

    mat_RR R;
    Cholesky(R, G); // <-- 使用我们新实现的 Cholesky 函数

    RR radius_sq = sqr(radius);
    vec_ZZ best_x;
    best_x.SetLength(n);
    long num_solutions = 0;

    SphereDecoder decoder(R, radius_sq, best_x, num_solutions);
    decoder.solve();

    if (num_solutions > 0) {
        r = best_x;
    } else {
        r.SetLength(n);
        clear(r);
    }
}

void SDSVP(vec_ZZ& r, const mat_RR& H) {
    RR max_rad = SqrRoot(to_RR(std::numeric_limits<double>::max()));
    SDSVP(r, H, max_rad);
}


// --- 主函数：示例和测试 ---
int main() {
    RR::SetPrecision(200);

    mat_RR H;
    H.SetDims(2, 2);
    H[0][0] = 999999; H[0][1] = -367880; 
    H[1][0] = 0; H[1][1] = 1; 
    
    
    cout << "基矩阵 H (列向量是基):" << endl;
    cout << H << endl;

    vec_ZZ shortest_vector_coords;
    
    SDSVP(shortest_vector_coords, H);

    cout << "\n找到的最短向量的坐标 r:" << endl;
    cout << shortest_vector_coords << endl;

    // --- 修正 #3: 手动将 mat_RR 转换为 mat_ZZ ---
    mat_ZZ H_zz;
    H_zz.SetDims(H.NumRows(), H.NumCols());
    for (long i = 0; i < H.NumRows(); ++i) {
        for (long j = 0; j < H.NumCols(); ++j) {
            H_zz[i][j] = RoundToZZ(H[i][j]);
        }
    }

    // 计算并输出实际的格向量 v = H * r
    // 这是基向量的整数线性组合
    vec_ZZ v = H_zz * shortest_vector_coords;

    cout << "\n对应的格向量 v = H * r:" << endl;
    cout << v << endl;

    return 0;
}