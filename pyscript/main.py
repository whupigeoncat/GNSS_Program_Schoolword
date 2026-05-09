import numpy as np
import matplotlib.pyplot as plt


# Windows 中文显示
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ============================================================
#  配置区
# ============================================================
INPUT_FILE = "RTK_CSH/SPP_Result.txt"   # 从 C++ 程序输出的结果文件
OUTPUT_PNG  = "enu_error.png"    # 输出的散点图
# 文件格式（每行）：
#   SecOfWeek  X  Y  Z  GPSSatNum  BDSSatNum  PDOP

# ============================================================
#  常量（与 C++ 保持一致）
# ============================================================
R_WGS84 = 6378137.0
F_WGS84 = 1.0 / 298.257223563
PAI     = 3.1415926535898
Rad     = PAI / 180.0


# ============================================================
#  坐标转换（与 C++ Coord_trans.cpp 逻辑一致）
# ============================================================
def xyz_to_blh(x, y, z, A=R_WGS84, F=F_WGS84):
    """空间直角坐标 → 大地坐标（弧度）"""
    e2 = 2 * F - F * F
    p = np.sqrt(x * x + y * y)
    if x * x + y * y + z * z < 1e-10:
        return 0.0, 0.0, -6371.0 / 2.0

    blh_b = np.arctan2(z, (x * x + y * y) * (1 - e2))
    blh_l = np.arctan2(y, x)
    blh_h = 0.0

    delta_z_0 = e2 * z
    for _ in range(10):
        N = A / np.sqrt(1 - e2 * np.sin(blh_b) ** 2)
        blh_b = np.arctan2(z + delta_z_0, p)
        blh_h = np.sqrt(p * p + (z + delta_z_0) ** 2) - N
        delta_z_1 = N * e2 * np.sin(blh_b)
        if abs(delta_z_0 - delta_z_1) < 1e-10:
            break
        delta_z_0 = delta_z_1

    N = A / np.sqrt(1 - e2 * np.sin(blh_b) ** 2)
    blh_b = np.arctan2(z + delta_z_0, p)
    blh_h = np.sqrt(x * x + y * y + (z + delta_z_1) ** 2) - N
    return blh_b, blh_l, blh_h


def blh_to_neu_matrix(b, l):
    """BLH（弧度）→ ENU 转换矩阵 R"""
    R = np.array([
        [-np.sin(l),              np.cos(l),              0],
        [-np.sin(b) * np.cos(l), -np.sin(b) * np.sin(l),  np.cos(b)],
        [ np.cos(b) * np.cos(l),  np.cos(b) * np.sin(l),  np.sin(b)],
    ])
    return R


def comp_enu_error(x0, y0, z0, xs, ys, zs):
    """计算 (xs,ys,zs) 相对于参考点 (x0,y0,z0) 的 ENU 误差"""
    b, l, _ = xyz_to_blh(x0, y0, z0)
    R = blh_to_neu_matrix(b, l)
    dxyz = np.array([xs - x0, ys - y0, zs - z0])
    dE, dN, dU = R @ dxyz
    return dE, dN, dU


# ============================================================
#  主程序
# ============================================================
def main():
    # 1. 读数据
    #    格式: SecOfWeek  X  Y  Z  GPSSatNum  BDSSatNum  PDOP
    data = np.loadtxt(INPUT_FILE)
    t     = data[:, 0]      # 历元时间 (s)
    x     = data[:, 1]
    y     = data[:, 2]
    z     = data[:, 3]
    gps_n = data[:, 4].astype(int)
    bds_n = data[:, 5].astype(int)
    pdop  = data[:, 6]

    n_epochs = len(t)
    print(f"历元总数: {n_epochs}")

    # 2. 用平均值作为参考真值
    x0 = np.mean(x)
    y0 = np.mean(y)
    z0 = np.mean(z)
    print(f"参考坐标: X={x0:.3f}  Y={y0:.3f}  Z={z0:.3f}")

    # 3. 计算每个历元的 ENU 误差
    e_err = np.zeros(n_epochs)
    n_err = np.zeros(n_epochs)
    u_err = np.zeros(n_epochs)

    for i in range(n_epochs):
        e_err[i], n_err[i], u_err[i] = comp_enu_error(x0, y0, z0, x[i], y[i], z[i])

    # 4. 统计量
    print("\n========== 精度统计 (ENU) ==========")
    for name, err in [("E方向 (东)", e_err), ("N方向 (北)", n_err), ("U方向 (天)", u_err)]:
        mean = np.mean(err)
        std  = np.std(err, ddof=1)
        rms  = np.sqrt(np.mean(err ** 2))
        print(f"{name}:")
        print(f"  均值 = {mean:.3f} m")
        print(f"  标准差 = {std:.3f} m")
        print(f"  RMS  = {rms:.3f} m")

    # 水平 RMS
    h_rms = np.sqrt(np.mean(e_err ** 2 + n_err ** 2))
    print(f"\n水平 RMS = {h_rms:.3f} m")
    print(f"三维 RMS = {np.sqrt(np.mean(e_err**2 + n_err**2 + u_err**2)):.3f} m")

    # 5. 散点图
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    fig.suptitle("SPP 定位误差分析 (GPS+BDS 双频)", fontsize=14)

    # 左列：时间序列散点图
    ax = axes[0, 0]
    ax.scatter(t, e_err, s=1, c='royalblue', alpha=0.6)
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_ylabel("E 误差 (m)")
    ax.set_xlabel("SecOfWeek (s)")

    ax = axes[1, 0]
    ax.scatter(t, n_err, s=1, c='forestgreen', alpha=0.6)
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_ylabel("N 误差 (m)")
    ax.set_xlabel("SecOfWeek (s)")

    ax = axes[2, 0]
    ax.scatter(t, u_err, s=1, c='crimson', alpha=0.6)
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_ylabel("U 误差 (m)")
    ax.set_xlabel("SecOfWeek (s)")

    # 右列：ENU 散点分布
    ax = axes[0, 1]
    ax.scatter(e_err, n_err, s=2, c='purple', alpha=0.5)
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel("E 误差 (m)")
    ax.set_ylabel("N 误差 (m)")
    ax.set_aspect('equal')

    ax = axes[1, 1]
    ax.scatter(np.sqrt(e_err**2 + n_err**2), u_err, s=2, c='darkorange', alpha=0.5)
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlabel("水平误差 (m)")
    ax.set_ylabel("U 误差 (m)")

    ax = axes[2, 1]
    ax.scatter(t, pdop, s=1, c='black', alpha=0.5)
    ax.set_ylabel("PDOP")
    ax.set_xlabel("SecOfWeek (s)")

    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=200)
    print(f"\n散点图已保存: {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
