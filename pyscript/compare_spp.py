import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ============================================================
#  常量
# ============================================================
R_WGS84 = 6378137.0
F_WGS84 = 1.0 / 298.257223563
PAI     = 3.1415926535898


def xyz_to_blh(x, y, z):
    e2 = 2 * F_WGS84 - F_WGS84 * F_WGS84
    p = np.sqrt(x * x + y * y)
    if x * x + y * y + z * z < 1e-10:
        return 0.0, 0.0, -6371.0 / 2.0
    blh_b = np.arctan2(z, (x * x + y * y) * (1 - e2))
    blh_l = np.arctan2(y, x)
    delta_z = e2 * z
    for _ in range(10):
        N = R_WGS84 / np.sqrt(1 - e2 * np.sin(blh_b) ** 2)
        blh_b = np.arctan2(z + delta_z, p)
        blh_h = np.sqrt(p * p + (z + delta_z) ** 2) - N
        delta_z1 = N * e2 * np.sin(blh_b)
        if abs(delta_z - delta_z1) < 1e-10:
            break
        delta_z = delta_z1
    N = R_WGS84 / np.sqrt(1 - e2 * np.sin(blh_b) ** 2)
    blh_b = np.arctan2(z + delta_z, p)
    blh_h = np.sqrt(x * x + y * y + (z + delta_z) ** 2) - N
    return blh_b, blh_l, blh_h


def blh_to_neu_matrix(b, l):
    return np.array([
        [-np.sin(l),  np.cos(l),  0],
        [-np.sin(b)*np.cos(l), -np.sin(b)*np.sin(l),  np.cos(b)],
        [ np.cos(b)*np.cos(l),  np.cos(b)*np.sin(l),  np.sin(b)],
    ])


def comp_enu(x0, y0, z0, xs, ys, zs):
    b, l, _ = xyz_to_blh(x0, y0, z0)
    R = blh_to_neu_matrix(b, l)
    dxyz = np.array([xs - x0, ys - y0, zs - z0])
    return R @ dxyz  # [dE, dN, dU]


def load_result(filename):
    """每行: SecOfWeek  X  Y  Z  GPSSatNum  BDSSatNum  PDOP"""
    data = np.loadtxt(filename)
    return {
        't':     data[:, 0],
        'x':     data[:, 1],
        'y':     data[:, 2],
        'z':     data[:, 3],
        'gps_n': data[:, 4].astype(int),
        'bds_n': data[:, 5].astype(int),
        'pdop':  data[:, 6],
        'name':  os.path.splitext(os.path.basename(filename))[0],
    }


def compute_stats(d):
    """计算 ENU 误差和统计量，参考坐标取平均值"""
    x0, y0, z0 = d['x'].mean(), d['y'].mean(), d['z'].mean()
    n = len(d['t'])
    e, n_, u = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        de, dn, du = comp_enu(x0, y0, z0, d['x'][i], d['y'][i], d['z'][i])
        e[i], n_[i], u[i] = de, dn, du

    stats = {
        'E_mean': e.mean(), 'E_std': e.std(ddof=1), 'E_rms': np.sqrt((e**2).mean()),
        'N_mean': n_.mean(), 'N_std': n_.std(ddof=1), 'N_rms': np.sqrt((n_**2).mean()),
        'U_mean': u.mean(), 'U_std': u.std(ddof=1), 'U_rms': np.sqrt((u**2).mean()),
        'H_rms':  np.sqrt((e**2 + n_**2).mean()),
        '3D_rms': np.sqrt((e**2 + n_**2 + u**2).mean()),
        'avg_sat': d['gps_n'].mean() + d['bds_n'].mean(),
    }
    return e, n_, u, stats


# ============================================================
#  读取数据
# ============================================================
DATA_DIR = "C:/Users/asus/Documents/Program/Program_GNSS/RTK_CSH/RTK_CSH"

# 试验1：系统对比
files_sys = [
    "SPP_Result_GPSonly.txt",
    "SPP_Result_BDSonly.txt",
    "SPP_Result.txt",
]

# 试验2：高度截止角
files_elv = [
    "SPP_Result_Elv10.txt",
    "SPP_Result_Elv20.txt",
    "SPP_Result_Elv30.txt",
]

# 试验3：单频vs双频
files_freq = [
    "SPP_Result_SingleFreq.txt",
    "SPP_Result_DualFreq.txt",
]

# 试验4：对流层
files_tropo = [
    "SPP_Result_NoTropo.txt",
    "SPP_Result_WithTropo.txt",
]


def process_group(files, labels, title, png_name):
    results = []
    stats_list = []
    for f in files:
        try:
            d = load_result(os.path.join(DATA_DIR, f))
            e, n, u, s = compute_stats(d)
            results.append((d, e, n, u, s))
            stats_list.append(s)
            print(f"  {f}: HRMS={s['H_rms']:.3f}  VRMS={s['U_rms']:.3f}  "
                  f"avg卫星={s['avg_sat']:.1f}")
        except FileNotFoundError:
            print(f"  {f}: 文件未找到，跳过")

    if len(results) < 2:
        print(f"  跳过（文件不足）")
        return

    # 画图
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=14)
    colors = ['#2196F3', '#FF5722', '#4CAF50', '#9C27B0']
    markers = ['o', 's', '^', 'D']

    # 左上：水平散点
    ax = axes[0, 0]
    for idx, (d, e, n, u, s) in enumerate(results):
        ax.scatter(e, n, s=1, c=colors[idx], alpha=0.4, label=f"{labels[idx]} RMS={s['H_rms']:.2f}")
    ax.axhline(0, color='gray', lw=0.5)
    ax.axvline(0, color='gray', lw=0.5)
    ax.set_xlabel('E 误差 (m)'); ax.set_ylabel('N 误差 (m)')
    ax.set_aspect('equal'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # 右上：U 误差时间序列
    ax = axes[0, 1]
    for idx, (d, e, n, u, s) in enumerate(results):
        ax.scatter(d['t'], u, s=1, c=colors[idx], alpha=0.4, label=f"{labels[idx]}")
    ax.axhline(0, color='gray', lw=0.5)
    ax.set_xlabel('时间 (s)'); ax.set_ylabel('U 误差 (m)')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # 左下：RMS 柱状图
    ax = axes[1, 0]
    x = np.arange(len(results))
    width = 0.25
    for i, metric in enumerate(['E_rms', 'N_rms', 'U_rms']):
        vals = [s[metric] for s in stats_list]
        bars = ax.bar(x + i * width, vals, width, label={'E_rms':'E', 'N_rms':'N', 'U_rms':'U'}[metric])
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                    f'{v:.2f}', ha='center', va='bottom', fontsize=7)
    ax.set_xticks(x + width)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('RMS (m)')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3, axis='y')

    # 右下：平均卫星数 + PDOP
    ax = axes[1, 1]
    for idx, (d, e, n, u, s) in enumerate(results):
        ax.scatter(s['avg_sat'], d['pdop'].mean(), s=80,
                   c=colors[idx], marker=markers[idx], label=f"{labels[idx]}")
    ax.set_xlabel('平均卫星数'); ax.set_ylabel('平均 PDOP')
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(png_name, dpi=200)
    print(f"  图表已保存: {png_name}\n")


# ============================================================
#  执行
# ============================================================
if __name__ == "__main__":
    print("=" * 60)

    print("\n[试验1] 系统对比: GPS only vs BDS only vs GPS+BDS")
    process_group(files_sys,
                  ["GPS only", "BDS only", "GPS+BDS"],
                  "SPP 系统对比试验",
                  "compare_system.png")

    print("\n[试验2] 高度截止角: 10° vs 20° vs 30°")
    process_group(files_elv,
                  ["10°", "20°", "30°"],
                  "SPP 高度截止角试验",
                  "compare_elevation.png")

    print("\n[试验3] 单频 vs 双频 IF 组合")
    process_group(files_freq,
                  ["单频 B1/L1", "双频 IF"],
                  "SPP 单频 vs 双频试验",
                  "compare_freq.png")

    print("\n[试验4] 对流层改正: 加 vs 不加")
    process_group(files_tropo,
                  ["无对流层", "有对流层"],
                  "SPP 对流层改正试验",
                  "compare_tropo.png")

    print("=" * 60)
    print("所有图表生成完毕。")
