from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


ROOT = Path(__file__).resolve().parents[1]
INPUT_FILE = ROOT / "SPP_Result.txt"
OUT_DIR = ROOT

R_WGS84 = 6378137.0
F_WGS84 = 1.0 / 298.257223563
BESTPOS_XYZ = np.array([-2267335.6694, 5008649.1555, 3222374.9736])


def load_spp_result(path: Path) -> tuple[np.ndarray, int]:
    """Read SPP result rows and skip incomplete tail lines."""
    rows: list[list[float]] = []
    skipped = 0
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if len(parts) != 7:
                skipped += 1
                continue
            try:
                rows.append([float(v) for v in parts])
            except ValueError:
                skipped += 1
    if not rows:
        raise ValueError(f"No valid SPP result rows found in {path}")
    return np.asarray(rows, dtype=float), skipped


def unwrap_sec_of_week(t: np.ndarray) -> np.ndarray:
    t_cont = np.asarray(t, dtype=float).copy()
    offset = 0.0
    for i in range(1, len(t_cont)):
        if t_cont[i] + offset < t_cont[i - 1]:
            offset += 604800.0
        t_cont[i] += offset
    return t_cont


def xyz_to_blh(x: float, y: float, z: float) -> tuple[float, float, float]:
    e2 = 2 * F_WGS84 - F_WGS84 * F_WGS84
    p = np.hypot(x, y)
    if x * x + y * y + z * z < 1e-10:
        return 0.0, 0.0, -6371.0 / 2.0

    b = np.arctan2(z, p * (1 - e2))
    l = np.arctan2(y, x)
    dz = e2 * z
    h = 0.0
    for _ in range(10):
        n = R_WGS84 / np.sqrt(1 - e2 * np.sin(b) ** 2)
        b = np.arctan2(z + dz, p)
        h = np.sqrt(p * p + (z + dz) ** 2) - n
        dz_next = n * e2 * np.sin(b)
        if abs(dz - dz_next) < 1e-10:
            dz = dz_next
            break
        dz = dz_next

    n = R_WGS84 / np.sqrt(1 - e2 * np.sin(b) ** 2)
    h = np.sqrt(p * p + (z + dz) ** 2) - n
    return b, l, h


def enu_matrix(b: float, l: float) -> np.ndarray:
    return np.array(
        [
            [-np.sin(l), np.cos(l), 0.0],
            [-np.sin(b) * np.cos(l), -np.sin(b) * np.sin(l), np.cos(b)],
            [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)],
        ]
    )


def compute_enu(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ref = BESTPOS_XYZ.copy()
    b, l, _ = xyz_to_blh(*ref)
    r = enu_matrix(b, l)
    dxyz = np.column_stack([x - ref[0], y - ref[1], z - ref[2]])
    enu = dxyz @ r.T
    return ref, enu


def rms(v: np.ndarray) -> float:
    return float(np.sqrt(np.mean(v * v)))


def percentile_abs(v: np.ndarray, p: float) -> float:
    return float(np.percentile(np.abs(v), p))


def robust_ylim(v: np.ndarray, floor: float = 0.5) -> float:
    return max(float(np.percentile(np.abs(v), 99.5) * 1.25), rms(v) * 3.0, floor)


def save_enu_timeseries(t_hour: np.ndarray, enu: np.ndarray) -> Path:
    fig_path = OUT_DIR / "report_fig_enu_timeseries.png"
    colors = ["#1f5fbf", "#198754", "#c0392b"]
    labels = ["E方向误差 / m", "N方向误差 / m", "U方向误差 / m"]
    names = ["E", "N", "U"]

    fig, axes = plt.subplots(3, 1, figsize=(11.5, 7.2), sharex=True)
    fig.suptitle("GPS/BDS双频SPP定位误差时间序列", fontsize=16, fontweight="bold")
    for i, ax in enumerate(axes):
        err = enu[:, i]
        ax.scatter(t_hour, err, s=1.0, color=colors[i], alpha=0.45, edgecolors="none")
        ax.axhline(0, color="#737373", lw=0.8)
        ax.set_ylabel(labels[i], fontsize=11)
        ax.set_ylim(-robust_ylim(err), robust_ylim(err))
        ax.grid(True, color="#dddddd", lw=0.6)
        ax.text(
            0.985,
            0.88,
            f"{names[i]} RMS = {rms(err):.3f} m",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#bdbdbd", alpha=0.95),
        )
    axes[-1].set_xlabel("观测时间 / h", fontsize=11)
    axes[-1].set_xlim(t_hour[0], t_hour[-1])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path


def save_enu_distribution(enu: np.ndarray) -> Path:
    fig_path = OUT_DIR / "report_fig_enu_distribution.png"
    e, n, u = enu.T
    h = np.hypot(e, n)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2))
    fig.suptitle("GPS/BDS双频SPP定位误差分布", fontsize=16, fontweight="bold")

    ax = axes[0]
    ax.scatter(e, n, s=1.6, color="#5b45a6", alpha=0.35, edgecolors="none")
    ax.axhline(0, color="#737373", lw=0.8)
    ax.axvline(0, color="#737373", lw=0.8)
    lim = max(float(np.percentile(h, 99.5) * 1.20), 1.0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.set_xlabel("E方向误差 / m")
    ax.set_ylabel("N方向误差 / m")
    ax.text(
        0.03,
        0.96,
        f"HRMS = {rms(h):.3f} m",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#bdbdbd", alpha=0.95),
    )

    ax = axes[1]
    bins = 80
    ax.hist(e, bins=bins, density=True, color="#1f5fbf", alpha=0.52, label=f"E RMS={rms(e):.3f} m")
    ax.hist(n, bins=bins, density=True, color="#198754", alpha=0.52, label=f"N RMS={rms(n):.3f} m")
    ax.hist(u, bins=bins, density=True, color="#c0392b", alpha=0.42, label=f"U RMS={rms(u):.3f} m")
    ax.axvline(0, color="#737373", lw=0.8)
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.set_xlabel("误差 / m")
    ax.set_ylabel("概率密度")
    ax.legend(fontsize=9, frameon=True)
    ax.text(
        0.97,
        0.96,
        f"3D RMS = {rms(np.sqrt(e * e + n * n + u * u)):.3f} m",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#bdbdbd", alpha=0.95),
    )

    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path


def save_sat_pdop(t_hour: np.ndarray, gps: np.ndarray, bds: np.ndarray, pdop: np.ndarray) -> Path:
    fig_path = OUT_DIR / "report_fig_satellite_pdop.png"
    total = gps + bds
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 6.4), sharex=True)
    fig.suptitle("参与解算卫星数与PDOP随时间变化", fontsize=16, fontweight="bold")

    ax = axes[0]
    ax.plot(t_hour, total, color="#333333", lw=1.0, label="总卫星数")
    ax.plot(t_hour, gps, color="#1f5fbf", lw=0.85, alpha=0.9, label="GPS")
    ax.plot(t_hour, bds, color="#198754", lw=0.85, alpha=0.9, label="BDS")
    ax.set_ylabel("卫星数 / 颗")
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(ncol=3, fontsize=9, frameon=True, loc="upper right")
    ax.text(
        0.015,
        0.90,
        f"平均总卫星数={total.mean():.2f}，GPS={gps.mean():.2f}，BDS={bds.mean():.2f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#bdbdbd", alpha=0.95),
    )

    ax = axes[1]
    ax.plot(t_hour, pdop, color="#8a4b08", lw=0.9)
    ax.axhline(pdop.mean(), color="#737373", lw=0.8, ls="--", label=f"平均PDOP={pdop.mean():.3f}")
    ax.set_ylabel("PDOP")
    ax.set_xlabel("观测时间 / h")
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(fontsize=9, frameon=True, loc="upper right")
    ax.set_xlim(t_hour[0], t_hour[-1])

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path


def save_error_pdop_relation(enu: np.ndarray, total_sat: np.ndarray, pdop: np.ndarray) -> Path:
    fig_path = OUT_DIR / "report_fig_error_pdop_relation.png"
    e, n, u = enu.T
    h = np.hypot(e, n)
    abs_u = np.abs(u)

    sample_step = max(len(pdop) // 12000, 1)
    idx = np.arange(0, len(pdop), sample_step)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.1))
    fig.suptitle("定位误差与几何强度指标关系", fontsize=16, fontweight="bold")

    ax = axes[0]
    ax.scatter(pdop[idx], h[idx], s=2.0, color="#1f5fbf", alpha=0.24, edgecolors="none", label="水平误差")
    ax.scatter(pdop[idx], abs_u[idx], s=2.0, color="#c0392b", alpha=0.18, edgecolors="none", label="|U误差|")
    ax.set_xlabel("PDOP")
    ax.set_ylabel("误差绝对值 / m")
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(fontsize=9, frameon=True)
    ax.set_xlim(max(0, np.percentile(pdop, 0.5) * 0.95), np.percentile(pdop, 99.5) * 1.05)
    ax.set_ylim(0, max(np.percentile(h, 99.3), np.percentile(abs_u, 99.3)) * 1.20)

    ax = axes[1]
    sat_values = np.unique(total_sat.astype(int))
    sat_centers, h_median, u_median = [], [], []
    for sv in sat_values:
        mask = total_sat.astype(int) == sv
        if mask.sum() < 30:
            continue
        sat_centers.append(sv)
        h_median.append(np.median(h[mask]))
        u_median.append(np.median(abs_u[mask]))
    ax.plot(sat_centers, h_median, marker="o", ms=3.8, color="#1f5fbf", label="水平误差中位数")
    ax.plot(sat_centers, u_median, marker="s", ms=3.8, color="#c0392b", label="|U误差|中位数")
    ax.set_xlabel("参与解算总卫星数 / 颗")
    ax.set_ylabel("误差中位数 / m")
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(fontsize=9, frameon=True)

    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path


def save_horizontal_up_scatter(enu: np.ndarray) -> Path:
    fig_path = OUT_DIR / "report_fig_horizontal_up_scatter.png"
    e, n, u = enu.T
    h = np.hypot(e, n)

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    ax.scatter(h, u, s=2.2, color="#ff8c00", alpha=0.62, edgecolors="none")
    ax.axhline(0, color="#9e9e9e", lw=0.8)
    ax.grid(True, axis="y", color="#dddddd", lw=0.6)
    ax.set_xlabel("水平误差 / m")
    ax.set_ylabel("U方向误差 / m")
    ax.set_title("水平误差与U方向误差散点图", fontsize=15, fontweight="bold")
    ax.set_xlim(0, max(np.percentile(h, 99.7) * 1.12, 1.0))
    ax.set_ylim(-robust_ylim(u), robust_ylim(u))
    ax.text(
        0.98,
        0.94,
        f"HRMS={rms(h):.3f} m, U RMS={rms(u):.3f} m",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#bdbdbd", alpha=0.94),
    )
    fig.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path


def save_error_cdf(enu: np.ndarray) -> Path:
    fig_path = OUT_DIR / "report_fig_error_cdf.png"
    e, n, u = enu.T
    h = np.hypot(e, n)
    d3 = np.sqrt(e * e + n * n + u * u)

    fig, ax = plt.subplots(figsize=(8.7, 5.3))
    for vals, label, color in [
        (np.abs(e), "|E|", "#1f5fbf"),
        (np.abs(n), "|N|", "#198754"),
        (np.abs(u), "|U|", "#c0392b"),
        (h, "水平误差", "#5b45a6"),
        (d3, "三维误差", "#333333"),
    ]:
        xs = np.sort(vals)
        ys = np.linspace(0, 100, len(xs), endpoint=True)
        ax.plot(xs, ys, lw=1.3, color=color, label=f"{label} 95%={np.percentile(vals, 95):.2f} m")
    ax.set_xlabel("误差绝对值 / m")
    ax.set_ylabel("累计概率 / %")
    ax.set_xlim(0, np.percentile(d3, 99.5) * 1.08)
    ax.set_ylim(0, 100)
    ax.grid(True, color="#dddddd", lw=0.6)
    ax.legend(fontsize=9, frameon=True, loc="lower right")
    ax.set_title("定位误差累计分布", fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fig_path


def write_stats(path: Path, skipped: int, t_hour: np.ndarray, ref: np.ndarray, enu: np.ndarray, gps: np.ndarray, bds: np.ndarray, pdop: np.ndarray) -> Path:
    e, n, u = enu.T
    h = np.hypot(e, n)
    d3 = np.sqrt(e * e + n * n + u * u)
    stats_path = OUT_DIR / "report_precision_stats.md"
    lines = [
        "# SPP结果精度统计",
        "",
        f"- 有效历元数：{len(t_hour)}",
        f"- 跳过残缺/非法行数：{skipped}",
        f"- 连续观测时长：{t_hour[-1] - t_hour[0]:.3f} h",
        "- 参考坐标来源：OEM719 BESTPOSB / RINEX APPROX POSITION XYZ",
        f"- 参考坐标：X={ref[0]:.4f} m，Y={ref[1]:.4f} m，Z={ref[2]:.4f} m",
        f"- 平均GPS卫星数：{gps.mean():.2f}",
        f"- 平均BDS卫星数：{bds.mean():.2f}",
        f"- 平均总卫星数：{(gps + bds).mean():.2f}",
        f"- 平均PDOP：{pdop.mean():.3f}",
        "",
        "| 方向 | 均值/m | 标准差/m | RMS/m | 95%绝对误差/m | 最大绝对误差/m |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, vals in [("E", e), ("N", n), ("U", u)]:
        lines.append(
            f"| {name} | {vals.mean():.3f} | {vals.std(ddof=1):.3f} | {rms(vals):.3f} | "
            f"{percentile_abs(vals, 95):.3f} | {np.max(np.abs(vals)):.3f} |"
        )
    lines += [
        f"| 平面 | - | - | {rms(h):.3f} | {np.percentile(h, 95):.3f} | {h.max():.3f} |",
        f"| 三维 | - | - | {rms(d3):.3f} | {np.percentile(d3, 95):.3f} | {d3.max():.3f} |",
        "",
        "## 生成图件",
        "",
        "- report_fig_enu_timeseries.png",
        "- report_fig_enu_distribution.png",
        "- report_fig_satellite_pdop.png",
        "- report_fig_error_pdop_relation.png",
        "- report_fig_horizontal_up_scatter.png",
        "- report_fig_error_cdf.png",
    ]
    stats_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return stats_path


def main() -> None:
    data, skipped = load_spp_result(INPUT_FILE)
    t = data[:, 0]
    x, y, z = data[:, 1], data[:, 2], data[:, 3]
    gps = data[:, 4]
    bds = data[:, 5]
    pdop = data[:, 6]
    t_cont = unwrap_sec_of_week(t)
    t_hour = (t_cont - t_cont[0]) / 3600.0

    ref, enu = compute_enu(x, y, z)

    outputs = [
        save_enu_timeseries(t_hour, enu),
        save_enu_distribution(enu),
        save_sat_pdop(t_hour, gps, bds, pdop),
        save_error_pdop_relation(enu, gps + bds, pdop),
        save_horizontal_up_scatter(enu),
        save_error_cdf(enu),
        write_stats(INPUT_FILE, skipped, t_hour, ref, enu, gps, bds, pdop),
    ]

    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
