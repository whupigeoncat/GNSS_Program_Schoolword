#include <iostream>
#include <iomanip>
#include <cmath>
#include "RTK_Structs.h"

using namespace std;

static bool NearlyEqual(double a, double b, double tol)
{
	return fabs(a - b) <= tol;
}

int main()
{
	cout << fixed << setprecision(10);

	int pass = 0, fail = 0;
	auto CHECK = [&](bool ok, const char* name) {
		if (ok) { cout << "[PASS] " << name << "\n"; ++pass; }
		else { cout << "[FAIL] " << name << "\n"; ++fail; }
		};

	// ==================== 时间转换测试（按图中参考值） ====================
	cout << "---------------- 时间转换测试 ----------------\n";

	COMMONTIME ct{};
	ct.Year = 2022; ct.Month = 9; ct.Day = 15;
	ct.Hour = 10; ct.Minute = 5; ct.Second = 30.0;

	MJDTIME mjd{};
	GPSTIME gt{};
	COMMONTIME ct_back{};

	bool ok1 = CommonTimeToMJDTime(&ct, &mjd);
	bool ok2 = MJDTimeToGPSTime(&mjd, &gt);
	bool ok3 = GPSTimeToCommonTime(&gt, &ct_back);

	cout << "通用时: " << ct.Year << " " << ct.Month << " " << ct.Day << " "
		<< ct.Hour << " " << ct.Minute << " " << ct.Second << "\n";
	cout << "简化儒略日: " << mjd.Days << " " << mjd.FracDay << "\n";
	cout << "GPST: " << gt.Week << " " << gt.SecOfWeek << "\n";
	cout << "通用时: " << ct_back.Year << " " << ct_back.Month << " " << ct_back.Day << " "
		<< ct_back.Hour << " " << ct_back.Minute << " " << ct_back.Second << "\n";

	// 图中参考：MJD=59837, 0.4204861111; GPST=2227, 381930
	CHECK(ok1 && ok2 && ok3, "时间函数调用成功");
	CHECK(mjd.Days == 59837, "MJD Days == 59837");
	CHECK(NearlyEqual(mjd.FracDay, 0.4204861111, 1e-10), "MJD FracDay");
	CHECK(gt.Week == 2227, "GPST Week == 2227");
	CHECK(NearlyEqual(gt.SecOfWeek, 381930.0, 1e-6), "GPST SecOfWeek");
	CHECK(ct_back.Year == ct.Year && ct_back.Month == ct.Month && ct_back.Day == ct.Day
		&& ct_back.Hour == ct.Hour && ct_back.Minute == ct.Minute
		&& NearlyEqual(ct_back.Second, ct.Second, 1e-6), "通用时往返一致");

	// ==================== 坐标转换测试（按图中参考值） ====================
	cout << "\n---------------- 坐标转换测试 ----------------\n";

	XYZ xyz{};
	xyz.x = -2267807.853;
	xyz.y = 5009320.431;
	xyz.z = 3221020.875;

	BLH blh{};
	XYZ xyz_back{};

	bool ok4 = XYZToBLH(&xyz, &blh, R_WGS84, F_WGS84);
	bool ok5 = BLHToXYZ(&blh, &xyz_back, R_WGS84, F_WGS84);

	cout << "X/m: " << xyz.x << " Y/m: " << xyz.y << " Z/m: " << xyz.z << "\n";
	cout << "B/rad: " << blh.b << " L/rad: " << blh.l << " H/m: " << blh.h << "\n";
	cout << "X/m: " << xyz_back.x << " Y/m: " << xyz_back.y << " Z/m: " << xyz_back.z << "\n";

	// 图中参考：B=0.532827, L=1.995908, H=39.917
	CHECK(ok4 && ok5, "坐标函数调用成功");
	CHECK(NearlyEqual(blh.b, 0.532827, 1e-6), "B/rad");
	CHECK(NearlyEqual(blh.l, 1.995908, 1e-6), "L/rad");
	CHECK(NearlyEqual(blh.h, 39.917, 1e-3), "H/m");
	CHECK(NearlyEqual(xyz_back.x, xyz.x, 1e-3)
		&& NearlyEqual(xyz_back.y, xyz.y, 1e-3)
		&& NearlyEqual(xyz_back.z, xyz.z, 1e-3), "XYZ往返一致");

	// 汇总
	cout << "\n========== SUMMARY ==========\n";
	cout << "PASS = " << pass << "\n";
	cout << "FAIL = " << fail << "\n";
	cout << "=============================\n";

	return (fail == 0) ? 0 : 1;
}