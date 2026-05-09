#include "RTK_Structs.h"

//对流层改正函数
double Hopfield(const double H, const double Elev)
{
	double p, T, e, hw, hd, Kw, Kd, delta_Trop;
	if (H < 0 || H>3e4)	return 0;
	//2. 输入与输出参数
	//	气象参数、测站高度、高度角
	double RH = RH0_Hop * exp(-0.0006396 * (H - H0_Hop));//测站相对湿度
	p = p0_Hop * pow((1 - 0.0000226 * (H - H0_Hop)), 5.225);//测站气压
	T = T0_Hop - 0.0065 * (H - H0_Hop);//测站干温
	e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);
	hw = 11000;
	hd = 40136 + 148.72 * (T0_Hop - 273.16);
	Kw = 155.2 * 1e-7 * (4810 / (T * T)) * e * (hw - H);
	Kd = 155.2 * 1e-7 * (p / T) * (hd - H);
	delta_Trop = Kd / (sin(sqrt(Elev * Elev + 6.25) * Rad)) + Kw / (sin(sqrt(Elev * Elev + 2.25) * Rad));//对流层延迟

	return delta_Trop;
	//	3. 返回值
	//	对流层改正值
	//	4. 内存与分配
	//	5. 容错处理
	//	➢测站高度不在对流层范围，直接输出为0
}

//线性组合探测粗差函数
void DetectOutlier(EPOCHOBS* Obs)
{
	if (Obs == nullptr) return;

	const double GF_THRES = 0.05; // m
	const double MW_THRES = 3.0;  // m

	int i, j;
	// 边界保护，防止 SatNum 异常导致越界
	int satCount = Obs->SatNum;
	if (satCount < 0) satCount = 0;
	if (satCount > MAXCHANNUM) satCount = MAXCHANNUM;

	MWGF CurComObs[MAXCHANNUM] = {};//存放当前历元计算结果
	memset(CurComObs, 0, sizeof(MWGF) * MAXCHANNUM);

	for (i = 0; i < satCount; i++)
	{
		SATOBS& sat = Obs->SatObs[i];
		sat.Valid = false;

		// 1) 观测完整性检查
		if (fabs(sat.P[0]) < 1e-3 || fabs(sat.P[1]) < 1e-3 || fabs(sat.L[0]) < 1e-3 || fabs(sat.L[1]) < 1e-3)
		{
			sat.Valid = false;
			continue;
		}

		// 2) 计算当前历元组合观测值
		CurComObs[i].Sys = sat.System;
		CurComObs[i].Prn = sat.Prn;
		CurComObs[i].GF = sat.L[0] - sat.L[1];
		CurComObs[i].n = 1;

		if (CurComObs[i].Sys == GPS)
		{
			CurComObs[i].MW = (FG1_GPS * sat.L[0] - FG2_GPS * sat.L[1]) / (FG1_GPS - FG2_GPS)
				- (FG1_GPS * sat.P[0] + FG2_GPS * sat.P[1]) / (FG1_GPS + FG2_GPS);
			CurComObs[i].PIF = (FG1_GPS * FG1_GPS * sat.P[0] - FG2_GPS * FG2_GPS * sat.P[1])
				/ (FG1_GPS * FG1_GPS - FG2_GPS * FG2_GPS);
		}
		else if (CurComObs[i].Sys == BDS)
		{
			CurComObs[i].MW = (FG1_BDS * sat.L[0] - FG3_BDS * sat.L[1]) / (FG1_BDS - FG3_BDS)
				- (FG1_BDS * sat.P[0] + FG3_BDS * sat.P[1]) / (FG1_BDS + FG3_BDS);
			CurComObs[i].PIF = (FG1_BDS * FG1_BDS * sat.P[0] - FG3_BDS * FG3_BDS * sat.P[1])
				/ (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS);//伪距IF组合
		}
		else
		{
			continue;
		}

		// 3) 与上一历元同星比较
		bool foundPrev = false; // 是否找到上一历元同星（不代表通过阈值)
		for (j = 0; j < MAXCHANNUM; j++)
		{
			if (Obs->ComObs[j].Prn == sat.Prn && Obs->ComObs[j].Sys == sat.System)
			{
				foundPrev = true;
				break;
			}
		}

		if (foundPrev == false) continue;

		double dGF, dMW;
		dGF = CurComObs[i].GF - Obs->ComObs[j].GF;
		dMW = CurComObs[i].MW - Obs->ComObs[j].MW;

		if (fabs(dGF) > GF_THRES || fabs(dMW) > MW_THRES)
		{
			sat.Valid = false;
		}
		else
		{
			sat.Valid = true;
			CurComObs[i].MW = (Obs->ComObs[j].n * Obs->ComObs[j].MW + CurComObs[i].MW) / (Obs->ComObs[j].n + 1);
			CurComObs[i].n = Obs->ComObs[j].n + 1;
		}
	}

	memcpy(Obs->ComObs, CurComObs, sizeof(MWGF) * MAXCHANNUM);
}
