#include "RTK_Structs.h"
#include <cstdio>
#include <cstring>

//Locktime，判断历元间相位数据连续性，发生中断时有周跳
static double g_BasLock[MAXGPSNUM + MAXBDSNUM][2] = {};
static double g_RovLock[MAXGPSNUM + MAXBDSNUM][2] = {};
static bool   g_HasLock[MAXGPSNUM + MAXBDSNUM] = {};

static double g_PreMW[MAXGPSNUM + MAXBDSNUM] = {};
static double g_PreGF[MAXGPSNUM + MAXBDSNUM] = {};
static bool g_HasComb[MAXGPSNUM + MAXBDSNUM] = {};

static int SatKey(GNSSSys sys, int prn)
{
	return(sys == GPS) ? (prn - 1) : (MAXGPSNUM + prn - 1);
}

//站间单差
void FormSDEpochObs(const EPOCHOBS* Bas_Epk, const EPOCHOBS* Rov_Epk, SDEPOCHOBS* SDObs)
{
	SDObs->SatNum = 0;
	for (int k = 0; k < MAXCHANNUM; k++) SDObs->SdSatObs[k] = SDSATOBS();

	int nBas = Bas_Epk->SatNum;
	if (nBas > MAXCHANNUM) nBas = MAXCHANNUM;
	int nRov = Rov_Epk->SatNum; 
	if (nRov > MAXCHANNUM) nRov = MAXCHANNUM;

	for (int i = 0; i < nRov; i++)
	{
		const SATOBS& rov = Rov_Epk->SatObs[i];
		if (rov.Prn <= 0 || (rov.System != GPS && rov.System != BDS)) continue;
		int jb = -1;
		for (int j = 0; j < nBas; j++)
		{
			if (Bas_Epk->SatObs[j].System == rov.System && Bas_Epk->SatObs[j].Prn == rov.Prn)
			{
				jb = j;
				break;
			}
		}
		if (jb < 0) continue;//不是共同星
		const SATOBS& bas = Bas_Epk->SatObs[jb];//拷贝当前共同星观测值

		SDSATOBS sd;
		sd.Prn = rov.Prn;
		sd.System = rov.System;
		sd.nRov = (short)i;
		sd.nBas = (short)jb;

		bool any = false, phaseOK = true;
		for (int f = 0; f < 2; f++)
		{
			if (fabs(bas.P[f]) < 1e-3 || fabs(bas.L[f]) < 1e-3 || fabs(rov.P[f]) < 1e-3 || fabs(rov.L[f]) < 1e-3)
			{
				continue;
			}

			sd.dP[f] = rov.P[f] - bas.P[f];
			sd.dL[f] = rov.L[f] - bas.L[f];

			any = true;//存在有效单差频率

			if (bas.half[f] == 0 || rov.half[f] == 0)//Parity半周
			{
				phaseOK = false;
			}

			int key = SatKey(sd.System, sd.Prn);
			if(g_HasLock[key])
			{
				if (rov.LockTime[f] < g_RovLock[key][f] - 1e-6) phaseOK = false;
				if (bas.LockTime[f] < g_BasLock[key][f] - 1e-6) phaseOK = false;
			}
			//载波锁定时间，判断失锁
			g_RovLock[key][f] = rov.LockTime[f];
			g_BasLock[key][f] = bas.LockTime[f];
			g_HasLock[key] = true;
		}

		if (!any) continue;
		sd.Valid = phaseOK ? 1 : 0;
		SDObs->SdSatObs[SDObs->SatNum] = sd;
		SDObs->SatNum++;
		
		SDObs->Time = Rov_Epk->Time;
	}
}


//周跳探测
void DetectCycleSlip(SDEPOCHOBS* Obs)
{
	for (int i = 0; i < Obs->SatNum; ++i)
	{
		double f1 = 0.0;
		double f2 = 0.0;

		if (Obs->SdSatObs[i].System != GPS && Obs->SdSatObs[i].System != BDS) continue;

		if (Obs->SdSatObs[i].System == GPS && (Obs->SdSatObs[i].Prn < 1 || Obs->SdSatObs[i].Prn > MAXGPSNUM))
		{
			continue;
		}

		if (Obs->SdSatObs[i].System == BDS && (Obs->SdSatObs[i].Prn < 1 || Obs->SdSatObs[i].Prn > MAXBDSNUM))
		{
			continue;
		}
		
		int key = SatKey(Obs->SdSatObs[i].System, Obs->SdSatObs[i].Prn);//卫星唯一编码

		if (Obs->SdSatObs[i].Valid != 1)
		{
			g_HasComb[key] = false;
			continue;
		}

		if (Obs->SdSatObs[i].System == GPS)
		{
			f1 = FG1_GPS;
			f2 = FG2_GPS;
		}

		else if (Obs->SdSatObs[i].System == BDS)
		{
			f1 = FG1_BDS;
			f2 = FG3_BDS;
		}

		double dMW = 0;
		double dGF = 0;

		if (fabs(Obs->SdSatObs[i].dP[0]) < 1e-3 ||
			fabs(Obs->SdSatObs[i].dP[1]) < 1e-3 ||
			fabs(Obs->SdSatObs[i].dL[0]) < 1e-3 ||
			fabs(Obs->SdSatObs[i].dL[1]) < 1e-3)
		{
			g_HasComb[key] = false;
			continue;//双频数据不全
		}


		Obs->SdCObs[i].MW = (Obs->SdSatObs[i].dL[0] * f1 - Obs->SdSatObs[i].dL[1] * f2) / (f1 - f2) - (Obs->SdSatObs[i].dP[0] * f1 + Obs->SdSatObs[i].dP[1] * f2) / (f1 + f2);
		Obs->SdCObs[i].GF = Obs->SdSatObs[i].dL[0] - Obs->SdSatObs[i].dL[1];

		//第一次遇到这颗卫星（上个历元或之前历元没有），直接建立基准
		if (!g_HasComb[key])
		{
			g_PreMW[key] = Obs->SdCObs[i].MW;
			g_PreGF[key] = Obs->SdCObs[i].GF;
			g_HasComb[key] = true;
			continue;
		}

		dMW = Obs->SdCObs[i].MW - g_PreMW[key];
		dGF = Obs->SdCObs[i].GF - g_PreGF[key];

		if (fabs(dMW) > 3.0 || fabs(dGF) > 0.05)
		{
			Obs->SdSatObs[i].Valid = 0;
		}

		//更新基准
		g_PreMW[key] = Obs->SdCObs[i].MW;
		g_PreGF[key] = Obs->SdCObs[i].GF;
	}
}

//基准星选择
void DetRefSat(const EPOCHOBS* Bas_Epk, const EPOCHOBS* Rov_Epk, SDEPOCHOBS* SDObs, DDCOBS* DDObs)
{
	if (Bas_Epk == nullptr || Rov_Epk == nullptr || SDObs == nullptr ||DDObs == nullptr) return;

	//每个历元参考星重选
	DDObs->RefPrn[0] = -1;  
	DDObs->RefPrn[1] = -1;  

	DDObs->RefPos[0] = -1;  
	DDObs->RefPos[1] = -1;  

	//最高高度角初始化
	double maxElev[2] = {
	-DBL_MAX,   // GPS
	-DBL_MAX    // BDS
	};

	for (int i = 0; i < SDObs->SatNum; i++)
	{
		const SDSATOBS& sd = SDObs->SdSatObs[i];
		int sysIndex = -1;

		//系统确定
		if (sd.System == GPS)
		{
			sysIndex = 0;

			if (sd.Prn < 1 || sd.Prn > MAXGPSNUM)
			{
				continue;
			}
		}
		else if (sd.System == BDS)
		{
			sysIndex = 1;

			if (sd.Prn < 1 || sd.Prn > MAXBDSNUM)
			{
				continue;
			}
		}
		else
		{
			//只处理GPS和BDS
			continue;
		}

		if ((sd.System != BDS && sd.System != GPS)||(sd.Valid!=1)||(sd.nBas < 0 ||sd.nBas >= MAXCHANNUM)||(sd.nRov < 0 || sd.nRov >= MAXCHANNUM))
		{
			continue;
		}

		const SATOBS& bas_sat = Bas_Epk->SatObs[sd.nBas];//观测数据
		const SATOBS& rov_sat = Rov_Epk->SatObs[sd.nRov];

		const SATMIDRES& bas_pos = Bas_Epk->SatPVT[sd.nBas];//卫星位置
		const SATMIDRES& rov_pos = Rov_Epk->SatPVT[sd.nRov];

		if (bas_sat.System != sd.System || bas_sat.Prn != sd.Prn || rov_sat.System != sd.System || rov_sat.Prn != sd.Prn || !bas_sat.Valid || !rov_sat.Valid || !bas_pos.Valid || !rov_pos.Valid)
		{
			continue;
		}

		bool qualityOK = true;


		for (int j = 0; j < 2; j++)
		{
			//两站的伪距和载波都必须存在
			if (!std::isfinite(bas_sat.P[j]) || !std::isfinite(bas_sat.L[j]) || !std::isfinite(rov_sat.P[j]) || !std::isfinite(rov_sat.L[j]) ||
				fabs(bas_sat.P[j]) < 1e-3 || fabs(bas_sat.L[j]) < 1e-3 || fabs(rov_sat.P[j]) < 1e-3 || fabs(rov_sat.L[j]) < 1e-3)
			{
				qualityOK = false;
				break;
			}

			//确认单差值是有限数。
			if (!std::isfinite(sd.dP[j]) || !std::isfinite(sd.dL[j]))
			{
				qualityOK = false;
				break;
			}

			//防止半周
			if (bas_sat.half[j] == 0 || rov_sat.half[j] == 0)
			{
				qualityOK = false;
				break;
			}

			//使用流动站两个频率的CN0作为筛选条件
			if (!std::isfinite(rov_sat.cn0[j]) || rov_sat.cn0[j] <= 35.0)
			{
				qualityOK = false;
				break;
			}

			//站间单差载波需要两站都连续跟踪
			if (!std::isfinite(bas_sat.LockTime[j]) ||
				bas_sat.LockTime[j] <= 6.0 ||
				!std::isfinite(rov_sat.LockTime[j]) ||
				rov_sat.LockTime[j] <= 6.0)
			{
				qualityOK = false;
				break;
			}
		}

		if (!qualityOK)
		{
			continue;
		}

		//高度角最大
		if (rov_pos.Elevation > maxElev[sysIndex])
		{
			maxElev[sysIndex] = rov_pos.Elevation;
			DDObs->RefPrn[sysIndex] = sd.Prn;

			// 保存的是SDObs->SdSatObs[]中的下标i
			DDObs->RefPos[sysIndex] = i;
		}
	}
}

//单个同步历元的RTK处理：SPP初值→基站卫星位置→参考星→浮点解
//返回true表示浮点解成功，结果存放在Raw->FloatRes
bool RTKProcEpoch(RAWDAT* Raw, const ROVERCFGINFO* Cfg, PPRESULT* RovRes)
{
	if (Raw == nullptr || RovRes == nullptr) return false;
	if (!Raw->BasPosValid) return false;//基站坐标无效

	//清理卫星位置
	for (int i = 0; i < MAXCHANNUM; i++)
	{
		Raw->BasEpk.SatPVT[i].Valid = false;
		Raw->RovEpk.SatPVT[i].Valid = false;
	}

	DetectOutlier(&Raw->BasEpk);
	DetectOutlier(&Raw->RovEpk);

	//SPP提供流动站初值，并解出流动站历元卫星位置
	//短基线下，用已知基站坐标作为流动站SPP的近似初值
	for (int k = 0; k < 3; k++)
	{
		Raw->RovEpk.Pos[k] = Raw->BasPos.xyz[k];
	}

	// SPP继续求解流动站的实际位置
	if (!SPP(&Raw->RovEpk, Raw, RovRes)) return false;

	//逐星计算基站历元卫星位置
	for (int i = 0; i < Raw->BasEpk.SatNum && i < MAXCHANNUM; i++)
	{
		const SATOBS& bs = Raw->BasEpk.SatObs[i];
		if (!bs.Valid || (bs.System != GPS && bs.System != BDS))
		{
			continue;
		}

		if (bs.P[0] <= 1e-6) continue;//伪距不可用
		ComputeSatOrbitAtSignalTrans(&Raw->BasEpk, Raw->GpsEph, Raw->BdsEph, Raw->BasPos.xyz, &Raw->BasEpk.SatPVT[i], i, bs.P[0]);
	}

	//参考星选取
	DetRefSat(&Raw->BasEpk, &Raw->RovEpk, &Raw->SdObs, &Raw->DDObs);
	if (Raw->DDObs.RefPrn[0] < 0 && Raw->DDObs.RefPrn[1] < 0) return false;//参考星选取失败

	//相对定位浮点解：1=EKF滤波，2=最小二乘（rtk.cfg第2项控制）
	bool ok = false;
	if (Cfg != nullptr && Cfg->RTKProcMode == 1)
	{
		ok = RTKFilter(&Raw->Ekf, Raw, RovRes, Cfg);
	}
	else
	{
		ok = RTKFloat(Raw, nullptr, RovRes, Cfg);
	}
	//模糊度固定：失败不影响浮点解，bFixed保持false
	if (ok) RTKFix(Raw, Cfg);
	return ok;
}

//单差卫星频率能否构成双差（两侧数据健康）
static bool DDFreqOK(const SDSATOBS& sd, const SDSATOBS& ref, int f)
{
	if (fabs(sd.dP[f]) < 1e-3 || fabs(sd.dL[f]) < 1e-3) return false;
	if (fabs(ref.dP[f]) < 1e-3 || fabs(ref.dL[f]) < 1e-3) return false;
	if (!std::isfinite(sd.dP[f]) || !std::isfinite(sd.dL[f])) return false;
	if (!std::isfinite(ref.dP[f]) || !std::isfinite(ref.dL[f])) return false;
	return true;
}

//相对定位浮点解，双差、最小二乘
bool RTKFloat(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov, const ROVERCFGINFO* Cfg)
{
	(void)Base;

	if (Raw == nullptr || Rov == nullptr) return false;
	if (!Raw->BasPosValid) return false;	//基站已知坐标未就绪
	if (!Rov->IsSuccess) return false;		//SPP失败，无流动站初值

	SDEPOCHOBS* SD = &Raw->SdObs;
	DDCOBS* DD = &Raw->DDObs;
	RTKFLOATRES* Res = &Raw->FloatRes;

	Res->Valid = false;

	int RefIdx[2] = { DD->RefPos[0],DD->RefPos[1] };//参考星在单差中序号
	if (RefIdx[0] < 0 && RefIdx[1] < 0) return false;//参考星未选出

	const double LAM[2][2] = { {WL1_GPS,WL2_GPS},{WL1_BDS,WL3_BDS} };

	double SigCode = 0.3, SigPhase = 0.01;//默认噪声
	if (Cfg != nullptr)//读取配置噪声
	{
		if (Cfg->CodeNoise > 1e-6)
		{
			SigCode = Cfg->CodeNoise;
		}
		if (Cfg->CPNoise > 1e-6)
		{
			SigPhase = Cfg->CPNoise;
		}
	}

	double Xb[3] = { Raw->BasPos.xyz[0],Raw->BasPos.xyz[1],Raw->BasPos.xyz[2] };
	int nSD = SD->SatNum;
	//if (nSD > MAXCHANNUM) nSD = MAXCHANNUM;

	//几何可用性、基站几何距离
	bool geoOK[MAXCHANNUM];
	double rhoB[MAXCHANNUM], rhoR[MAXCHANNUM];
	for (int i = 0; i < MAXCHANNUM; i++)
	{
		geoOK[i] = false;
		rhoB[i] = rhoR[i] = 0.0;
	}

	//逐卫星检查卫星位置
	for (int i = 0; i < nSD; i++)
	{
		const SDSATOBS& sd = SD->SdSatObs[i];
		if (sd.nBas < 0 || sd.nBas >= MAXCHANNUM || sd.nRov < 0 || sd.nRov >= MAXCHANNUM) continue;
		const SATMIDRES& bas = Raw->BasEpk.SatPVT[sd.nBas];
		const SATMIDRES& rov = Raw->RovEpk.SatPVT[sd.nRov];
		if (!bas.Valid || !rov.Valid) continue;
		geoOK[i] = true;

		double dx = Xb[0] - bas.SatPos[0];
		double dy = Xb[1] - bas.SatPos[1];
		double dz = Xb[2] - bas.SatPos[2];
		rhoB[i] = sqrt(dx * dx + dy * dy + dz * dz);//基站站星几何距离
	}

	int AmbCol[MAXCHANNUM][2], PRow[MAXCHANNUM][2], LRow[MAXCHANNUM][2];//模糊度X列号、伪距相位BW行号初始化
	for (int i = 0; i < MAXCHANNUM; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			AmbCol[i][j] = PRow[i][j] = LRow[i][j] = -1;
		}
	}

	int nP[4] = { 0 }, nL[4] = { 0 }, nAmbGrp[4] = { 0 };//0=GPS-L1 1=GPS-L2 2=BDS-B1I 3=BDS-B3I
	int nDDSat[2] = { 0,0 };//各系统双差卫星数

	//统计各系统各频率的双差观测数和模糊度参数数
	for (int i = 0; i < nSD; i++)
	{
		if (!geoOK[i]) continue;
		const SDSATOBS& sd = SD->SdSatObs[i];
		if (sd.Valid != 1) continue;
		int sys = (sd.System == GPS) ? 0 : ((sd.System == BDS) ? 1 : -1);
		if (sys < 0) continue;
		if (RefIdx[sys] < 0 || i == RefIdx[sys]) continue;//参考星不参与

		const SDSATOBS& ref = SD->SdSatObs[RefIdx[sys]];//参考星
		bool used = false;
		for (int j = 0; j < 2; j++)//双频
		{
			if (!DDFreqOK(sd, ref, j)) continue;//构不成双差
			int g = sys * 2 + j;
			nP[g]++;
			nL[g]++;
			nAmbGrp[g]++;
			used = true;
		}
		if (used) nDDSat[sys]++;//这可非参考星至少有一个频率构成双差
	}
	DD->DDSatNum[0] = nDDSat[0];
	DD->DDSatNum[1] = nDDSat[1];
	DD->Sats = nDDSat[0] + nDDSat[1];

	//模糊度列GPS-L1、GPS-L2、BDS-B1I、BDS-B3I顺序，未知矩阵X前3列为xyz修正量，第4列开始是模糊都
	int AmbBase[4];
	AmbBase[0] = 3;
	for (int i = 1; i < 4; i++)
	{
		AmbBase[i] = AmbBase[i - 1] + nAmbGrp[i - 1];//各频率模糊度列开始列数
	}

	int nObs = 0;
	for (int i = 0; i < 4; i++) nObs += nP[i] + nL[i];//双差观测方程数
	int nAmb = AmbBase[3] + nAmbGrp[3] - 3;//双差模糊度参数总数
	int nX = 3 + nAmb;//未知数个数3个坐标改正值+双差模糊度

	if (nX > 3 + MAXCHANNUM * 2) return false;//超限
	if (nObs <= nX) return false;//双差观测方程数量大于未知数，则求解；否则相对定位不成功，函数退出
	
	//分配行号与模糊度列号
	int nextAmb[4] = { 0,0,0,0 };//4频模糊度计数
	int RowBlk[MAXCHANNUM * 4];//每行观测方程属于哪个权阵块
	int nBlk[8] = { 0 };//每个权阵块行行数
	double AmbN0[MAXCHANNUM * 2];//模糊度初值（周）
	short AmbSysArr[MAXCHANNUM * 2], AmbPrnArr[MAXCHANNUM * 2], AmbFreqArr[MAXCHANNUM * 2];
	
	int row = 0;
	for (int i = 0; i < nSD; i++)
	{
		if (!geoOK[i]) continue;
		const SDSATOBS& sd = SD->SdSatObs[i];
		if (sd.Valid != 1) continue;
		int sys = (sd.System == GPS) ? 0 : ((sd.System == BDS) ? 1 : -1);
		if (sys < 0) continue;
		if (RefIdx[sys] < 0 || i == RefIdx[sys]) continue;

		const SDSATOBS& ref = SD->SdSatObs[RefIdx[sys]];
		for (int j = 0; j < 2; j++)
		{
			if (!DDFreqOK(sd, ref, j)) continue;
			int g = sys * 2 + j;

			PRow[i][j] = row;
			RowBlk[row] = g * 2 + 0;
			nBlk[g * 2 + 0]++;
			row++;

			LRow[i][j] = row;
			RowBlk[row] = g * 2 + 1;
			nBlk[g * 2 + 1]++;
			row++;

			//模糊度初值：双差相位-双差伪距在/波长
			int c = AmbBase[g] + nextAmb[g]++;
			AmbCol[i][j] = c;
			AmbN0[c - 3] = ((sd.dL[j] - ref.dL[j]) - (sd.dP[j] - ref.dP[j])) / LAM[sys][j];
			AmbSysArr[c - 3] = sd.System;
			AmbPrnArr[c - 3] = sd.Prn;
			AmbFreqArr[c - 3] = (short)j;
		}
	}

	//组装权阵P
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nObs, nX);
	Eigen::VectorXd W = Eigen::VectorXd::Zero(nObs);
	Eigen::MatrixXd Pmat = Eigen::MatrixXd::Zero(nObs, nObs);
	Eigen::VectorXd Namb = Eigen::VectorXd::Zero(nAmb);//当前浮点模糊度
	for (int i = 0; i < nAmb; i++) Namb(i) = AmbN0[i];

	for (int i = 0; i < nObs; i++)
	{
		for (int j = 0; j < nObs; j++)
		{
			if (RowBlk[i] != RowBlk[j]) continue;
			int blk = RowBlk[i];
			double sig = (blk % 2 == 0) ? SigCode : SigPhase;//伪距/相位噪声
			double n = (double)nBlk[blk];
			Pmat(i, j) = ((i == j) ? n : -1.0) / (2.0 * sig * sig * (n + 1.0));
		}
	}

	//迭代最小二乘
	double Xr[3] = { Rov->Position[0],Rov->Position[1],Rov->Position[2] };

	const int MAX_ITER = 10;
	const double POS_THRES = 1e-3;//收敛阈值(m)

	int iter = 0;
	bool converged = false;
	Eigen::VectorXd vLast = Eigen::VectorXd::Zero(nObs);//最终残差
	Eigen::MatrixXd NLast = Eigen::MatrixXd::Zero(nX, nX);//最终法方程

	for (iter = 0; iter < MAX_ITER; iter++)
	{
		//计算流动站到参考星的几何距离+方向余弦
		double cosR[MAXCHANNUM][3];
		for (int i = 0; i < nSD; i++)
		{
			if (!geoOK[i]) continue;
			const SATMIDRES& rov = Raw->RovEpk.SatPVT[SD->SdSatObs[i].nRov];
			double dx = Xr[0] - rov.SatPos[0];
			double dy = Xr[1] - rov.SatPos[1];
			double dz = Xr[2] - rov.SatPos[2];
			rhoR[i] = sqrt(dx * dx + dy * dy + dz * dz);
			if (rhoR[i] < 1.0) return false;//几何异常
			cosR[i][0] = dx / rhoR[i];
			cosR[i][1] = dy / rhoR[i];
			cosR[i][2] = dz / rhoR[i];
		}

		//组装B与W
		B.setZero();
		W.setZero();
		for (int i = 0; i < nSD; i++)
		{
			if (!geoOK[i]) continue;
			const SDSATOBS& sd = SD->SdSatObs[i];
			int sys = (sd.System == GPS) ? 0 : ((sd.System == BDS) ? 1 : -1);
			if (sys < 0) continue;
			if (RefIdx[sys] < 0 || i == RefIdx[sys]) continue;
			int ri = RefIdx[sys];
			const SDSATOBS& ref = SD->SdSatObs[ri];

			double rhoDD = (rhoR[i] - rhoB[i]) - (rhoR[ri] - rhoB[ri]);//双差几何距离

			for (int j = 0; j < 2; j++)
			{
				int rowP = PRow[i][j];
				if (rowP < 0) continue;
				int rowL = LRow[i][j];
				int colA = AmbCol[i][j];

				for (int k = 0; k < 3; k++)
				{
					double coef = cosR[i][k] - cosR[ri][k];//流动站方向余弦差
					B(rowP, k) = coef;
					B(rowL, k) = coef;
				}
				B(rowL, colA) = LAM[sys][j];//相位行加模糊度序列

				W(rowP) = (sd.dP[j] - ref.dP[j]) - rhoDD;
				W(rowL) = (sd.dL[j] - ref.dL[j]) - rhoDD - LAM[sys][j] * Namb(colA - 3);
			}
		}

		//最小二乘N=B'PB,L=B'PW
		NLast = B.transpose() * Pmat * B;
		Eigen::VectorXd Lvec = B.transpose() * Pmat * W;
		Eigen::LDLT<Eigen::MatrixXd> ldlt(NLast);
		if (ldlt.info() != Eigen::Success) return false;//法方程奇异
		Eigen::VectorXd dx = ldlt.solve(Lvec);
		if (!dx.allFinite()) return false;

		vLast = W - B * dx;//残差

		//更新流动站位置和双差模糊度
		Xr[0] += dx(0);
		Xr[1] += dx(1);
		Xr[2] += dx(2);
		for (int c = 0; c < nAmb; c++) Namb(c) += dx(3 + c);

		//收敛判断
		if (sqrt(dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2)) < POS_THRES)
		{
			converged = true;
			iter++;
			break;
		}
	}
	if (!converged) return false;//迭代未收敛

	//精度评定
	double vtpv = (vLast.transpose() * Pmat * vLast)(0);
	double sigma0 = sqrt(vtpv / (double)(nObs - nX));//单位权中误差
	Eigen::MatrixXd Qxx = NLast.inverse();
	if (!Qxx.allFinite()) return false;

	//保存
	Res->ObsNum = nObs;
	Res->AmbNum = nAmb;
	for (int g = 0; g < 4; g++) Res->AmbNumGrp[g] = nAmbGrp[g];
	for (int c = 0; c < nAmb; c++)
	{
		Res->AmbSys[c] = AmbSysArr[c];
		Res->AmbPrn[c] = AmbPrnArr[c];
		Res->AmbFreq[c] = AmbFreqArr[c];
		Res->FloatAmb[c] = Namb(c);
	}

	Res->dPos[0] = Xr[0] - Xb[0];
	Res->dPos[1] = Xr[1] - Xb[1];
	Res->dPos[2] = Xr[2] - Xb[2];

	XYZ bXYZ{ Xb[0],Xb[1],Xb[2] };
	XYZ rXYZ{ Xr[0],Xr[1],Xr[2] };
	NEU dNeu;
	if (CompEnudPos(&bXYZ, &rXYZ, &dNeu))
	{
		Res->dNeu[0] = dNeu.dN;
		Res->dNeu[1] = dNeu.dE;
		Res->dNeu[2] = dNeu.dU;
	}

	Res->Sigma0 = sigma0;
	Res->PosRMS[0] = sigma0 * sqrt(Qxx(0, 0));
	Res->PosRMS[1] = sigma0 * sqrt(Qxx(1, 1));
	Res->PosRMS[2] = sigma0 * sqrt(Qxx(2, 2));

	for (int i = 0; i < nX; i++)
	{
		for (int j = 0; j < nX; j++)
		{
			Res->Qxx[i * (3 + MAXCHANNUM * 2) + j] = Qxx(i, j);
		}
	}

	Res->IterNum = iter;
	Res->Time = SD->Time;
	Res->Valid = true;

	DD->dPos[0] = Res->dPos[0];
	DD->dPos[1] = Res->dPos[1];
	DD->dPos[2] = Res->dPos[2];
	DD->bFixed = false;//浮点解

	return true;
}

//卡尔曼滤波
static const double KF_SIG_POS = 8.0;//流动站位置的先验标准差
static const double KF_SIG_AMB = 4.0;//模糊度初值标准差m
static const double KF_SIG_AMB_RW = 1e-10;//没有周跳时，模糊度过程噪声设为1E-10

//前一历元模糊度状态中查列号，失败返回-1
static int KFFindAmb(const RTKEKF* Kf, GNSSSys sys, int prn, int freq, int nAmbPrev)
{
	for (int k = 0; k < nAmbPrev; k++)
	{
		if (Kf->AmbSys[k] == (short)sys && Kf->AmbPrn[k] == (short)prn && Kf->AmbFreq[k] == (short)freq) return k;
	}
	return -1;
}

//卡尔曼滤波浮点解，结果存入Raw->FloatRes
bool RTKFilter(RTKEKF* Kf, RAWDAT* Raw, PPRESULT* Rov, const ROVERCFGINFO* Cfg)
{
	const int ND = 3 + MAXCHANNUM * 2;//状态最大维

	if (Kf == nullptr || Raw == nullptr || Rov == nullptr) return false;
	if (!Raw->BasPosValid) return false;
	if (!Rov->IsSuccess) return false;//SPP失败

	SDEPOCHOBS* SD = &Raw->SdObs;
	DDCOBS* DD = &Raw->DDObs;
	RTKFLOATRES* Res = &Raw->FloatRes;
	Res->Valid = false;

	int RefIdx[2] = { DD->RefPos[0],DD->RefPos[1] };//参考星单差序号
	if (RefIdx[0] < 0 && RefIdx[1] < 0) return false;//参考星未选出

	const double LAM[2][2] = { { WL1_GPS,WL2_GPS },{ WL1_BDS,WL3_BDS } };

	double SigCode = 0.3, SigPhase = 0.01;//默认噪声
	if (Cfg != nullptr)
	{
		if (Cfg->CodeNoise > 1e-6) SigCode = Cfg->CodeNoise;
		if (Cfg->CPNoise > 1e-6) SigPhase = Cfg->CPNoise;
	}

	double Xb[3] = { Raw->BasPos.xyz[0],Raw->BasPos.xyz[1],Raw->BasPos.xyz[2] };
	int nSD = SD->SatNum;
	if (nSD > MAXCHANNUM) nSD = MAXCHANNUM;

	//复用RTKFloat,几何可用性+站星几何距离
	bool geoOK[MAXCHANNUM] = {};
	double rhoB[MAXCHANNUM] = {};
	for (int i = 0; i < nSD; i++)
	{
		const SDSATOBS& sd = SD->SdSatObs[i];
		if (sd.nBas < 0 || sd.nBas >= MAXCHANNUM || sd.nRov < 0 || sd.nRov >= MAXCHANNUM) continue;
		const SATMIDRES& bas = Raw->BasEpk.SatPVT[sd.nBas];
		const SATMIDRES& rov = Raw->RovEpk.SatPVT[sd.nRov];
		if (!bas.Valid || !rov.Valid) continue;
		geoOK[i] = true;

		double dx = Xb[0] - bas.SatPos[0];
		double dy = Xb[1] - bas.SatPos[1];
		double dz = Xb[2] - bas.SatPos[2];
		rhoB[i] = sqrt(dx * dx + dy * dy + dz * dz);//基站坐标固定
	}

	//复用RTKFloat，两遍扫描，统计双差观测数+模糊度参数数，行号+模糊度序列号+初值
	int AmbCol[MAXCHANNUM][2], PRow[MAXCHANNUM][2], LRow[MAXCHANNUM][2];
	for (int i = 0; i < MAXCHANNUM; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			AmbCol[i][j] = PRow[i][j] = LRow[i][j] = -1;
		}
	}

	int nP[4] = { 0 }, nL[4] = { 0 }, nAmbGrp[4] = { 0 };
	int nDDSat[2] = { 0,0 };

	for (int i = 0; i < nSD; i++)
	{
		if (!geoOK[i]) continue;
		const SDSATOBS& sd = SD->SdSatObs[i];
		if (sd.Valid != 1) continue;
		int sys = (sd.System == GPS) ? 0 : ((sd.System == BDS) ? 1 : -1);
		if (sys < 0) continue;
		if (RefIdx[sys] < 0 || i == RefIdx[sys]) continue;//参考星不参与

		const SDSATOBS& ref = SD->SdSatObs[RefIdx[sys]];
		bool used = false;
		for (int j = 0; j < 2; j++)
		{
			if (!DDFreqOK(sd, ref, j)) continue;
			int g = sys * 2 + j;
			nP[g]++;
			nL[g]++;
			nAmbGrp[g]++;
			used = true;
		}
		if (used) nDDSat[sys]++;
	}
	DD->DDSatNum[0] = nDDSat[0];
	DD->DDSatNum[1] = nDDSat[1];
	DD->Sats = nDDSat[0] + nDDSat[1];

	int AmbBase[4];
	AmbBase[0] = 3;
	for (int i = 1; i < 4; i++)
	{
		AmbBase[i] = AmbBase[i - 1] + nAmbGrp[i - 1];
	}

	int nObs = 0;
	for (int i = 0; i < 4; i++) nObs += nP[i] + nL[i];
	int nAmb = AmbBase[3] + nAmbGrp[3] - 3;
	int nX = 3 + nAmb;//未知数：3个坐标改正值+双差模糊度

	if (nX > 3 + MAXCHANNUM * 2) return false;//保险
	if (nObs <= nX) return false;//双差观测方程数量大于未知数，则求解；否则相对定位不成功，函数退出

	//分配行号与模糊度列号
	int nextAmb[4] = { 0,0,0,0 };
	int RowBlk[MAXCHANNUM * 4];//每行所属权阵块
	int nBlk[8] = { 0 };
	double AmbN0[MAXCHANNUM * 2];//模糊度初值（周）
	short AmbSysArr[MAXCHANNUM * 2], AmbPrnArr[MAXCHANNUM * 2], AmbFreqArr[MAXCHANNUM * 2];

	int row = 0;
	for (int i = 0; i < nSD; i++)
	{
		if (!geoOK[i]) continue;
		const SDSATOBS& sd = SD->SdSatObs[i];
		if (sd.Valid != 1) continue;
		int sys = (sd.System == GPS) ? 0 : ((sd.System == BDS) ? 1 : -1);
		if (sys < 0) continue;
		if (RefIdx[sys] < 0 || i == RefIdx[sys]) continue;

		const SDSATOBS& ref = SD->SdSatObs[RefIdx[sys]];
		for (int j = 0; j < 2; j++)
		{
			if (!DDFreqOK(sd, ref, j)) continue;
			int g = sys * 2 + j;

			PRow[i][j] = row;
			RowBlk[row] = g * 2 + 0;
			nBlk[g * 2 + 0]++;
			row++;

			LRow[i][j] = row;
			RowBlk[row] = g * 2 + 1;
			nBlk[g * 2 + 1]++;
			row++;

			//模糊度初值：双差相位-双差伪距在/波长
			int c = AmbBase[g] + nextAmb[g]++;
			AmbCol[i][j] = c;
			AmbN0[c - 3] = ((sd.dL[j] - ref.dL[j]) - (sd.dP[j] - ref.dP[j])) / LAM[sys][j];
			AmbSysArr[c - 3] = sd.System;
			AmbPrnArr[c - 3] = sd.Prn;
			AmbFreqArr[c - 3] = (short)j;
		}
	}

	//滤波时间更新，构造预报值X0/P0，数据中断超30s按首次滤波处理（重新初始化）
	bool reinit = !Kf->IsInit;
	if (!reinit && Kf->Time.Week > 0)
	{
		if (fabs(GetDiffTime(&SD->Time, &Kf->Time)) > 30.0) reinit = true;
	}
	int nAmbPrev = reinit ? 0 : Kf->nAmbs;

	Eigen::VectorXd X0 = Eigen::VectorXd::Zero(nX);
	Eigen::MatrixXd P0 = Eigen::MatrixXd::Zero(nX, nX);

	for (int k = 0; k < 3; k++)
	{
		X0(k) = Rov->Position[k];
		P0(k, k) = KF_SIG_POS * KF_SIG_POS;
	}

	Eigen::VectorXd Xa = Eigen::VectorXd::Zero(nAmbPrev);
	Eigen::MatrixXd Pa = Eigen::MatrixXd::Zero(nAmbPrev, nAmbPrev);
	for (int k = 0; k < nAmbPrev; k++)
	{
		Xa(k) = Kf->X[3 + k];
		for (int j = 0; j < nAmbPrev; j++) Pa(k, j) = Kf->P[(3 + k) * ND + (3 + j)];
	}

	Eigen::MatrixXd T = Eigen::MatrixXd::Zero(nAmb, nAmbPrev);
	for (int c = 0; c < nAmb; c++)
	{
		int sys = (AmbSysArr[c] == GPS) ? 0 : 1;
		int io = KFFindAmb(Kf, (GNSSSys)AmbSysArr[c], AmbPrnArr[c], AmbFreqArr[c], nAmbPrev);
		bool refChanged = (Kf->DDObs.RefPrn[sys] < 0 || Kf->DDObs.RefPrn[sys] != DD->RefPrn[sys]);

		if (!refChanged)
		{
			if (io >= 0) { T(c, io) = 1.0; continue; }//连续跟踪：模糊度状态不变
		}
		else
		{
			//基准星变化：N(j-r') = N(j-r) - N(r'-r)，协方差经T重整（课件"模糊度及其协方差的重整"）
			int ir = (DD->RefPrn[sys] > 0) ? KFFindAmb(Kf, (GNSSSys)AmbSysArr[c], DD->RefPrn[sys], AmbFreqArr[c], nAmbPrev) : -1;
			if (io >= 0 && ir >= 0) { T(c, io) = 1.0; T(c, ir) = -1.0; continue; }
		}
		//走到这里=新升卫星/周跳后回归/无法重整 → 保持T该行全0，下面按初值重新初始化
	}

	Eigen::VectorXd Xa2 = Eigen::VectorXd::Zero(nAmb);
	Eigen::MatrixXd Pa2 = Eigen::MatrixXd::Zero(nAmb, nAmb);
	if (nAmbPrev > 0)
	{
		Xa2 = T * Xa;
		Pa2 = T * Pa * T.transpose();
	}
	for (int c = 0; c < nAmb; c++)
	{
		int sys = (AmbSysArr[c] == GPS) ? 0 : 1;
		double lam = LAM[sys][AmbFreqArr[c]];
		if (nAmbPrev == 0 || T.row(c).isZero(1e-12))//新初始化列
		{
			Xa2(c) = AmbN0[c];
			for (int d = 0; d < nAmb; d++) { Pa2(c, d) = Pa2(d, c) = 0.0; }
			Pa2(c, c) = (KF_SIG_AMB / lam) * (KF_SIG_AMB / lam);//(3~6m/λ)^2，注意模糊度单位是周
		}
		else
		{
			Pa2(c, c) += KF_SIG_AMB_RW;//无周跳：过程噪声1e-10（课件）
		}
	}

	for (int c = 0; c < nAmb; c++) X0(3 + c) = Xa2(c);
	P0.bottomRightCorner(nAmb, nAmb) = Pa2;//位置与模糊度间协方差置0（位置由SPP独立预报）

	//滤波测量更新
	//预报位置处的流动站几何距离与方向余弦（EKF在预报点线性化，不迭代）
	double rhoR[MAXCHANNUM] = {};
	double cosR[MAXCHANNUM][3] = {};
	for (int i = 0; i < nSD; i++)
	{
		if (!geoOK[i]) continue;
		const SATMIDRES& rov = Raw->RovEpk.SatPVT[SD->SdSatObs[i].nRov];
		double dx = X0(0) - rov.SatPos[0];
		double dy = X0(1) - rov.SatPos[1];
		double dz = X0(2) - rov.SatPos[2];
		rhoR[i] = sqrt(dx * dx + dy * dy + dz * dz);
		if (rhoR[i] < 1.0) return false;//几何异常
		cosR[i][0] = dx / rhoR[i];
		cosR[i][1] = dy / rhoR[i];
		cosR[i][2] = dz / rhoR[i];
	}

	//组装B、W（与RTKFloat区别：模糊度用预报值X0）
	Eigen::MatrixXd B = Eigen::MatrixXd::Zero(nObs, nX);
	Eigen::VectorXd W = Eigen::VectorXd::Zero(nObs);
	for (int i = 0; i < nSD; i++)
	{
		if (!geoOK[i]) continue;
		const SDSATOBS& sd = SD->SdSatObs[i];
		int sys = (sd.System == GPS) ? 0 : ((sd.System == BDS) ? 1 : -1);
		if (sys < 0) continue;
		if (RefIdx[sys] < 0 || i == RefIdx[sys]) continue;
		int ri = RefIdx[sys];
		const SDSATOBS& ref = SD->SdSatObs[ri];

		double rhoDD = (rhoR[i] - rhoB[i]) - (rhoR[ri] - rhoB[ri]);//双差几何距离

		for (int j = 0; j < 2; j++)
		{
			int rowP = PRow[i][j];
			if (rowP < 0) continue;
			int rowL = LRow[i][j];
			int colA = AmbCol[i][j];

			for (int k = 0; k < 3; k++)
			{
				double coef = cosR[i][k] - cosR[ri][k];//流动站方向余弦差
				B(rowP, k) = coef;
				B(rowL, k) = coef;
			}
			B(rowL, colA) = LAM[sys][j];//相位行加模糊度列

			W(rowP) = (sd.dP[j] - ref.dP[j]) - rhoDD;
			W(rowL) = (sd.dL[j] - ref.dL[j]) - rhoDD - LAM[sys][j] * X0(colA);
		}
	}

	//观测噪声阵R：每块cov=2σ²(I+J)
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(nObs, nObs);
	for (int i = 0; i < nObs; i++)
	{
		for (int j = 0; j < nObs; j++)
		{
			if (RowBlk[i] != RowBlk[j]) continue;
			double sig = (RowBlk[i] % 2 == 0) ? SigCode : SigPhase;//偶数块=伪距，奇数块=相位
			R(i, j) = 2.0 * sig * sig * ((i == j) ? 2.0 : 1.0);
		}
	}

	//滤波测量更新（课件P33-34）
	Eigen::MatrixXd S = B * P0 * B.transpose() + R;
	Eigen::LDLT<Eigen::MatrixXd> ldltS(S);
	if (ldltS.info() != Eigen::Success) return false;
	Eigen::MatrixXd K = ldltS.solve(Eigen::MatrixXd(B * P0)).transpose();//K = P0·Bᵀ·S⁻¹，nX×nObs
	if (!K.allFinite()) return false;
	Eigen::VectorXd X = X0 + K * W;
	Eigen::MatrixXd A = Eigen::MatrixXd::Identity(nX, nX) - K * B;
	Eigen::MatrixXd Pm = A * P0 * A.transpose() + K * R * K.transpose();
	if (!X.allFinite() || !Pm.allFinite()) return false;

	//后验残差质量指标
	Eigen::VectorXd v = W - B * (X - X0);
	Eigen::MatrixXd Rinv = R.inverse();
	double sigma0 = sqrt((v.transpose() * Rinv * v)(0) / (double)(nObs - nX));

	//保存后验状态，供下一历元时间更新（失败提前return）
	for (int k = 0; k < nX; k++) Kf->X[k] = X(k);
	for (int k = 0; k < nX; k++)
		for (int j = 0; j < nX; j++)
			Kf->P[k * ND + j] = Pm(k, j);
	for (int c = 0; c < nAmb; c++)
	{
		Kf->AmbSys[c] = AmbSysArr[c];
		Kf->AmbPrn[c] = AmbPrnArr[c];
		Kf->AmbFreq[c] = AmbFreqArr[c];
	}
	Kf->nAmbs = nAmb;
	Kf->DDObs = *DD;//记录本历元参考星，供下历元判断基准星变化
	Kf->Time = SD->Time;
	Kf->IsInit = true;

	//输出（与RTKFloat同格式）
	Res->ObsNum = nObs;
	Res->AmbNum = nAmb;
	for (int g = 0; g < 4; g++) Res->AmbNumGrp[g] = nAmbGrp[g];
	for (int c = 0; c < nAmb; c++)
	{
		Res->AmbSys[c] = AmbSysArr[c];
		Res->AmbPrn[c] = AmbPrnArr[c];
		Res->AmbFreq[c] = AmbFreqArr[c];
		Res->FloatAmb[c] = X(3 + c);
	}

	Res->dPos[0] = X(0) - Xb[0];
	Res->dPos[1] = X(1) - Xb[1];
	Res->dPos[2] = X(2) - Xb[2];

	XYZ bXYZ{ Xb[0],Xb[1],Xb[2] };
	XYZ rXYZ{ X(0),X(1),X(2) };
	NEU dNeu;
	if (CompEnudPos(&bXYZ, &rXYZ, &dNeu))
	{
		Res->dNeu[0] = dNeu.dN;
		Res->dNeu[1] = dNeu.dE;
		Res->dNeu[2] = dNeu.dU;
	}

	Res->Sigma0 = sigma0;
	Res->PosRMS[0] = sqrt(Pm(0, 0));//滤波协方差直接给出形式精度，应逐历元收缩
	Res->PosRMS[1] = sqrt(Pm(1, 1));
	Res->PosRMS[2] = sqrt(Pm(2, 2));

	for (int i = 0; i < nX; i++)
		for (int j = 0; j < nX; j++)
			Res->Qxx[i * ND + j] = Pm(i, j);

	Res->IterNum = 1;//KF无迭代
	Res->Time = SD->Time;
	Res->Valid = true;

	DD->dPos[0] = Res->dPos[0];
	DD->dPos[1] = Res->dPos[1];
	DD->dPos[2] = Res->dPos[2];
	DD->bFixed = false;//浮点解

	return true;
}


//进行最小二乘或KF解算，得到基线分量和浮点模糊度
//② 采用Z变换，对模糊度搜索空间进行重新参数化，降低浮点模糊度之间的相
//关性（整数变换、降相关）
//③ 联合使用序贯条件最小二乘平差和一个离散搜索策略来估算整数模糊度
//④ 通过逆变换Z
//- 1将模糊度重新转换回原来的模糊度空间（所给出基线向量的
//空间）
//⑤ 将整数模糊度作为整数固定，确定最终的基线分量
bool RTKFix(RAWDAT* Raw, const ROVERCFGINFO* Cfg)
{
	const int ND = 3 + MAXCHANNUM * 2;
	
	if (Raw == nullptr) return false;
	RTKFLOATRES* Res = &Raw->FloatRes;
	DDCOBS* DD = &Raw->DDObs;

	DD->bFixed = false;
	if (!Res->Valid) return false;

	int n = Res->AmbNum;
	if (n<4 || n>MAXCHANNUM * 2) return false;//参数太少，不固定

	//浮点模糊度、协方差块
	Eigen::VectorXd a(n);
	Eigen::MatrixXd Qa = Eigen::MatrixXd::Zero(n, n);
	Eigen::MatrixXd Qba = Eigen::MatrixXd::Zero(3, n);
	Eigen::MatrixXd Qb = Eigen::MatrixXd::Zero(3, 3);

	for (int i = 0; i < n; i++)
	{
		a(i) = Res->FloatAmb[i];
		for (int j = 0; j < n; j++)
		{
			Qa(i, j) = Res->Qxx[(3 + i) * ND + (3 + j)];
		}
		for (int k = 0; k < 3; k++)
		{
			Qba(k, i) = Res->Qxx[k * ND + (3 + i)];
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Qb(i, j) = Res->Qxx[i * ND + j];
		}
	}

	//LAMBDA搜索（F列主序：前n最优、后n次优；s为对应残差平方和）
	double F[MAXCHANNUM * 4] = {};
	double s[2] = { 0.0,0.0 };
	if (lambda(n, 2, a.data(), Qa.data(), F, s) != 0) return false;//搜索失败

	//Ratio检验（s[0]=最优、s[1]=次优，内部升序排序）
	double Ratio = (s[0] > 1e-12) ? (s[1] / s[0]) : 0.0;
	DD->Ratio = Ratio;

	double thres = (Cfg != nullptr && Cfg->RatioThres > 1.0) ? Cfg->RatioThres : 30;
	if (Ratio < thres) return false;//检验失败，保持浮点解

	//固定解
	Eigen::VectorXd d = a - Eigen::Map<Eigen::VectorXd>(F, n);
	Eigen::LDLT<Eigen::MatrixXd> ldltQa(Qa);
	if (ldltQa.info() != Eigen::Success) return false;
	Eigen::VectorXd corr = Qba * ldltQa.solve(d);
	Eigen::MatrixXd QbFix = Qb - Qba * ldltQa.solve(Qba.transpose());
	if (!corr.allFinite() || !QbFix.allFinite()) return false;

	//保存最优、次优整数模糊度
	for (int i = 0; i < n; i++)
	{
		DD->FixedAmb[i] = F[i];
		DD->FixedAmb[n + i] = F[n + i];
	}
	DD->ResAmb[0] = s[0];
	DD->ResAmb[1] = s[1];

	DD->dPos[0] = Res->dPos[0] - corr(0);
	DD->dPos[1] = Res->dPos[1] - corr(1);
	DD->dPos[2] = Res->dPos[2] - corr(2);

	//固定解ENU
	XYZ bXYZ{ Raw->BasPos.xyz[0],Raw->BasPos.xyz[1],Raw->BasPos.xyz[2] };
	XYZ rXYZ{ bXYZ.x + DD->dPos[0], bXYZ.y + DD->dPos[1], bXYZ.z + DD->dPos[2] };
	NEU dNeu;
	if (CompEnudPos(&bXYZ, &rXYZ, &dNeu))
	{
		DD->dNeu[0] = dNeu.dN;
		DD->dNeu[1] = dNeu.dE;
		DD->dNeu[2] = dNeu.dU;
	}
	BLH bBLH;
	if (XYZToBLH(&bXYZ, &bBLH, R_WGS84, F_WGS84))
	{
		Eigen::Matrix3d RNeu;
		if (BLHToNEUMat(&bBLH, &RNeu))
		{
			Eigen::Matrix3d Qneu = RNeu * QbFix * RNeu.transpose();
			DD->FixRMS[0] = (float)sqrt(Qneu(0, 0) + Qneu(1, 1));//水平
			DD->FixRMS[1] = (float)sqrt(Qneu(2, 2));//垂直
		}
	}

	DD->bFixed = true;
	return true;
}