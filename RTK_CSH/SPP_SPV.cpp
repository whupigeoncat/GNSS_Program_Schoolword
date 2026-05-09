#include"RTK_Structs.h"
#include <vector>
#include <algorithm>

//卫星时间修正
static void NormalizedGPSTime(GPSTIME* t)
{
	if (t == nullptr) return;

	while (t->SecOfWeek >= 604800.0)
	{
		t->SecOfWeek -= 604800.0;
		++t->Week;
	}
	while (t->SecOfWeek < 0.0)
	{
		t->SecOfWeek += 604800.0;
		--t->Week;
	}
}

//信号发射时刻卫星位置（每卫星）
bool ComputeSatOrbitAtSignalTrans(const EPOCHOBS* Epk, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, double RcvPos[3], SATMIDRES* MidRes, int satIndex, double Pobs)
{
	if (Epk == nullptr || GPSEph == nullptr || BDSEph == nullptr) return false;
	if (MidRes == nullptr) return false;
	if (satIndex < 0 || satIndex >= Epk->SatNum || satIndex >= MAXCHANNUM) return false;
	MidRes->Valid = false;

	const SATOBS* sat = &Epk->SatObs[satIndex];//单个卫星

	if (sat == nullptr) return false;

	if (sat->System == GPS && (sat->Prn < 1 || sat->Prn > MAXGPSNUM)) return false;
	if (sat->System == BDS && (sat->Prn < 1 || sat->Prn > MAXBDSNUM)) return false;

	if (Pobs <= 1e-6) return false;//伪距观测值

	double dtj = 0.0;
	GPSTIME ts = Epk->Time;
	SATMIDRES clkMid{};
	for (int k = 0; k < 2; k++)
	{
		ts = Epk->Time;
		ts.SecOfWeek -= (Pobs / C_Light + dtj);
		NormalizedGPSTime(&ts);

		if (!CompSatClkOff(sat->Prn, sat->System, &ts, GPSEph, BDSEph, &clkMid)) return false;

		dtj = clkMid.SatClkOft;
		ts.SecOfWeek -= dtj;
		NormalizedGPSTime(&ts);
	}
	// 4) 在最终发射时刻 ts 计算卫星位置、速度、钟差
	bool ok = false;
	if (sat->System == GPS)
	{
		ok = CompGPSSatPVT(sat->Prn, &ts, &GPSEph[sat->Prn - 1], MidRes);
	}
	if (sat->System == BDS)
	{
		ok = CompBDSSatPVT(sat->Prn, &ts, &BDSEph[sat->Prn - 1], MidRes);
	}
	if (!ok) return false;

	// 5) 地球自转改正（Sagnac）
	const double dx = MidRes->SatPos[0] - RcvPos[0];
	const double dy = MidRes->SatPos[1] - RcvPos[1];
	const double dz = MidRes->SatPos[2] - RcvPos[2];

	const double rho = sqrt(dx * dx + dy * dy + dz * dz);//几何距离
	const double dt1 = rho / C_Light;//信号传播时间
	const double ang = Omega_WGS * dt1;

	const double x = MidRes->SatPos[0];
	const double y = MidRes->SatPos[1];

	MidRes->SatPos[0] = cos(ang) * x + sin(ang) * y;
	MidRes->SatPos[1] = -sin(ang) * x + cos(ang) * y;

	MidRes->Valid = true;

	return true;
}

bool SPP(EPOCHOBS* Epoch, RAWDAT* Raw, PPRESULT* Result)
{
	if (Epoch == nullptr || Raw == nullptr || Result == nullptr) return false;

	//上一历元结果作为本次的初值

	//清空输出结果的地址，保留历元时间
	*Result = PPRESULT();
	Result->Time = Epoch->Time;

	//双系统分别估算钟差
	double bGPS = 0.0;
	double bBDS = 0.0;

	//位置初值（优先上一历元定位结果，其次接收机粗定位）
	Eigen::Vector3d xr(0.0, 0.0, 0.0);//初始位置设定
	if (Result->IsSuccess)
	{
		xr << Result->Position[0], Result->Position[1], Result->Position[2];
		bGPS = C_Light * Result->RcvClkOft[0];
		bBDS = C_Light * Result->RcvClkOft[1];
	}
	else if (fabs(Epoch->Pos[0]) + fabs(Epoch->Pos[1]) + fabs(Epoch->Pos[2]) > 1.0)
	{
		xr << Epoch->Pos[0], Epoch->Pos[1], Epoch->Pos[2];
	}

	const int MAX_ITER = 40;
	const double POS_EPS = 1e-4;
	const double CLK_EPS = 1e-4;

	for (int iter = 0; iter < MAX_ITER; ++iter)
	{
		int nObs = 0;
		int ngps = 0;
		int nbds = 0;
		int nall = 0;

		Eigen::MatrixXd B(MAXCHANNUM, 5);
		Eigen::VectorXd w(MAXCHANNUM);
		B.setZero();
		w.setZero();

		for (int i = 0; i < MAXCHANNUM; ++i)//按卫星一个一个处理观测值
		{
			SATOBS& sat = Epoch->SatObs[i];
			if (sat.System != GPS && sat.System != BDS) continue;
			if (sat.Prn <= 0) continue;

			double rcvPos[3] = { xr(0), xr(1), xr(2) };
			SATMIDRES mid;
			//双频无电离层组合伪距观测方程
			double Pobs = 0.0, Trop = 0.0, tgd_m = 0.0;
			if (sat.System == GPS)
				Pobs = (FG12R2 * sat.P[0] - sat.P[1]) / (FG12R2 - 1.0);
			else
				Pobs = (FC13R2 * sat.P[0] - sat.P[1]) / (FC13R2 - 1.0);

			if (ComputeSatOrbitAtSignalTrans(Epoch, Raw->GpsEph, Raw->BdsEph, rcvPos, &mid, i, Pobs))
			{
				Epoch->SatPVT[i] = mid;//卫星位置速度

				XYZ xr_xyz{ xr(0),xr(1),xr(2) };
				XYZ xs_xyz{ mid.SatPos[0],mid.SatPos[1] ,mid.SatPos[2] };
				BLH xr_blh;
				double elev, az;

				XYZToBLH(&xr_xyz, &xr_blh, R_WGS84, F_WGS84);

				if (CompSatEIAz(&xr_xyz, &xs_xyz, &elev, &az))
				{
					double elev_deg = elev / Rad;
					Trop = Hopfield(xr_blh.h, elev_deg);//对流层延迟
				}
			}
			else
			{
				Epoch->SatPVT[i].Valid = false;
				continue;
			}

			//下一步线性化
			if (mid.Valid)//tgd群延迟
			{
				++nall;
				if (sat.System == GPS)
				{
					tgd_m = 0.0; // 按 PPT 简化
					++ngps;
				}
				else if (sat.System == BDS)
				{
					tgd_m = C_Light * (FG1_BDS * FG1_BDS) / (FG1_BDS * FG1_BDS - FG3_BDS * FG3_BDS) * mid.Tgd1;
					++nbds;
				}
				else
				{
					tgd_m = 0.0;
				}
			}
			else continue;

			double dx = xr(0) - mid.SatPos[0];
			double dy = xr(1) - mid.SatPos[1];
			double dz = xr(2) - mid.SatPos[2];
			double rho = sqrt(dx * dx + dy * dy + dz * dz);//几何距离
			if (rho < 1.0) continue;

			double l = dx / rho;
			double m = dy / rho;
			double n = dz / rho;

			double sat_dti_m = C_Light * mid.SatClkOft;//卫星钟差
			double bSys = (sat.System == GPS) ? bGPS : bBDS;//接收机钟差

			double wi = Pobs - rho + sat_dti_m - bSys - Trop - tgd_m;//残差

			B(nObs, 0) = l;
			B(nObs, 1) = m;
			B(nObs, 2) = n;
			B(nObs, 3) = (sat.System == GPS) ? 1.0 : 0.0;
			B(nObs, 4) = (sat.System == BDS) ? 1.0 : 0.0;
			w(nObs) = wi;
			++nObs;
		}
		//统计参与定位的各系统卫星数和所有卫星数
		Result->GPSSatNum = ngps;
		Result->BDSSatNum = nbds;
		Result->AllSatNum = nall;

		Eigen::MatrixXd B_head = B.topRows(nObs);      // nObs × 5
		Eigen::VectorXd w_head = w.head(nObs);         // nObs × 1

		// 权矩阵只需 nObs × nObs
		Eigen::MatrixXd P = Eigen::MatrixXd::Identity(nObs, nObs);
		Eigen::MatrixXd N = B_head.transpose() * P * B_head;
		Eigen::VectorXd W = B_head.transpose() * P * w_head;

		//如果GPS或BDS卫星数量为0，重构法方程矩阵N和Y
		Eigen::MatrixXd N_new;
		Eigen::VectorXd W_new;
		double s;
		// 只保留 GPS=0 的情况举例
		if (ngps == 0 && nbds > 0)
		{
			N_new.resize(4, 4);
			W_new.resize(4);

			N_new.block(0, 0, 3, 3) = N.block(0, 0, 3, 3);
			N_new.block(0, 3, 3, 1) = N.col(4).head(3);
			N_new.block(3, 0, 1, 3) = N.row(4).head(3);
			N_new(3, 3) = N(4, 4);

			W_new << W.head(3), W(4);
			s = 4;
		}
		else if (nbds == 0 && ngps > 0)
		{
			N_new.resize(4, 4);
			W_new.resize(4);

			N_new = N.block(0, 0, 4, 4);
			W_new = W.head(4);
			s = 4;
		}
		else if (nbds > 0 && ngps > 0)
		{
			N_new = N;
			W_new = W;
			s = 5;
		}
		else return false;

		//最小二乘求解
		Eigen::VectorXd X = N_new.ldlt().solve(W_new);

		Eigen::MatrixXd A_new;
		int n_params = 0;

		if (ngps == 0 && nbds > 0)
		{
			n_params = 4;
			A_new.resize(nObs, 4);
			A_new.block(0, 0, nObs, 3) = B_head.block(0, 0, nObs, 3);
			A_new.col(3) = B_head.col(4);
		}
		else if (nbds == 0 && ngps > 0)
		{
			n_params = 4;
			A_new = B_head.block(0, 0, nObs, 4);
		}
		else
		{
			n_params = 5;
			A_new = B_head;
		}

		Eigen::VectorXd v = A_new * X - w_head;

		//卫星总数是否大于未知参数数量，如果卫星数不足，直接返回定位失bai
		if (nObs <= n_params) return false;

		double sigma0 = sqrt((v.transpose() * P * v)(0) / (nObs - n_params));

		//检查𝑥ො是否收敛，如果没有收敛，将初始位置设定为
		double pos_norm = X.head(3).norm();
		double clk_norm = 0.0;
		int idx = 3;
		if (ngps > 0) clk_norm = (std::max)(clk_norm, std::abs(X(idx++)));
		if (nbds > 0) clk_norm = (std::max)(clk_norm, std::abs(X(idx)));

		xr(0) += X(0);
		xr(1) += X(1);
		xr(2) += X(2);
		idx = 3;
		if (ngps > 0) { bGPS += X(idx++); }
		if (nbds > 0) { bBDS += X(idx); }

		if (pos_norm < POS_EPS && clk_norm < CLK_EPS)
		{
			//检查𝑥ො是否收敛，如果没有收敛，将初始位置设定为
			Eigen::MatrixXd Q = N_new.inverse();
			double PDOP = sqrt(Q(0, 0) + Q(1, 1) + Q(2, 2));

			Result->Position[0] = xr(0);
			Result->Position[1] = xr(1);
			Result->Position[2] = xr(2);
			Result->PDOP = PDOP;
			Result->RcvClkOft[0] = bGPS / C_Light;
			Result->RcvClkOft[1] = bBDS / C_Light;
			Result->SigmaPos = sigma0 * sqrt(Q(0,0) + Q(1,1) + Q(2,2));
			Result->IsSuccess = true;
			return true;
		}
	}
	//未收敛
	Result->IsSuccess = false;
	return false;
}

//单点测速
void SPV(EPOCHOBS* Epoch, PPRESULT* Result)
{
	if (Epoch == nullptr || Result == nullptr) return;
	if (!Result->IsSuccess) return;

	Result->Velocity[0] = Result->Velocity[1] = Result->Velocity[2] = 0.0;
	Result->RcvClkSft = 0.0;
	Result->SigmaVel = 999.9;

	Eigen::MatrixXd B(MAXCHANNUM, 4);
	Eigen::VectorXd w(MAXCHANNUM);

	B.setZero();
	w.setZero();

	int nObs = 0;

	for (int i = 0; i < Epoch->SatNum && i < MAXCHANNUM; ++i)
	{
		SATOBS& sat = Epoch->SatObs[i];
		SATMIDRES& mid = Epoch->SatPVT[i];

		if (sat.System != GPS && sat.System != BDS) continue;
		if (sat.Prn <= 0) continue;
		if (!mid.Valid) continue;

		double Dobs = sat.D[0];
		if (fabs(Dobs) < 1e-6) continue;

		XYZ xr_xyz, xs_xyz;
		xr_xyz.x = Result->Position[0];
		xr_xyz.y = Result->Position[1];
		xr_xyz.z = Result->Position[2];

		xs_xyz.x = mid.SatPos[0];
		xs_xyz.y = mid.SatPos[1];
		xs_xyz.z = mid.SatPos[2];

		double elev = 0.0, az = 0.0;
		if (!CompSatEIAz(&xr_xyz, &xs_xyz, &elev, &az)) continue;
		if (elev < 10.0 * Rad) continue;

		double dx = Result->Position[0] - mid.SatPos[0];
		double dy = Result->Position[1] - mid.SatPos[1];
		double dz = Result->Position[2] - mid.SatPos[2];
		double rho = sqrt(dx * dx + dy * dy + dz * dz);
		if (rho < 1.0) continue;

		double l = dx / rho;
		double m = dy / rho;
		double n = dz / rho;

		double satProj = l * mid.SatVel[0] + m * mid.SatVel[1] + n * mid.SatVel[2];

		double satClkRate = C_Light * mid.SatClkSft;

		double wj = Dobs + satProj + satClkRate;

		B(nObs, 0) = l;
		B(nObs, 1) = m;
		B(nObs, 2) = n;
		B(nObs, 3) = 1.0;
		w(nObs) = wj;

		++nObs;
	}

	if (nObs < 4) return;

	Eigen::MatrixXd A = B.topRows(nObs);
	Eigen::VectorXd L = w.head(nObs);
	Eigen::MatrixXd P = Eigen::MatrixXd::Identity(nObs, nObs);

	Eigen::MatrixXd N = A.transpose() * P * A;
	Eigen::VectorXd U = A.transpose() * P * L;
	Eigen::VectorXd Xhat = N.ldlt().solve(U);

	Result->Velocity[0] = Xhat(0);
	Result->Velocity[1] = Xhat(1);
	Result->Velocity[2] = Xhat(2);
	Result->RcvClkSft = Xhat(3) / C_Light;

	if (nObs > 4)
	{
		Eigen::VectorXd v = A * Xhat - L;
		Eigen::MatrixXd Qxx = N.inverse();
		double sigma0 = sqrt((v.transpose() * P * v)(0) / (nObs - 4));
		Result->SigmaVel = sigma0 * sqrt(Qxx(0, 0) + Qxx(1, 1) + Qxx(2, 2));
	}
}