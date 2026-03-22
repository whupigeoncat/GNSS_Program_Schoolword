#include"RTK_Structs.h"


bool CompSatClkOff(const int Prn, const GNSSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATMIDRES* Mid)
{
	return true;
}

bool CompGPSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATMIDRES* Mid)
{
	//星历过期判断，星历健康标记Health是否为0，如果过期或者不健康，返回false，否则返回true
	if (Eph->SVHealth == 1)
	{
		Mid->Valid = false;
		return false;
	}

	double n0, tk, n, Mk, Ek, vk, FIk, delta_uk, delta_rk, delta_ik, uk, rk, ik, xk_, yk_, OMEGAk, F, delta_tr, tkc;

	double A = Eph->SqrtA * Eph->SqrtA;
	n0 = sqrt(GM_GPS / (A * A * A));
	tk = GetDiffTime(t, &Eph->TOE);
	n = n0 + Eph->DeltaN;
	Mk = Eph->M0 + n * tk;
	Ek = Mk;

	for (int i = 0; i < 100; i++)
	{
		const double f = Ek - Eph->e * sin(Ek) - Mk;
		const double fp = 1.0 - Eph->e * cos(Ek);
		const double dE = -f / fp;
		Ek += dE;
		if (fabs(dE) < 1e-12) break;
	}

	vk = atan2((sqrt(1 - Eph->e * Eph->e) * sin(Ek)) / (1 - Eph->e * cos(Ek)), (cos(Ek) - Eph->e) / (1 - Eph->e * cos(Ek)));
	FIk = vk + Eph->omega;
	delta_uk = Eph->Cus * sin(2 * FIk) + Eph->Cuc * cos(2 * FIk);
	delta_rk = Eph->Crs * sin(2 * FIk) + Eph->Crc * cos(2 * FIk);
	delta_ik = Eph->Cis * sin(2 * FIk) + Eph->Cic * cos(2 * FIk);
	uk = FIk + delta_uk;
	rk = A * (1 - Eph->e * cos(Ek)) + delta_rk;
	ik = Eph->i0 + delta_ik + Eph->iDot * tk;
	xk_ = rk * cos(uk);
	yk_ = rk * sin(uk);
	OMEGAk = Eph->OMEGA + (Eph->OMEGADot - Omega_WGS) * tk - Omega_WGS * Eph->TOE.SecOfWeek;
	Mid->SatPos[0] = xk_ * cos(OMEGAk) - yk_ * cos(ik) * sin(OMEGAk);
	Mid->SatPos[1] = xk_ * sin(OMEGAk) + yk_ * cos(ik) * cos(OMEGAk);
	Mid->SatPos[2] = yk_ * sin(ik);

	F = -(2 * sqrt(GM_GPS)) / (C_Light * C_Light);
	delta_tr = F * Eph->e * sqrt(A) * sin(Ek);
	tkc = GetDiffTime(t, &Eph->TOC);
	Mid->SatClkOft = Eph->ClkBias + Eph->ClkDrift * tkc + Eph->ClkDriftRate * tkc * tkc + delta_tr;
	Mid->SatClkSft = Eph->ClkDrift + 2.0 * Eph->ClkDriftRate * tkc;
	Mid->Tgd1 = Eph->TGD1;
	Mid->Tgd2 = Eph->TGD2;

	double Ekdot, FIkdot, ukdot, rkdot, Ikdot, OMEGAkdot;
	Ekdot = n / (1 - Eph->e * cos(Ek));
	FIkdot = (sqrt(1 - Eph->e * Eph->e) / (1 - Eph->e * cos(Ek))) * Ekdot;


	//注意BDS星历的参考时间为BDT，GPST - BDT = 14s

	//计算成功，Mid->Valid = true; 否则赋值为false
}