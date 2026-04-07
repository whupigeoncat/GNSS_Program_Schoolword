#include"RTK_Structs.h"

bool CompSatClkOff(const int Prn, const GNSSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATMIDRES* Mid)
{
	int Tlimt = 7500.0;
	GPSTIME T_curt = *t;
	GPSEPHREC* Eph;
	if (Sys == GPS) Eph = GPSEph + Prn - 1;
	else if (Sys == BDS)
	{
		Eph = GPSEph + Prn - 1;
		T_curt.SecOfWeek -= 14;
		T_curt.Week -= 1356;
		Tlimt = 3900.0;
	}
	else return false;

	if (Eph->Sys != Sys || Eph->PRN != Prn) return false;

	double tkc = GetDiffTime(t, &Eph->TOC);
	if (fabs(tkc) > Tlimt || Eph->SVHealth != 0) return false;

	Mid->SatClkOft = Eph->ClkBias + Eph->ClkDrift * tkc + Eph->ClkDriftRate * tkc * tkc;
	Mid->SatClkSft = Eph->ClkDrift + 2 * Eph->ClkDriftRate * tkc;
	Mid->Valid = true;

	return true;
}

bool CompGPSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATMIDRES* Mid)
{
	//星历过期判断，星历健康标记Health是否为0，如果过期或者不健康，返回false，否则返回true
	if (Eph->SVHealth != 0)
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
	Mid->Tgd1 = Eph->TGD1;
	Mid->Tgd2 = Eph->TGD2;

	double Ekdot, FIkdot, ukdot, rkdot, Ikdot, OMEGAkdot;
	Ekdot = n / (1 - Eph->e * cos(Ek));
	FIkdot = (sqrt(1 - Eph->e * Eph->e) / (1 - Eph->e * cos(Ek))) * Ekdot;
	ukdot = 2 * (Eph->Cus * cos(2 * FIk) - Eph->Cuc * sin(2 * FIk)) * FIkdot + FIkdot;
	rkdot = A * Eph->e * sin(Ek) * Ekdot + 2 * (Eph->Crs * cos(2 * FIk) - Eph->Crc * sin(2 * FIk)) * FIkdot;
	Ikdot = Eph->iDot + 2 * (Eph->Cis * cos(2 * FIk) - Eph->Cic * sin(2 * FIk)) * FIkdot;
	OMEGAkdot = Eph->OMEGADot - Omega_WGS;

	Eigen::Matrix<double, 3, 4> R_dot;
	R_dot << cos(OMEGAk), -sin(OMEGAk) * cos(ik), -(xk_ * sin(OMEGAk) + yk_ * cos(OMEGAk) * cos(ik)),  yk_ * sin(OMEGAk)* sin(ik),
			 sin(OMEGAk),  cos(OMEGAk) * cos(ik),  (xk_ * cos(OMEGAk) - yk_ * sin(OMEGAk) * cos(ik)), -yk_ * cos(OMEGAk) * sin(ik),
			 0,			   sin(ik),					0,												   yk_ * cos(ik);

	double xk_dot, yk_dot;
	xk_dot = rkdot * cos(uk) - rk * ukdot * sin(uk);
	yk_dot = rkdot * sin(uk) + rk * ukdot * cos(uk);

	Eigen::Vector4d V;
	V << xk_dot, yk_dot, OMEGAkdot, Ikdot;
	Eigen::Vector3d v = R_dot * V;

	Mid->SatVel[0] = v(0);
	Mid->SatVel[1] = v(1);
	Mid->SatVel[2] = v(2);
	
	double delta_trdot;
	delta_trdot = F * Eph->e * sqrt(A) * cos(Ek) * Ekdot;
	Mid->SatClkSft = Eph->ClkDrift + 2.0 * Eph->ClkDriftRate * tkc + delta_trdot;
	Mid->Valid = true;
	return true;

}

bool CompBDSSatPVT(const int Prn, const GPSTIME* t, const GPSEPHREC* Eph, SATMIDRES* Mid)
{

}