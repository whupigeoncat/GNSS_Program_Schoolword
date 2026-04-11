#include"RTK_Structs.h"

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


void ComputeGPSSatOrbitAtSignalTrans(const EPOCHOBS* Epk, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, double RcvPos[3], SATMIDRES* MidRes)
{
	if (Epk == nullptr || GPSEph == nullptr || BDSEph == nullptr) return;
	

	const SATOBS* sat = nullptr;
	for (int i = 0; i < Epk->SatNum && i < MAXCHANNUM; i++)
	{
		if (Epk->SatObs[i].System == GPS || Epk->SatObs[i].System == BDS)
		{
			if (Epk->SatObs[i].P[0] > 1e-6 || Epk->SatObs[i].P[1] > 1e-6)
			{
				sat = &Epk->SatObs[i];
				break;
			}
		}
	}
	if (sat == nullptr) return;

	double P = 0.0;
	if (sat->P[0] > 1e-6)
	{
		P = sat->P[0];
	}
	else if (sat->P[1] > 1e-6)
	{
		P = sat->P[1];
	}
	else return;

	double dtj = 0.0;
	SATMIDRES tmp{};
	for (int k = 0; k < 2; k++)
	{
		GPSTIME ts = Epk->Time;
		ts.SecOfWeek -= (P / C_Light + dtj);
		NormalizedGPSTime(&ts);

		bool ok = false;
		if (sat->System == GPS)
		{
			ok = CompGPSSatPVT(sat->Prn, &ts, &GPSEph[sat->Prn - 1], &tmp);
		}
		else
		{
			ok = CompBDSSatPVT(sat->Prn, &ts, &BDSEph[sat->Prn - 1], &tmp);
		}

		if (!ok) return;

		dtj = tmp.SatClkOft;
		*MidRes = tmp;
	}
	MidRes->Valid = true;
}