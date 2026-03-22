#include"RTK_Structs.h"

//通用时-》简化儒略日
bool CommonTimeToMJDTime(const COMMONTIME* CT, MJDTIME* MJD)
{
	if (CT == nullptr || MJD == nullptr || CT->Second < 0)
	{
		return false;
	}

	int y, m;

	if (CT->Month <= 2)
	{
		y = CT->Year - 1;
		m = CT->Month + 12;
	}
	else
	{
		y = CT->Year;
		m = CT->Month;
	}

	double ut = CT->Hour + CT->Minute / 60.0 + CT->Second / 3600.0;

	MJD->Days = floor(365.25 * y) + floor(30.6001 * (m + 1)) + CT->Day - 679019;
	MJD->FracDay = ut / 24.0;

	return true;
}

//简化儒略日-》通用时
bool MJDTimeToCommonTime(const MJDTIME* MJD, COMMONTIME* CT)
{
	if (MJD == nullptr || CT == nullptr || MJD->FracDay < 0.0 || MJD->FracDay >= 1.0)
	{
		return false;
	}
	int a, b, c, d, e;

	double jd = MJD->Days + MJD->FracDay + 2400000.5;
	a = (int)floor(jd + 0.5);
	b = a + 1537;
	c = floor((b - 122.1) / 365.25);
	d = floor(365.25 * c);
	e = floor((b - d) / 30.6001);
	CT->Day = b - d - floor(30.6001 * e);
	CT->Month = e - 1 - 12 * floor(e / 14);
	CT->Year = c - 4715 - floor((7 + CT->Month) / 10);
	CT->Hour = floor(MJD->FracDay * 24);
	CT->Minute = floor((MJD->FracDay * 24 - CT->Hour) * 60);
	CT->Second = ((MJD->FracDay * 24 - CT->Hour) * 60 - CT->Minute) * 60;

	return true;
}

//简化儒略日-》GPS时
bool MJDTimeToGPSTime(const MJDTIME* MJD, GPSTIME* GT)
{
	if (MJD == nullptr || GT == nullptr) return false;
	if (MJD->FracDay < 0.0 || MJD->FracDay >= 1.0) return false;

	GT->Week = floor((MJD->Days + MJD->FracDay - 44244) / 7);
	GT->SecOfWeek = (MJD->Days + MJD->FracDay - 44244 - GT->Week * 7) * 86400;

	//容错处理
	if (GT->SecOfWeek >= 604800.0) 
	{ 
		GT->SecOfWeek -= 604800.0; GT->Week += 1; 
	}

	if (GT->SecOfWeek < 0.0) 
	{ 
		GT->SecOfWeek += 604800.0; GT->Week -= 1; 
	}

	return true;
}

//GPS时-》简化儒略日——————待优化
bool GPSTimeToMJDTime(const GPSTIME* GT, MJDTIME* MJD)
{
	if (MJD == nullptr || GT == nullptr) return false;

	MJD->Days = 44244 + GT->Week * 7 + floor(GT->SecOfWeek / 86400);
	MJD->FracDay = GT->SecOfWeek / 86400 - floor(GT->SecOfWeek / 86400);

	return true;
}

//通用时-》GPS时
bool CommonTimeToGPSTime(const COMMONTIME* CT, GPSTIME* GT)
{
	if (CT == nullptr || GT == nullptr) return false;

	MJDTIME MJD;

	if (!CommonTimeToMJDTime(CT, &MJD)) return false;
	if (!MJDTimeToGPSTime(&MJD, GT)) return false;

	return true;
}

//GPS时-》通用时
bool GPSTimeToCommonTime(const GPSTIME* GT, COMMONTIME* CT)
{
	if (CT == nullptr || GT == nullptr) return false;

	MJDTIME MJD;

	if (!GPSTimeToMJDTime(GT, &MJD)) return false;
	if (!MJDTimeToCommonTime(&MJD, CT)) return false;

	return true;
}