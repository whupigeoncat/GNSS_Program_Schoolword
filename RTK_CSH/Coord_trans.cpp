#include"RTK_Structs.h"


//大地转空间直角
bool BLHToXYZ(const BLH *blh, XYZ *xyz, const double A, const double F)
{
	double e2, N;
	double B = blh->b;
	double L = blh->l;
	double H = blh->h;

	e2 = 2 * F - (F * F);
	N = A / sqrt(1 - (e2 * (sin(B) * sin(B))));

	if(fabs(B)<=PAI/2&&fabs(L)<=PAI)	//经纬度限制
	{
		xyz->x = (N + H) * cos(B) * cos(L);
		xyz->y = (N + H) * cos(B) * sin(L);
		xyz->z = (N * (1 - e2) + H) * sin(B);
		return true;
	}
	return false;
}

//直角转大地
bool XYZToBLH(const XYZ* xyz, BLH* blh, const double A, const double F)
{
	double e2, N, delta_Z[2];
	double X = xyz->x;
	double Y = xyz->y;
	double Z = xyz->z;
	double p = sqrt(X * X + Y * Y);

	e2 = 2 * F - (F * F);
	delta_Z[0] = e2 * Z;

	if (X * X + Y * Y < 1e-10)
	{
		blh->b = 0;
		blh->l = 0;
		blh->h = -6371.0 / 2.0;
		return true;
	}

	blh->b = atan2(Z, (X * X + Y * Y) * (1 - e2));
	blh->l = atan2(Y, X);

	for (int i = 0; i < 10; i++)
	{
		N = A / sqrt(1 - (e2 * (sin(blh->b) * sin(blh->b))));
		blh->b = atan2(Z + delta_Z[0], p);
		blh->h = sqrt(p + (Z + delta_Z[0]) * (Z + delta_Z[0])) - N;
		delta_Z[1] = N * e2 * sin(blh->b);
		if (fabs(delta_Z[0] - delta_Z[1]) < 1e-10) break;
		delta_Z[0] = delta_Z[1];
	}
	if (fabs(delta_Z[0] - delta_Z[1]) > 1e-10)
	{
		return false;
	}
	else
	{
		N = A / sqrt(1 - (e2 * (sin(blh->b) * sin(blh->b))));
		blh->b = atan2(Z + delta_Z[0], sqrt(X * X + Y * Y));
		blh->h = sqrt(X * X + Y * Y + (Z + delta_Z[1]) * (Z + delta_Z[1])) - N;
		return true;
	}
}

//测站地平坐标转换矩阵R
bool BLHToNEUMat(const BLH* blh,Eigen::Matrix3d* R)
{
	if (blh == nullptr || R == nullptr) return false;

	double B = blh->b;
	double L = blh->l;

	*R << -sin(L),			 cos(L),		  0, 
		  -sin(B) * cos(L), -sin(B) * sin(L), cos(B),
		   cos(B) * cos(L),	 cos(B) * sin(L), sin(B);

	return true;
}

//定位误差计算
bool CompEnudPos(const XYZ* xyz0, const XYZ* xyzs, NEU* dneu)
{
	if (xyz0 == nullptr || xyzs == nullptr) return false;

	BLH blh0;
	double X0 = xyz0->x;
	double Y0 = xyz0->y;
	double Z0 = xyz0->z;
	double Xs = xyzs->x;
	double Ys = xyzs->y;
	double Zs = xyzs->z;

	Eigen::Matrix3d R;
	Eigen::Vector3d dEnuVec;
	Eigen::Vector3d dxyzVec;

	if (!XYZToBLH(xyz0, &blh0, R_WGS84, F_WGS84)) return false;
	if (!BLHToNEUMat(&blh0, &R)) return false;
	dxyzVec << Xs - X0, Ys - Y0, Zs - Z0;
	dEnuVec = R * dxyzVec;

	dneu->dE = dEnuVec(0);
	dneu->dN = dEnuVec(1);
	dneu->dU = dEnuVec(2);

	return true;
}

//卫星高度角方位角计算
bool CompSatEIAz(const XYZ* xyz0, const XYZ* xyzs, double* Elev, double* Azim)
{
	if (xyz0 == nullptr || xyzs == nullptr || Elev == nullptr || Azim == nullptr) return false;

	NEU dneu;

	if (!CompEnudPos(xyz0, xyzs, &dneu)) return false;

	*Elev = atan2(dneu.dU, sqrt(dneu.dN * dneu.dN + dneu.dE * dneu.dE));
	*Azim = atan2(dneu.dE, dneu.dN);
	if (*Azim < 0) *Azim += 2 * PAI;

	return true;
}