#include<cmath>
#include <cstring>
/****************************************************************************
  MatrixInv
  
  目的：矩阵求逆,采用全选主元高斯-约当法  
  
  编号：01016

  参数:
  n      M1的行数和列数
  a      输入矩阵
  b      输出矩阵   b=inv(a)
  返回值：1=正常，0=致命错误, -1=求逆失败

****************************************************************************/

int MatrixInv(int n, double a[], double b[])
{
    int i,j,k,l,u,v,is[60],js[60];   /* matrix dimension <= 60 */
	double d, p;

	/* 将输入矩阵赋值给输出矩阵b，下面对b矩阵求逆，a矩阵不变 */
	memcpy(b, a, n * n * sizeof(double));
	for (k = 0; k < n; k++)
	{
		d = 0.0;
		for (i = k; i < n; i++) /* 查找右下角方阵中主元素的位置 */
		{
			for (j = k; j < n; j++)
			{
				l = n * i + j;
				p = fabs(b[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (d < 1E-15) /* 主元素接近于0，矩阵不可逆 */
		{
			return (-1);
		}

		if (is[k] != k) /* 对主元素所在的行与右下角方阵的首行进行调换 */
		{
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		if (js[k] != k) /* 对主元素所在的列与右下角方阵的首列进行调换 */
		{
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}

		l = k * n + k;
		b[l] = 1.0 / b[l]; /* 初等行变换 */
		for (j = 0; j < n; j++)
		{
			if (j != k)
			{
				u = k * n + j;
				b[u] = b[u] * b[l];
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				for (j = 0; j < n; j++)
				{
					if (j != k)
					{
						u = i * n + j;
						b[u] = b[u] - b[i * n + k] * b[k * n + j];
					}
				}
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				u = i * n + k;
				b[u] = -b[u] * b[l];
			}
		}
	}

	for (k = n - 1; k >= 0; k--) /* 将上面的行列调换重新恢复 */
	{
		if (js[k] != k)
		{
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				u = i * n + k;
				v = is[k] + i * n;
				p = b[u];
				b[u] = b[v];
				b[v] = p;
			}
		}
	}

	return (1);
}
