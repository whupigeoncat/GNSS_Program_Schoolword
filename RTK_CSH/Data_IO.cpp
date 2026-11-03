#include "RTK_Structs.h"
#include <cstdio>
#include <cstring>


//解析配置文件
bool ReadRTKConfigInfo(const char FName[], ROVERCFGINFO& ANInfo)
{
	int i;
	FILE* Frin;
	char Line[256];
#define CFGITEMS 12
	const char CfgStr[CFGITEMS][32] = {
		"FROM FILE OR COM",
		"BASE IP ADDRESS AND PORT",
		"ROVER IP ADDRESS AND PORT",
		"ROVER COM SETUP",
		"BASE OBSDATA SOURCE FILE",
		"ROVER OBSDATA SOURCE FILE",
		"POSITION RESULT FILE",
		"POSITION DIFF FILE",
		"CODE AND CARRIER PHASE NOISE",
		"THRESHOLD FOR ELEVATION MASK",
		"RATIO FOR DD FIXED SOLUTION",
		"RTK PROCESSING MODE"
	};

	if (fopen_s(&Frin, FName, "rt") != 0 || Frin == NULL)
	{
		printf("ERROR: Can not Open SATODS configuration Info File: %s\n", FName);
		return false;
	}

	while (fgets(Line, sizeof(Line), Frin))
	{
		char* pCom = strchr(Line, '!');      //截断'!'后的行内注释
		if (pCom) *pCom = '\0';

		for (i = 0; i < CFGITEMS; i++)
		{
			if (strncmp(Line, CfgStr[i], strlen(CfgStr[i])) != 0) continue;
			char* p = Line + strlen(CfgStr[i]);//跳过键名，再跳过冒号和空格
			while (*p == ' ' || *p == ':') p++;

			switch (i)
			{
			case 0:  sscanf_s(p, "%hd", &ANInfo.IsFileData); break;
			case 1:  sscanf_s(p, "%19s %hd", ANInfo.BasNetIP, (unsigned)sizeof(ANInfo.BasNetIP), &ANInfo.BasNetPort); break;
			case 2:  sscanf_s(p, "%19s %hd", ANInfo.RovNetIP, (unsigned)sizeof(ANInfo.RovNetIP), &ANInfo.RovNetPort); break;
			case 3:  sscanf_s(p, "%d %d", &ANInfo.RovPort, &ANInfo.RovBaud); break;
			case 4:  sscanf_s(p, "%255s", ANInfo.BasObsDatFile, (unsigned)sizeof(ANInfo.BasObsDatFile)); break;
			case 5:  sscanf_s(p, "%255s", ANInfo.RovObsDatFile, (unsigned)sizeof(ANInfo.RovObsDatFile)); break;
			case 6:  sscanf_s(p, "%255s", ANInfo.ResFile, (unsigned)sizeof(ANInfo.ResFile)); break;
			case 7:  sscanf_s(p, "%255s", ANInfo.DiffFile, (unsigned)sizeof(ANInfo.DiffFile)); break;
			case 8:  sscanf_s(p, "%lf %lf", &ANInfo.CodeNoise, &ANInfo.CPNoise); break;
			case 9:  sscanf_s(p, "%lf", &ANInfo.ElevThreshold); break;
			case 10: sscanf_s(p, "%lf", &ANInfo.RatioThres); break;
			case 11: sscanf_s(p, "%hd", &ANInfo.RTKProcMode); break;
			}
			break;
		}
	}
	fclose(Frin);
	return true;
}

//getdifftime配对判定，thr同步阈值
static int PairEpoch(RAWDAT* Raw, double thr)
{
	if (Raw->BasEpk.Time.Week == 0 && Raw->BasEpk.Time.SecOfWeek == 0.0) return -1;//基站历元未收到，追基站
	double dt = GetDiffTime(&Raw->RovEpk.Time, &Raw->BasEpk.Time);
	if (fabs(dt) <= thr) return 1;//同步成功
	if (dt < 0.0) return 0;//流动站早于基站，丢弃当前流动站历元取下一个
	return -1;//流动站晚于基站，继续读基站
}

//用select等待socket可读：1=有数据可读 0=超时 -1=通信错误
static int WaitReadable(SOCKET sock, int timeoutMs)
{
	fd_set fds;
	FD_ZERO(&fds);
	FD_SET(sock, &fds);

	struct timeval tv;
	tv.tv_sec = timeoutMs / 1000;
	tv.tv_usec = (timeoutMs % 1000) * 1000;

	//Windows下select第一个参数被忽略，超时内任一socket可读即返回>0
	int n = select(0, &fds, NULL, NULL, &tv);
	if (n == SOCKET_ERROR) return -1;
	return (n > 0) ? 1 : 0;
}

//从一条流读一块数据并解析存入RTK数据结构体
static int PumpStream(NOVDECODER* dec, EPOCHOBS* epk, GPSEPHREC geph[], GPSEPHREC beph[], PPRESULT* pos, FILE* fp, SOCKET* pSock)
{
	unsigned char tmp[MAXRAWLEN];
	int n = 0;

	if (dec->Len >= MAXRAWLEN)//缓冲区已满且无完整帧：整体丢弃重新同步，避免0长度读被误判为断连
	{
		dec->Len = 0;
		dec->nSyncLoss++;
	}

	if (fp) n = (int)fread(tmp, 1, MAXRAWLEN - dec->Len, fp);//读文件
	if (pSock) n = recv(*pSock, (char*)tmp, MAXRAWLEN - dec->Len, 0);//读数据流
	if (n <= 0) return -1;
	FeedNovData(dec, tmp, n);//新读到的数据追加进流缓冲区准备解析

	int id;
	while ((id = DecodeNovStream(dec, epk, geph, beph, pos))!= 0)//循环解码直到Msg=43
	{
		if (id == 43) return 1;//解到了rangeb
	}
	return 0;
}

//数据同步，文件模式
int GetSynGNSSObs(FILE* FBas, FILE* FRov, RAWDAT* Raw)
{
	if (!Raw || !FBas || !FRov) return -1;
	const double thr = SYN_THRESH_FILE;
	PPRESULT pos{};//基站BESTPOS临时结果
	for (;;)
	{
		int r = PumpStream(&Raw->RovDec, &Raw->RovEpk, Raw->GpsEph, Raw->BdsEph, nullptr, FRov, nullptr);
		if (r < 0) return -1;
		if (r == 0) continue;
		//与基站历元求差判定匹配
		int p = PairEpoch(Raw, thr);
		if (p == 1) return 1;//同步成功
		if (p == 0) continue;//丢弃该流动站历元，去下一个

		//流动站晚于基站，追基站
		for (;;)
		{
			r = PumpStream(&Raw->BasDec, &Raw->BasEpk, Raw->GpsEph, Raw->BdsEph, &pos, FBas, nullptr);
			if (r < 0) return -1;//基站文件结束，RTK处理结束
			if (r == 0) continue;
			if (pos.IsSuccess && !Raw->BasPosValid)//首个BESTPOS为固定基准站坐标rb
			{
				memcpy(Raw->BasPos.xyz, pos.Position, sizeof(double) * 3);
				Raw->BasPosValid = true;
			}
			p = PairEpoch(Raw, thr);
			if (p == 1) return 1;
			if (p == 0) break;//流动站过旧，同步不成功，回外层取新流动站历元
		}
	}
}

//实时模式，阈值放宽
int GetSynGNSSObs(SOCKET& BasSock, SOCKET& RovSock, RAWDAT* Raw)
{
	if (!Raw || BasSock == INVALID_SOCKET || RovSock == INVALID_SOCKET) return -1;
	const double thr = SYN_THRESH_RT;
	PPRESULT pos{};//基站BESTPOS临时结果
	for (;;)
	{
		int r = 0;
		while ((r = WaitReadable(RovSock, 100)) == 0)
		{
			if (WaitReadable(BasSock, 0) > 0)//基站有数据，先吸收掉
			{
				int rb = PumpStream(&Raw->BasDec, &Raw->BasEpk, Raw->GpsEph, Raw->BdsEph, &pos, nullptr, &BasSock);
				if (rb < 0) return -1;//通信中断
				if (pos.IsSuccess && !Raw->BasPosValid)//首次有效单点定位作为基站坐标
				{
					memcpy(Raw->BasPos.xyz, pos.Position, sizeof(double) * 3);
					Raw->BasPosValid = true;
				}
			}
		}
		if (r < 0) return -1;//流动站错误，中断
		int rr = PumpStream(&Raw->RovDec, &Raw->RovEpk, Raw->GpsEph, Raw->BdsEph, nullptr, nullptr, &RovSock);
		if (rr < 0) return -1;
		if (rr == 0) continue;//历元不完整

		//与基站历元求差判定匹配
		int p = PairEpoch(Raw, thr);
		if (p == 1) return 1;//同步成功
		if (p == 0) continue;//丢弃当前流动站历元，取下一个

		//流动站晚于基站，等基站,即p==-1
		for (;;)
		{
			int wb = WaitReadable(BasSock, 200);
			if (wb < 0) return -1;
			if (wb == 0) return 0;
			r = PumpStream(&Raw->BasDec, &Raw->BasEpk, Raw->GpsEph, Raw->BdsEph, &pos, nullptr, &BasSock);
			if (r < 0) return -1;
			if (r == 0) continue;
			if (pos.IsSuccess && !Raw->BasPosValid)
			{
				memcpy(Raw->BasPos.xyz, pos.Position, sizeof(double) * 3);
				Raw->BasPosValid = true;
			}
			p = PairEpoch(Raw, thr);
			if (p == 1) return 1;
			if (p == 0) break;
		}
	}
}