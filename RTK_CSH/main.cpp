#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>
#include <thread>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cctype>
#include "RTK_Structs.h"

using namespace std;

int mode = 0;

// ========== 试验配置 ==========
int    GPS_MODE = 2;      // 0=仅GPS  1=仅BDS  2=GPS+BDS
double ELV_CUTOFF = 0.0;   // 高度截止角（度）
int    FREQ_MODE = 2;      // 1=单频L1/B1  2=双频IF组合
bool   TROPO_ENABLE = true;  // true=加对流层  false=不加
// ==============================


int main()
{
    ofstream out("SPP_Result.txt");
    out << fixed << setprecision(4);

    cout << "输入模式\n0=文件读取(SPP)\n1=实时TCP(SPP)\n2=RTK文件模式\n3=RTK实时模式\n";

	if (!(cin >> mode))
	{
		cin.clear();
		cin.ignore((numeric_limits<streamsize>::max)(), '\n');
		cout << "无效输入";
		return 1;
	}
	if (mode == 0)
	{
        FILE* fp = nullptr;
        unsigned char Buff[MAXRAWLEN];
        int Len, LenRead, val;
        EPOCHOBS Obs;
        GPSEPHREC GpsEph[MAXGPSNUM], BdsEph[MAXBDSNUM];
        RAWDAT raw{};
        PPRESULT result{};

        if (fopen_s(&fp, "C:\\Users\\asus\\Documents\\Program\\Program_GNSS\\RTK_CSH\\Data\\oem719-202603111200.bin", "rb") != 0 || fp == nullptr)
        {
            printf("The file 'oem719-202603111200.bin' was not opened\n");
            return 0;
        }

        Len = 0;
        do {
            LenRead = fread(Buff + Len, sizeof(unsigned char), MAXRAWLEN - Len, fp);
            Len = Len + LenRead;

            val = DecodeNovOem7Dat(Buff, Len, &Obs, GpsEph, BdsEph);

            if (val == 43) // RANGEB观测
            {
                DetectOutlier(&Obs);
                memcpy(raw.GpsEph, GpsEph, sizeof(GpsEph));
                memcpy(raw.BdsEph, BdsEph, sizeof(BdsEph));

                if (SPP(&Obs, &raw, &result))
                {
                    cout << fixed << setprecision(3)
                        <<result.Time.SecOfWeek
                        << "SPP: X=" << result.Position[0]
                        << " Y=" << result.Position[1]
                        << " Z=" << result.Position[2]
                        << " PDOP=" << result.PDOP 
                    //SPV(&Obs, &result);
                    //cout << fixed << setprecision(3)
                    //    << "SPV: Vx=" << result.Velocity[0]
                    //    << " Vy=" << result.Velocity[1]
                    //    << " Vz=" << result.Velocity[2]
                        << " Satnum" << result.AllSatNum
                        << " Bdsnum" << result.BDSSatNum
                        << " Gpsnum" << result.GPSSatNum << "\n";

                    out << result.Time.SecOfWeek << " "
                        << result.Position[0] << " " << result.Position[1] << " " << result.Position[2] << " "
                        << result.GPSSatNum << " " << result.BDSSatNum << " " << result.PDOP << "\n";
                }
            }

        } while (!feof(fp));
        return 1;
	}

    if (mode == 1)
    {
        SOCKET NetGps;
        if (!OpenSocket(NetGps, "47.114.134.129", 7190))
        {
            cout << "tcp error\n";
            return 0;
        }

        unsigned char Buff[MAXRAWLEN]{};
        unsigned char tmp[MAXRAWLEN]{};
        int Len = 0;

        EPOCHOBS Obs{};
        GPSEPHREC GpsEph[MAXGPSNUM]{}, BdsEph[MAXBDSNUM]{};
        RAWDAT raw{};
        PPRESULT result{};

        while (true)
        {
            int lenR = recv(NetGps, (char*)tmp, MAXRAWLEN, 0);
            if (lenR <= 0)
            {
                cout << "tcp recv end\n";
                break;
            }

            if (Len + lenR > MAXRAWLEN) Len = 0; // 防溢出，简单丢弃

            memcpy(Buff + Len, tmp, lenR);
            Len += lenR;

            while (true)
            {
                int val = DecodeNovOem7Dat(Buff, Len, &Obs, GpsEph, BdsEph);
                if (val == 0) break;

                if (val == 43) // RANGEB
                {
                    DetectOutlier(&Obs);
                    memcpy(raw.GpsEph, GpsEph, sizeof(GpsEph));
                    memcpy(raw.BdsEph, BdsEph, sizeof(BdsEph));

                    if (SPP(&Obs, &raw, &result))
                    {
                        cout << fixed << setprecision(3)
                            << result.Time.SecOfWeek
                            << "SPP: X=" << result.Position[0]
                            << " Y=" << result.Position[1]
                            << " Z=" << result.Position[2]
                            << " PDOP=" << result.PDOP
                        //SPV(&Obs, &result);
                        //cout << fixed << setprecision(3)
                            //<< "SPV: Vx=" << result.Velocity[0]
                            //<< " Vy=" << result.Velocity[1]
                            //<< " Vz=" << result.Velocity[2]
                            << " Satnum" << result.AllSatNum
                            << " Bdsnum" << result.BDSSatNum
                            << " Gpsnum" << result.GPSSatNum << "\n";

                        out << result.Time.SecOfWeek << " "
                            << result.Position[0] << " " << result.Position[1] << " " << result.Position[2] << " "
                            << result.GPSSatNum << " " << result.BDSSatNum << " " << result.PDOP << "\n";
                    }
                }
            }
        }

        CloseSocket(NetGps);
        return 0;
    }

    if (mode == 2)
    {
        ROVERCFGINFO cfg;
        char cfgName[256] = "rtk.cfg";
        cout << "配置文件名（回车默认rtk.cfg）:";
        cin.ignore((numeric_limits<streamsize>::max)(), '\n');
        cin.getline(cfgName, 256);
        if (cfgName[0] == '\0') strcpy_s(cfgName, "rtk.cfg");
        if (!ReadRTKConfigInfo(cfgName, cfg))
        {
            return 1;
        }

        if (cfg.IsFileData != 1)
        {
            cout << "配置错误：文件模式需要 FROM FILE OR COM: 1\n";
            return 1;
        }

        // 将配置文件中的高度截止角用于解算，单位为度
        ELV_CUTOFF = cfg.ElevThreshold;

        FILE* FBas = nullptr, * FRov = nullptr;
        if (fopen_s(&FBas, cfg.BasObsDatFile, "rb") != 0 || !FBas || fopen_s(&FRov, cfg.RovObsDatFile, "rb") != 0 || !FRov)
        {
            printf("ERROR: can not open Base/Rover data file\n");
            return 1;
        }

        RAWDAT* Raw = new RAWDAT();
        InitNovDecoder(&Raw->BasDec, 0);
        InitNovDecoder(&Raw->RovDec, 1);
        
        ofstream rtkout("RTK_Sync_Result.txt");
        rtkout << fixed << setprecision(3);
        ofstream rtkpos("RTK_Result.txt");//浮点解结果
        rtkpos << fixed << setprecision(9);
        int nSync = 0, nFail = 0, nRtkOK = 0, nRtkFail = 0;

        

        while (1)
        {
            int s = GetSynGNSSObs(FBas, FRov, Raw);
            if (s < 0) break;
            if (s == 0)
            {
                nFail++;
                continue;
            }
            nSync++;

            // 1. 基准站和流动站形成单差
            FormSDEpochObs(&Raw->BasEpk,&Raw->RovEpk,&Raw->SdObs);

            // 2. 对单差载波进行 MW/GF 周跳探测
            DetectCycleSlip(&Raw->SdObs);
            

            double dt = GetDiffTime(&Raw->RovEpk.Time,&Raw->BasEpk.Time);

            cout << fixed << setprecision(3)
                << "Bas " << Raw->BasEpk.Time.Week << " " << Raw->BasEpk.Time.SecOfWeek
                << " sat=" << Raw->BasEpk.SatNum
                << " | Rov " << Raw->RovEpk.Time.Week << " " << Raw->RovEpk.Time.SecOfWeek
                << " sat=" << Raw->RovEpk.SatNum
                << " | dt=" << dt
                << " | BasePos " << (Raw->BasPosValid ? "OK" : "--");
            if (Raw->BasPosValid)
            {
                cout << " [" << Raw->BasPos.xyz[0] << " " << Raw->BasPos.xyz[1] << " " << Raw->BasPos.xyz[2] << "]";
            }
            cout << "\n";

            rtkout << Raw->RovEpk.Time.Week << " " << Raw->RovEpk.Time.SecOfWeek << " "
                << Raw->BasEpk.Time.SecOfWeek << " " << dt << " "
                << Raw->BasEpk.SatNum << " " << Raw->RovEpk.SatNum << "\n";

            //浮点解
            PPRESULT RovRes;
            if (RTKProcEpoch(Raw, &cfg, &RovRes))
            {
                RTKFLOATRES& fr = Raw->FloatRes;
                bool bFix = Raw->DDObs.bFixed;
                const double* pBest = bFix ? Raw->DDObs.dPos : fr.dPos;
                const double* nBest = bFix ? Raw->DDObs.dNeu : fr.dNeu;
                XYZ rXYZ{ Raw->BasPos.xyz[0] + pBest[0],Raw->BasPos.xyz[1] + pBest[1],Raw->BasPos.xyz[2] + pBest[2] };
                BLH rBLH;
                XYZToBLH(&rXYZ, &rBLH, R_WGS84, F_WGS84);

                nRtkOK++;
                cout << fixed << setprecision(3)
                    << "RTK " << fr.Time.SecOfWeek
                    << ": dN=" << nBest[0] << " dE=" << nBest[1] << " dU=" << nBest[2]
                    << " amb=" << fr.AmbNum << (bFix ? " FIX" : " float")
                    << " ratio=" << Raw->DDObs.Ratio
                    << " rms=" << (bFix ? Raw->DDObs.FixRMS[0] : fr.PosRMS[0])
                    << "/" << (bFix ? Raw->DDObs.FixRMS[1] : fr.PosRMS[1]) << "\n";

                rtkpos << fr.Time.Week << " " << fr.Time.SecOfWeek << " "
                    << rBLH.b * Deg << " " << rBLH.l * Deg << " " << rBLH.h << " "
                    << (bFix ? 1 : 0) << " "
                    << Raw->DDObs.DDSatNum[0] << " " << Raw->DDObs.DDSatNum[1] << " "
                    << fr.AmbNum << " " << fr.Sigma0 << " " << Raw->DDObs.Ratio << "\n";
            }
            else
            {
                nRtkFail++;
            }
        }
        printf("Bas流: RANGEB=%ld 星历=%ld/%ld BESTPOS=%ld CRC错=%ld 丢帧=%ld\n",
            Raw->BasDec.nRange, Raw->BasDec.nGeph, Raw->BasDec.nBeph,
            Raw->BasDec.nPos, Raw->BasDec.nBadCrc, Raw->BasDec.nSyncLoss);
        printf("Rov流: RANGEB=%ld 星历=%ld/%ld BESTPOS=%ld CRC错=%ld 丢帧=%ld\n",
            Raw->RovDec.nRange, Raw->RovDec.nGeph, Raw->RovDec.nBeph,
            Raw->RovDec.nPos, Raw->RovDec.nBadCrc, Raw->RovDec.nSyncLoss);
        printf("同步成功历元=%d 不同步轮=%d\n", nSync, nFail);

        printf("RTK浮点解成功=%d 失败=%d\n", nRtkOK, nRtkFail);

        delete Raw;
        fclose(FBas);
        fclose(FRov);
        return 0;
    }

    if (mode == 3)
    {
        ROVERCFGINFO cfg;
        char cfgName[256] = "rtk_rt.cfg";
        cout << "配置文件名(直接回车用rtk_rt.cfg)：";
        cin.ignore((numeric_limits<streamsize>::max)(), '\n');
        cin.getline(cfgName, 256);
        if (cfgName[0] == '\0') strcpy_s(cfgName, "rtk_rt.cfg");
        if (!ReadRTKConfigInfo(cfgName, cfg)) return 1;

        if (cfg.IsFileData != 2)
        {
            cout << "配置错误：实时模式需要 FROM FILE OR COM: 2\n";
            return 1;
        }

        // 将配置文件中的高度截止角用于解算，单位为度
        ELV_CUTOFF = cfg.ElevThreshold;

        ofstream rtkrt("RTK_RT.txt");
        if (!rtkrt)
        {
            cout << "无法创建实时结果文件 RTK_RT.txt\n";
            return 1;
        }
        rtkrt << fixed << setprecision(9);

        SOCKET BasSock, RovSock;
        if (!OpenSocket(BasSock, cfg.BasNetIP, cfg.BasNetPort) ||
            !OpenSocket(RovSock, cfg.RovNetIP, cfg.RovNetPort))
        {
            cout << "tcp error\n";
            return 0;
        }

        RAWDAT* Raw = new RAWDAT();
        InitNovDecoder(&Raw->BasDec, 0);
        InitNovDecoder(&Raw->RovDec, 1);

        while (true)
        {
            int s = GetSynGNSSObs(BasSock, RovSock, Raw);

            if (s < 0)
            {
                cout << "通信中断\n";
                break;
            }

            if (s == 0)
            {
                // 本轮超时或未同步，不能使用旧历元数据
                continue;
            }

            //单差
            FormSDEpochObs(&Raw->BasEpk,&Raw->RovEpk,&Raw->SdObs);

            //周跳探测
            DetectCycleSlip(&Raw->SdObs);

            double dt = GetDiffTime(&Raw->RovEpk.Time,&Raw->BasEpk.Time);

            cout << fixed << setprecision(3)
                << "Bas " << Raw->BasEpk.Time.SecOfWeek << " sat=" << Raw->BasEpk.SatNum
                << " | Rov " << Raw->RovEpk.Time.SecOfWeek << " sat=" << Raw->RovEpk.SatNum
                << " | dt=" << dt
                << " | BasePos " << (Raw->BasPosValid ? "OK" : "--") << "\n";
            
            //3.RTK浮点解
            PPRESULT RovRes;
            if (RTKProcEpoch(Raw, &cfg, &RovRes))
            {
                RTKFLOATRES& fr = Raw->FloatRes;
                DDCOBS& dd = Raw->DDObs;

                bool bFix = dd.bFixed;
                const double* pBest = bFix ? dd.dPos : fr.dPos;
                const double* nBest = bFix ? dd.dNeu : fr.dNeu;

                // 基线加基站坐标，得到当前RTK流动站绝对位置
                XYZ rXYZ{
                    Raw->BasPos.xyz[0] + pBest[0],
                    Raw->BasPos.xyz[1] + pBest[1],
                    Raw->BasPos.xyz[2] + pBest[2]
                };

                cout << fixed << setprecision(3)
                    << "RTK " << fr.Time.Week << " " << fr.Time.SecOfWeek
                    << (bFix ? " FIX" : " FLOAT")
                    << " ratio=" << dd.Ratio
                    << " Ref(G/C)=" << dd.RefPrn[0] << "/" << dd.RefPrn[1]
                    << " DD(G/C)=" << dd.DDSatNum[0] << "/" << dd.DDSatNum[1]
                    << " amb=" << fr.AmbNum << "\n";

                cout << "  Baseline XYZ(m): "
                    << pBest[0] << " " << pBest[1] << " " << pBest[2]
                    << " | NEU(m): "
                    << nBest[0] << " " << nBest[1] << " " << nBest[2] << "\n";

                cout << "  Rover XYZ(m): "
                    << rXYZ.x << " " << rXYZ.y << " " << rXYZ.z << "\n";

                BLH rBLH;
                if (XYZToBLH(&rXYZ, &rBLH, R_WGS84, F_WGS84))
                {
                    cout << "  Rover BLH(deg,deg,m): "
                        << setprecision(9)
                        << rBLH.b * Deg << " " << rBLH.l * Deg << " "
                        << setprecision(3) << rBLH.h << "\n";

                    // 与文件模式保持相同的列顺序，每个成功历元写入一行
                    rtkrt << fr.Time.Week << " " << fr.Time.SecOfWeek << " "
                        << rBLH.b * Deg << " " << rBLH.l * Deg << " " << rBLH.h << " "
                        << (bFix ? 1 : 0) << " "
                        << dd.DDSatNum[0] << " " << dd.DDSatNum[1] << " "
                        << fr.AmbNum << " " << fr.Sigma0 << " " << dd.Ratio << endl;
                    if (!rtkrt)
                    {
                        cout << "写入实时结果文件 RTK_RT.txt 失败\n";
                        break;
                    }
                }

                cout << setprecision(3);
                if (bFix)
                {
                    cout << "  sigmaH=" << dd.FixRMS[0]
                        << " sigmaU=" << dd.FixRMS[1] << " m";
                }
                else
                {
                    cout << "  sigmaXYZ="
                        << fr.PosRMS[0] << "/"
                        << fr.PosRMS[1] << "/"
                        << fr.PosRMS[2] << " m";
                }
                cout << endl;
            }
            else
            {
                cout << "RTK fail: BasPos=" << (Raw->BasPosValid ? "OK" : "--")
                    << " SPP=" << (RovRes.IsSuccess ? "OK" : "--")
                    << " nSD=" << Raw->SdObs.SatNum << "\n";   // 临时调试
            }
        }

        CloseSocket(BasSock);
        CloseSocket(RovSock);
        delete Raw;
        return 0;
    }
}
