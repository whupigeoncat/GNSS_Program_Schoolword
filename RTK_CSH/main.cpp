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


int main()
{
	cout << "输入模式\n0=文件读取\n1=实时TCP\n";
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
                    }
                }
            }
        }

        CloseSocket(NetGps);
        return 0;
    }
}