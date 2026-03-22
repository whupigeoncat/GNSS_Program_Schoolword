#include "RTK_Structs.h"
#include <cstdio>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdarg>

#define POLYCRC32   0xEDB88320u /* CRC32 polynomial */

unsigned int crc32(const unsigned char* buff, int len)
{
    int i, j;
    unsigned int crc = 0;

    for (i = 0; i < len; i++)
    {
        crc ^= buff[i];
        for (j = 0; j < 8; j++)
        {
            if (crc & 1) crc = (crc >> 1) ^ POLYCRC32;
            else crc >>= 1;
        }
    }
    return crc;
}

double  R8(unsigned char* p) { double  r; memcpy(&r, p, 8); return r; }
float   R4(unsigned char* p) { float  r; memcpy(&r, p, 4); return r; }
int     I4(unsigned char* p) { int  r; memcpy(&r, p, 4); return r; }
unsigned int UI4(unsigned char* p) { unsigned int r; memcpy(&r, p, 4); return r; }
short   I2(unsigned char* p) { short  r; memcpy(&r, p, 2); return r; }
unsigned short UI2(unsigned char* p) { unsigned short r; memcpy(&r, p, 2); return r; }

bool DataReceive(const std::string& filepath)
{
    FILE* fp = nullptr;
    if (fopen_s(&fp, filepath.c_str(), "rb") != 0 || fp == nullptr)
    {
        return false;
    }

    unsigned char ioBuf[MAXRAWLEN];
    unsigned char streamBuf[MAXRAWLEN * 4];
    int streamLen = 0;

    EPOCHOBS obs{};
    GPSEPHREC geph[MAXGPSNUM]{};
    GPSEPHREC beph[MAXBDSNUM]{};

    bool ok = true;

    while (true)
    {
        const size_t n = fread(ioBuf, 1, sizeof(ioBuf), fp);
        if (n == 0) break;

        if (streamLen + static_cast<int>(n) > static_cast<int>(sizeof(streamBuf)))
        {
            ok = false;
            break;
        }

        memcpy(streamBuf + streamLen, ioBuf, n);
        streamLen += static_cast<int>(n);

        DecodeNovOem7Dat(streamBuf, streamLen, &obs, geph, beph);
    }

    if (ferror(fp)) ok = false;
    fclose(fp);
    return ok;
}

int DecodeNovOem7Dat(unsigned char Buff[], int& Len, EPOCHOBS* obs, GPSEPHREC geph[], GPSEPHREC beph[])
{
    int i = 0;
    int MsgId, MsgLen;
    unsigned char Msg[10240], MsgHead[28];

    while (i < Len)
    {
        // 找同步头
        while (i <= Len - 3 &&
            !(Buff[i] == 0xAA && Buff[i + 1] == 0x44 && Buff[i + 2] == 0x12))
        {
            ++i;
        }

        if (i + 28 > Len) break;

        memcpy(MsgHead, Buff + i, 28);
        MsgId = UI2(MsgHead + 4);
        MsgLen = UI2(MsgHead + 8);

        const int frameLen = 28 + MsgLen + 4;
        if (i + frameLen > Len) break;
        if (frameLen > (int)sizeof(Msg)) { ++i; continue; }

        memcpy(Msg, Buff + i, frameLen);

        if (crc32(Msg, 28 + MsgLen) != UI4(Msg + 28 + MsgLen))
        {
            ++i;
            continue;
        }

        switch (MsgId)
        {
        case 43:   // RANGEB
            decode_rangeb_oem7(Msg, obs);
            break;
        case 7:    // GPSEPHEM(B)
            decode_gpsephem(Msg, geph);
            break;
        case 1696: // BDSEPHEMERIS(B)
            decode_bdsephem(Msg, beph);
            break;
        default:
            break;
        }

        i += frameLen;
    }

    // 循环结束后：把未处理尾巴前移到 Buff 开头，并更新 Len
    if (i < 0) i = 0;
    if (i > Len) i = Len;

    int remain = Len - i;          // 剩余不足一帧/未处理字节数
    if (remain > 0)
    {
        memmove(Buff, Buff + i, remain);  // 可能重叠，必须用 memmove
    }
    Len = remain;                  // 返回给上层（主循环）继续拼接
    return 1;
}

int decode_rangeb_oem7(unsigned char* buff, EPOCHOBS* obs)
{
    if (buff == nullptr || obs == nullptr) return 0;

    const unsigned short msgLen = UI2(buff + 8);
    if (msgLen < 4) return 0;

    const int obsItemLen = 44;
    int ObsNum = static_cast<int>(UI4(buff + 28));
    const int maxObsByLen = (static_cast<int>(msgLen) - 4) / obsItemLen;
    if (ObsNum < 0) ObsNum = 0;
    if (ObsNum > maxObsByLen) ObsNum = maxObsByLen;

    obs->Time.Week = UI2(buff + 14);
    obs->Time.SecOfWeek = UI4(buff + 16) / 1000.0;

    obs->SatNum = 0;
    memset(obs->SatObs, 0, MAXCHANNUM * sizeof(SATOBS));

    unsigned char* p = buff + 32; // 28-byte header + 4-byte ObsNum

    for (int i = 0; i < ObsNum; ++i, p += obsItemLen)
    {
        const unsigned int ch_tr_status = UI4(p + 40);
        const int phaseLockFlag = (ch_tr_status >> 10) & 0x01;
        const int parityFlag = (ch_tr_status >> 11) & 0x01;
        const int codeLockedFlag = (ch_tr_status >> 12) & 0x01;
        const int satSystem = (ch_tr_status >> 16) & 0x07;
        const int signalType = (ch_tr_status >> 21) & 0x1F;
        const int prn = UI2(p);

        GNSSSys sys = UNKS;
        int freq = -1;
        double wl = 0.0;

        if (satSystem == 0)
        {
            sys = GPS;
            if (signalType == 0) { freq = 0; wl = WL1_GPS; }
            else if (signalType == 9) { freq = 1; wl = WL2_GPS; }
            else continue;
        }
        else if (satSystem == 4)
        {
            sys = BDS;
            if (signalType == 0 || signalType == 4) { freq = 0; wl = WL1_BDS; }
            else if (signalType == 2 || signalType == 6) { freq = 1; wl = WL3_BDS; }
            else continue;
        }
        else
        {
            continue;
        }

        int m = -1;
        for (int j = 0; j < MAXCHANNUM; ++j)
        {
            if (obs->SatObs[j].System == sys && obs->SatObs[j].Prn == prn)
            {
                m = j;
                break;
            }
        }
        if (m < 0)
        {
            for (int j = 0; j < MAXCHANNUM; ++j)
            {
                if (obs->SatObs[j].System == UNKS && obs->SatObs[j].Prn == 0)
                {
                    m = j;
                    break;
                }
            }
        }
        if (m < 0) continue;

        obs->SatObs[m].Prn = static_cast<short>(prn);
        obs->SatObs[m].System = sys;
        obs->SatObs[m].P[freq] = (codeLockedFlag == 1) ? R8(p + 4) : 0.0;
        obs->SatObs[m].L[freq] = -wl * ((phaseLockFlag == 1) ? R8(p + 16) : 0.0);
        obs->SatObs[m].D[freq] = -wl * R4(p + 28);
        obs->SatObs[m].LockTime[freq] = R4(p + 36);
        obs->SatObs[m].cn0[freq] = R4(p + 32);
        obs->SatObs[m].half[freq] = static_cast<unsigned char>(parityFlag);
    }

    for (int k = 0; k < MAXCHANNUM; ++k)
    {
        if (obs->SatObs[k].System != UNKS && obs->SatObs[k].Prn != 0) obs->SatNum++;
    }

    return 1;
}

// BDSEPHEMERIS(B)
int decode_bdsephem(unsigned char* buff, GPSEPHREC* eph)
{
    const unsigned short msgLen = UI2(buff + 8);
    if (msgLen < 196) return 0;
    unsigned char* p = buff + 28;
    const unsigned int prn = UI4(p + 0);
    if (prn == 0 || prn > MAXBDSNUM) return 0;
    GPSEPHREC& e = eph[prn - 1];
    e.Sys = BDS;
    e.PRN = static_cast<unsigned short>(prn);

    const short week = static_cast<short>(UI4(p + 4));
    e.SVAccuracy = R8(p + 8);
    e.SVHealth = static_cast<short>(UI4(p + 16));

    e.TGD1 = R8(p + 20);
    e.TGD2 = R8(p + 28);

    e.IODC = static_cast<double>(UI4(p + 36));
    e.TOC.Week = week;
    e.TOC.SecOfWeek = static_cast<double>(UI4(p + 40));

    e.ClkBias = R8(p + 44);
    e.ClkDrift = R8(p + 52);
    e.ClkDriftRate = R8(p + 60);

    e.IODE = static_cast<double>(UI4(p + 68));
    e.TOE.Week = week;
    e.TOE.SecOfWeek = static_cast<double>(UI4(p + 72));

    e.SqrtA = R8(p + 76);
    e.e = R8(p + 84);
    e.omega = R8(p + 92);
    e.DeltaN = R8(p + 100);
    e.M0 = R8(p + 108);
    e.OMEGA = R8(p + 116);
    e.OMEGADot = R8(p + 124);
    e.i0 = R8(p + 132);
    e.iDot = R8(p + 140);

    e.Cuc = R8(p + 148);
    e.Cus = R8(p + 156);
    e.Crc = R8(p + 164);
    e.Crs = R8(p + 172);
    e.Cic = R8(p + 180);
    e.Cis = R8(p + 188);

    return 1;
}

//GPSEPHEM(B)
int decode_gpsephem(unsigned char* buff, GPSEPHREC* eph)
{
    const unsigned short msgLen = UI2(buff + 8);
    if (msgLen < 224) return 0; // 最高访问到 p+216 且读取8字节

    unsigned char* p = buff + 28;
    const unsigned int prn = UI4(p + 0);
    if (prn == 0 || prn > MAXGPSNUM) return 0;

    eph[prn - 1].Sys = GPS;

    eph[prn - 1].PRN = static_cast<unsigned short>(prn); // 存到对应星历槽位
    eph[prn - 1].SVHealth = static_cast<short>(UI4(p + 12));
    eph[prn - 1].IODE = UI4(p + 16);
    eph[prn - 1].TOE.Week = static_cast<short>(UI4(p + 24));
    eph[prn - 1].TOE.SecOfWeek = R8(p + 32);
    const double A = R8(p + 40);
    eph[prn - 1].SqrtA = (A > 0.0) ? std::sqrt(A) : 0.0;
    eph[prn - 1].DeltaN = R8(p + 48);
    eph[prn - 1].M0 = R8(p + 56);
    eph[prn - 1].e = R8(p + 64);
    eph[prn - 1].omega = R8(p + 72);
    eph[prn - 1].Cuc = R8(p + 80);
    eph[prn - 1].Cus = R8(p + 88);
    eph[prn - 1].Crc = R8(p + 96);
    eph[prn - 1].Crs = R8(p + 104);
    eph[prn - 1].Cic = R8(p + 112);
    eph[prn - 1].Cis = R8(p + 120);
    eph[prn - 1].i0 = R8(p + 128);
    eph[prn - 1].iDot = R8(p + 136);
    eph[prn - 1].OMEGA = R8(p + 144);
    eph[prn - 1].OMEGADot = R8(p + 152);
    eph[prn - 1].IODC = static_cast<double>(UI4(p + 160));
    eph[prn - 1].TOC.Week = eph[prn - 1].TOE.Week;
    eph[prn - 1].TOC.SecOfWeek = R8(p + 164);
    eph[prn - 1].TGD1 = R8(p + 172);
    eph[prn - 1].TGD2 = 0.0;
    eph[prn - 1].ClkBias = R8(p + 180);
    eph[prn - 1].ClkDrift = R8(p + 188);
    eph[prn - 1].ClkDriftRate = R8(p + 196);
    eph[prn - 1].SVAccuracy = R8(p + 216);

    return 1;
}

//还没完

int decode_psrpos(unsigned char* buff, PPRESULT* pos)
{
    if (buff == nullptr || pos == nullptr) return 0;

    const unsigned short MsgLen = UI2(buff + 8);
    if (MsgLen < 72) return 0;

    pos->Time.Week = UI2(buff + 14);
    pos->Time.SecOfWeek = UI4(buff + 16) / 1000.0;

    unsigned char* p = buff + 28;

    unsigned int sol_status = UI4(p + 0);
    pos->IsSuccess = (sol_status == 0) ? true : false;

    BLH blh;
    XYZ xyz;

    blh.b = R8(p + 8) * PAI / 180;
    blh.l = R8(p + 16) * PAI / 180;
    blh.h = R8(p + 24) + R4(p + 32);

    BLHToXYZ(&blh, &xyz, R_WGS84, F_WGS84);

    pos->Position[0] = xyz.x;
    pos->Position[1] = xyz.y;
    pos->Position[2] = xyz.z;



    //float lat_sigma = R4(p + 40);
    //float lon_sigma = R4(p + 44);
    //float hgt_sigma = R4(p + 48);
    //pos->SigmaPos = sqrt(lat_sigma * lat_sigma + lon_sigma * lon_sigma + hgt_sigma * hgt_sigma);
    //pos->AllSatNum = p[64];
    ////这里卫星数不确定
    //pos->GPSSatNum = p[65];
    //pos->BDSSatNum = p[65];

    return 1;
}