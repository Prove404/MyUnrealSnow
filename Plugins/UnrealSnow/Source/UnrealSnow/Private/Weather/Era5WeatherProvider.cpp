#include "Weather/Era5WeatherProvider.h"

#include "Misc/Paths.h"
#include "HAL/FileManager.h"

UEra5WeatherProvider::UEra5WeatherProvider()
{
}

static FString MakeRawFilename(const FString& Var, const FDateTime& T)
{
    return FString::Printf(TEXT("%s_%04d%02d%02d%02d.raw"), *Var, T.GetYear(), T.GetMonth(), T.GetDay(), T.GetHour());
}

void UEra5WeatherProvider::EnsureInitialized() const
{
    if (bInitialized)
    {
        return;
    }
    bInitialized = true;

    if (!RawFixturesDir.Path.IsEmpty())
    {
        const FString Dir = RawFixturesDir.Path;
        // Discover hourly times from one variable (e.g., t2m)
        ScanRawFixturesTimes(Dir, Var_T2M);
        if (AvailableTimes.Num() > 0)
        {
            EarliestTime = AvailableTimes[0];
            LatestTime = AvailableTimes.Last();
        }
    }
}

bool UEra5WeatherProvider::ScanRawFixturesTimes(const FString& Dir, const FString& VarName) const
{
    AvailableTimes.Reset();
    TArray<FString> Files;
    IFileManager::Get().FindFiles(Files, *(FPaths::Combine(Dir, FString::Printf(TEXT("%s_*.raw"), *VarName))), true, false);
    Files.Sort();
    for (const FString& File : Files)
    {
        FDateTime T;
        if (ParseTimeFromFilename(File, T))
        {
            AvailableTimes.Add(T);
        }
    }
    return AvailableTimes.Num() > 0;
}

bool UEra5WeatherProvider::ParseTimeFromFilename(const FString& Filename, FDateTime& OutTime)
{
    // Expect pattern Var_YYYYMMDDHH.raw
    int32 UnderscoreIdx;
    if (!Filename.FindChar(TEXT('_'), UnderscoreIdx)) return false;
    const int32 Start = UnderscoreIdx + 1;
    if (Filename.Len() < Start + 10) return false;
    const FString TS = Filename.Mid(Start, 10);
    const int32 Y = FCString::Atoi(*TS.Mid(0, 4));
    const int32 M = FCString::Atoi(*TS.Mid(4, 2));
    const int32 D = FCString::Atoi(*TS.Mid(6, 2));
    const int32 H = FCString::Atoi(*TS.Mid(8, 2));
    OutTime = FDateTime(Y, M, D, H);
    return true;
}

bool UEra5WeatherProvider::ReadRawTileAtTime(const FString& Dir, const FString& VarName, const FDateTime& Time, TArray<float>& Out) const
{
    if (GridSize.X <= 0 || GridSize.Y <= 0) return false;
    const FString Path = FPaths::Combine(Dir, MakeRawFilename(VarName, Time));
    TArray<uint8> Bytes;
    if (!FFileHelper::LoadFileToArray(Bytes, *Path))
    {
        return false;
    }
    const int64 Expected = static_cast<int64>(GridSize.X) * GridSize.Y * sizeof(float);
    if (Bytes.Num() != Expected)
    {
        return false;
    }
    Out.SetNumUninitialized(GridSize.X * GridSize.Y);
    FMemory::Memcpy(Out.GetData(), Bytes.GetData(), Expected);
    return true;
}

float UEra5WeatherProvider::SampleBilinear(const TArray<float>& Tile, double LatDeg, double LonDeg) const
{
    if (GridSize.X <= 1 || GridSize.Y <= 1) return 0.0f;
    const double fI = (LonDeg - LatLonOrigin.Y) / LatLonStep.Y; // i along X
    const double fJ = (LatDeg - LatLonOrigin.X) / LatLonStep.X; // j along Y
    const int32 i0 = FMath::Clamp((int32)FMath::FloorToInt(fI), 0, GridSize.X - 2);
    const int32 j0 = FMath::Clamp((int32)FMath::FloorToInt(fJ), 0, GridSize.Y - 2);
    const int32 i1 = i0 + 1;
    const int32 j1 = j0 + 1;
    const float tx = (float)FMath::Clamp(fI - i0, 0.0, 1.0);
    const float ty = (float)FMath::Clamp(fJ - j0, 0.0, 1.0);
    auto IDX = [this](int32 I, int32 J) { return J * GridSize.X + I; };
    const float v00 = Tile[IDX(i0, j0)];
    const float v10 = Tile[IDX(i1, j0)];
    const float v01 = Tile[IDX(i0, j1)];
    const float v11 = Tile[IDX(i1, j1)];
    const float v0 = FMath::Lerp(v00, v10, tx);
    const float v1 = FMath::Lerp(v01, v11, tx);
    return FMath::Lerp(v0, v1, ty);
}

bool UEra5WeatherProvider::FindTimeBounds(const FDateTime& T, FDateTime& T0, FDateTime& T1, double& Alpha) const
{
    if (AvailableTimes.Num() == 0) return false;
    if (T < AvailableTimes[0] || T > AvailableTimes.Last()) return false;
    for (int32 k = 0; k < AvailableTimes.Num() - 1; ++k)
    {
        if (T >= AvailableTimes[k] && T <= AvailableTimes[k + 1])
        {
            T0 = AvailableTimes[k];
            T1 = AvailableTimes[k + 1];
            const double dt = (T1 - T0).GetTotalSeconds();
            const double dt0 = (T - T0).GetTotalSeconds();
            Alpha = dt > 0.0 ? dt0 / dt : 0.0;
            return true;
        }
    }
    return false;
}

bool UEra5WeatherProvider::GetForTime_Implementation(const FDateTime& Time, FWeatherSample& Out) const
{
    // NOTE: Do not call Execute_GetForTime on 'this' to avoid recursion.
    EnsureInitialized();
    const FString Dir = RawFixturesDir.Path;
    FDateTime T0, T1;
    double A = 0.0;
    if (!FindTimeBounds(Time, T0, T1, A))
    {
        return false;
    }

    TArray<float> t0_t2m, t1_t2m;
    TArray<float> t0_u10, t1_u10;
    TArray<float> t0_v10, t1_v10;
    TArray<float> t0_tp,  t1_tp;
    TArray<float> t0_sp,  t1_sp;
    TArray<float> t0_ssrd, t1_ssrd;

    // Load required variables
    if (!ReadRawTileAtTime(Dir, Var_T2M, T0, t0_t2m) || !ReadRawTileAtTime(Dir, Var_T2M, T1, t1_t2m)) return false;
    if (!ReadRawTileAtTime(Dir, Var_U10, T0, t0_u10) || !ReadRawTileAtTime(Dir, Var_U10, T1, t1_u10)) return false;
    if (!ReadRawTileAtTime(Dir, Var_V10, T0, t0_v10) || !ReadRawTileAtTime(Dir, Var_V10, T1, t1_v10)) return false;
    if (!ReadRawTileAtTime(Dir, Var_TP,  T0, t0_tp)  || !ReadRawTileAtTime(Dir, Var_TP,  T1, t1_tp))  return false;
    if (!ReadRawTileAtTime(Dir, Var_SP,  T0, t0_sp)  || !ReadRawTileAtTime(Dir, Var_SP,  T1, t1_sp))  return false;
    // Optional SSRD
    const bool bHasSSRD = ReadRawTileAtTime(Dir, Var_SSRD, T0, t0_ssrd) && ReadRawTileAtTime(Dir, Var_SSRD, T1, t1_ssrd);

    const double Lat = SampleLatLon.X;
    const double Lon = SampleLatLon.Y;

    auto LerpAt = [&](const TArray<float>& A0, const TArray<float>& A1) -> float
    {
        const float v0 = SampleBilinear(A0, Lat, Lon);
        const float v1 = SampleBilinear(A1, Lat, Lon);
        return FMath::Lerp(v0, v1, (float)A);
    };

    Out.AirTemp_K = LerpAt(t0_t2m, t1_t2m);
    Out.U10_ms = LerpAt(t0_u10, t1_u10);
    Out.V10_ms = LerpAt(t0_v10, t1_v10);
    const float tp_native = LerpAt(t0_tp, t1_tp);
    const float sp_native = LerpAt(t0_sp, t1_sp);

    // Units
    Out.Pressure_Pa = bTP_Is_m ? sp_native : sp_native; // sp assumed Pa already
    if (bHasSSRD)
    {
        const float ssrd_native = LerpAt(t0_ssrd, t1_ssrd);
        Out.SWdown_Wm2 = bSSRD_Is_J_per_m2 ? (ssrd_native / 3600.0f) : ssrd_native;
    }
    else
    {
        Out.SWdown_Wm2 = 150.0f;
    }

    // ERA5 tp often m of water. Convert to mm/h
    Out.Precip_mmph = bTP_Is_m ? (tp_native * 1000.0f) : tp_native;
    Out.RH_frac = 0.8f; // optional if huss not present

    return true;
}


