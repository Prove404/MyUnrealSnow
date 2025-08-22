#pragma once

#include "CoreMinimal.h"
#include "Simulation/ISnowWeatherProvider.h"
#include "UObject/Object.h"
#include "Engine/EngineTypes.h"
#include "Era5WeatherProvider.generated.h"

/** ERA5 weather provider (fixtures-backed). */
UCLASS(BlueprintType, EditInlineNew, DefaultToInstanced)
class UNREALSNOW_API UEra5WeatherProvider : public UObject, public ISnowWeatherProvider
{
	GENERATED_BODY()

public:
	// Paths to NetCDF files (optional; if not available, uses RawFixturesDir with simple .raw tiles)
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5")
	TArray<FFilePath> NetCDFPaths;

	// Optional raw fixtures directory (binary float32 tiles per-hour): <Var>_YYYYMMDDHH.raw
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5")
	FDirectoryPath RawFixturesDir;

	// Variable names (ERA5 defaults)
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Variables")
	FString Var_SSRD = TEXT("ssrd");
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Variables")
	FString Var_T2M = TEXT("t2m");
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Variables")
	FString Var_U10 = TEXT("u10");
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Variables")
	FString Var_V10 = TEXT("v10");
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Variables")
	FString Var_TP = TEXT("tp");
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Variables")
	FString Var_SP = TEXT("sp");
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Variables")
	FString Var_HUSS = TEXT("huss");

	// Grid definition for raw fixtures (lat = Lat0 + j * dLat, lon = Lon0 + i * dLon)
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Grid")
	FIntPoint GridSize = FIntPoint(0, 0);
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Grid")
	FVector2D LatLonOrigin = FVector2D(0.0, 0.0); // (Lat0, Lon0)
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Grid")
	FVector2D LatLonStep = FVector2D(0.0, 0.0); // (dLat, dLon)

	// Single sample location (deg). Future: derive from world/landscape
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Sampling")
	FVector2D SampleLatLon = FVector2D(46.0, 7.8);

	// Time range discovered from fixtures
	UPROPERTY(VisibleAnywhere, Category="ERA5")
	mutable FDateTime EarliestTime = FDateTime();
	UPROPERTY(VisibleAnywhere, Category="ERA5")
	mutable FDateTime LatestTime = FDateTime();

	// Unit flags
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Units")
	bool bSSRD_Is_J_per_m2 = true; // divide by 3600 to W/m^2
	UPROPERTY(EditAnywhere, Category="UnrealSnow|ERA5|Units")
	bool bTP_Is_m = true; // convert to mm (per-hour accumulation assumed)

public:
	UEra5WeatherProvider();

	// ISnowWeatherProvider
	virtual bool GetForTime_Implementation(const FDateTime& Time, FWeatherSample& Out) const override;

private:
	mutable bool bInitialized = false;
	mutable TArray<FDateTime> AvailableTimes;

	void EnsureInitialized() const;
	bool ScanRawFixturesTimes(const FString& Dir, const FString& VarName) const;
	static bool ParseTimeFromFilename(const FString& Filename, FDateTime& OutTime);
	bool ReadRawTileAtTime(const FString& Dir, const FString& VarName, const FDateTime& Time, TArray<float>& Out) const;
	float SampleBilinear(const TArray<float>& Tile, double LatDeg, double LonDeg) const;
	bool FindTimeBounds(const FDateTime& T, FDateTime& T0, FDateTime& T1, double& Alpha) const;
};


