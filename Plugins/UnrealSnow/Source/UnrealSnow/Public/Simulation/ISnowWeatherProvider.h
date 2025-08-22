#pragma once

#include "CoreMinimal.h"
#include "UObject/Interface.h"
#include "ISnowWeatherProvider.generated.h"

USTRUCT(BlueprintType)
struct FWeatherSample
{
	GENERATED_BODY()

	UPROPERTY(EditAnywhere, BlueprintReadWrite) float SWdown_Wm2 = 0.f;
	UPROPERTY(EditAnywhere, BlueprintReadWrite) float AirTemp_K   = 273.15f;
	UPROPERTY(EditAnywhere, BlueprintReadWrite) float RH_frac     = 0.5f;
	UPROPERTY(EditAnywhere, BlueprintReadWrite) float U10_ms      = 0.f;
	UPROPERTY(EditAnywhere, BlueprintReadWrite) float V10_ms      = 0.f;
	UPROPERTY(EditAnywhere, BlueprintReadWrite) float Precip_mmph = 0.f;
	UPROPERTY(EditAnywhere, BlueprintReadWrite) float Pressure_Pa = 101325.f;
};

UINTERFACE(BlueprintType)
class USnowWeatherProvider : public UInterface
{
	GENERATED_BODY()
};

class ISnowWeatherProvider
{
	GENERATED_BODY()
public:
	// BlueprintNativeEvent so BPs or C++ can implement it
	UFUNCTION(BlueprintNativeEvent, BlueprintCallable, Category="UnrealSnow|Weather")
	bool GetForTime(const FDateTime& Time, FWeatherSample& Out) const;
};


