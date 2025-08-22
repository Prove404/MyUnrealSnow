#pragma once

#include "CoreMinimal.h"
#include "UObject/Object.h"
#include "UObject/ObjectMacros.h"
#include "ISnowWeatherProvider.h"
#include "StochasticWeatherProvider.generated.h"

/** Simple stochastic weather based on a 2-state (wet/dry) Markov chain. */
UCLASS(BlueprintType, Blueprintable)
class UNREALSNOW_API UStochasticWeatherProvider : public UObject, public ISnowWeatherProvider
{
    GENERATED_BODY()

public:
    // Markov chain transition probabilities
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float P_WetGivenWet = 0.6f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float P_WetGivenDry = 0.2f;

    // Precipitation when wet: Gamma(k, theta) with mean = k*theta [mm/h]
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float GammaShape_k = 2.0f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float GammaScale_theta = 1.0f;

    // Wind components ~ Normal(mu, sigma)
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float WindU_Mean = 2.0f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float WindU_Std = 1.0f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float WindV_Mean = 1.0f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float WindV_Std = 1.0f;

    // Temperature baseline: T(t) = BaseK + AmpK * sin(2pi * (doy/365)) + noise
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float TempBase_K = 268.15f; // -5C

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float TempAmplitude_K = 10.0f;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    float TempNoiseStd_K = 1.5f;

    // Seed
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Weather")
    int32 RandomSeed = 1337;

public:
    UStochasticWeatherProvider();

    // ISnowWeatherProvider
    virtual bool GetForTime_Implementation(const FDateTime& Time, FWeatherSample& Out) const override;

private:
    mutable FRandomStream RNG;
    mutable bool bIsWet = false;

    float SampleGamma(float ShapeK, float ScaleTheta) const;
    float SampleNormal(float Mean, float StdDev) const;
};


