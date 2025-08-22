#include "Weather/StochasticWeatherProvider.h"
#include "ISnowWeatherProvider.h"
#include "UObject/Object.h"
#include "UObject/ObjectMacros.h"
#include "Math/UnrealMathUtility.h"

UStochasticWeatherProvider::UStochasticWeatherProvider()
{
	RNG.Initialize(RandomSeed);
}

// Marsaglia and Tsang method for Gamma(k, theta), k>0, theta>0
float UStochasticWeatherProvider::SampleGamma(float ShapeK, float ScaleTheta) const
{
	const float k = FMath::Max(ShapeK, 1e-4f);
	const float theta = FMath::Max(ScaleTheta, 1e-6f);
	if (k < 1.0f)
	{
		// Boosting technique: Gamma(k, theta) = Gamma(k+1, theta) * U^(1/k)
		const float u = FMath::Clamp(RNG.FRand(), 1e-6f, 1.0f);
		return SampleGamma(k + 1.0f, theta) * FMath::Pow(u, 1.0f / k);
	}

	const float d = k - 1.0f / 3.0f;
	const float c = 1.0f / FMath::Sqrt(9.0f * d);
	while (true)
	{
		float x = SampleNormal(0.0f, 1.0f);
		float v = 1.0f + c * x;
		if (v <= 0.0f) continue;
		v = v * v * v;
		float u = RNG.FRand();
		if (u < 1.0f - 0.0331f * (x * x) * (x * x))
		{
			return theta * d * v;
		}
		if (FMath::Loge(u) < 0.5f * x * x + d * (1.0f - v + FMath::Loge(v)))
		{
			return theta * d * v;
		}
	}
}

float UStochasticWeatherProvider::SampleNormal(float Mean, float StdDev) const
{
	float u1 = FMath::Clamp(RNG.FRand(), 1e-6f, 1.0f);
	float u2 = FMath::Clamp(RNG.FRand(), 1e-6f, 1.0f);
	float z0 = FMath::Sqrt(-2.0f * FMath::Loge(u1)) * FMath::Cos(2.0f * PI * u2);
	return Mean + StdDev * z0;
}

bool UStochasticWeatherProvider::GetForTime_Implementation(const FDateTime& Time, FWeatherSample& Out) const
{
    // NOTE: Do not call Execute_GetForTime on 'this' to avoid recursion.
	// 2-state Markov chain update
	const bool bPrevWet = bIsWet;
	const float pWet = bPrevWet ? P_WetGivenWet : P_WetGivenDry;
	bIsWet = (RNG.FRand() < pWet);

	// Precipitation [mm/h]
	Out.Precip_mmph = bIsWet ? SampleGamma(GammaShape_k, GammaScale_theta) : 0.0f;

	// Wind components [m/s]
	Out.U10_ms = SampleNormal(WindU_Mean, WindU_Std);
	Out.V10_ms = SampleNormal(WindV_Mean, WindV_Std);

	// Temperature [K] seasonal baseline + noise
	int32 Year = Time.GetYear();
	int32 Month = Time.GetMonth();
	int32 Day = Time.GetDay();
	FDateTime StartOfYear(Year, 1, 1);
	int32 DOY = (Time.GetTicks() - StartOfYear.GetTicks()) / ETimespan::TicksPerDay + 1;
	float Phase = 2.0f * PI * (static_cast<float>(DOY) / 365.0f);
	Out.AirTemp_K = TempBase_K + TempAmplitude_K * FMath::Sin(Phase) + SampleNormal(0.0f, TempNoiseStd_K);

	// Other fields
	Out.SWdown_Wm2 = 150.0f; // placeholder
	Out.RH_frac = FMath::Clamp(0.6f + SampleNormal(0.0f, 0.1f), 0.0f, 1.0f);
	Out.Pressure_Pa = 80000.0f + SampleNormal(0.0f, 200.0f);

	return true;
}


