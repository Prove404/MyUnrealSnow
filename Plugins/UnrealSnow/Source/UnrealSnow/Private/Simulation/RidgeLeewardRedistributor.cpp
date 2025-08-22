#include "Simulation/RidgeLeewardRedistributor.h"
#include "Simulation/ISnowRedistributor.h"

static inline float Clamp01(float V) { return FMath::Clamp(V, 0.0f, 1.0f); }

void URidgeLeewardRedistributor::Apply_Implementation(FIntPoint GridSize,
	TArray<float>& InOut_SWE,
	const TArray<float>& Slope_deg,
	const TArray<float>& Aspect_deg,
	FVector2D WindDirDeg,
	float WindSpeed)
{
	const int32 W = GridSize.X;
	const int32 H = GridSize.Y;
	if (W <= 0 || H <= 0) return;
	const int32 Num = W * H;
	if (InOut_SWE.Num() != Num || Slope_deg.Num() != Num || Aspect_deg.Num() != Num) return;

	TArray<float> Loss; Loss.Init(0.0f, Num);
	TArray<float> Gain; Gain.Init(0.0f, Num);

	// Precompute wind direction unit vector (from-direction in deg -> to-direction add 180)
	const float WindFromDeg = WindDirDeg.X;
	const float WindToDeg = FMath::Fmod(WindFromDeg + 180.0f + 360.0f, 360.0f);
	const float WindToRad = FMath::DegreesToRadians(WindToDeg);
	const float wx = FMath::Sin(WindToRad);
	const float wy = FMath::Cos(WindToRad);

	const float SMin = MinSlopeDeg;
	const float SMax = MaxSlopeDeg;
	const float SRangeInv = (SMax > SMin) ? 1.0f / (SMax - SMin) : 0.0f;

	// Small 3x3 stencil biased downwind
	const int32 Offs[9][2] = {
		{ 0, 0}, { 1, 0}, {-1, 0}, { 0, 1}, { 0,-1},
		{ 1, 1}, { 1,-1}, {-1, 1}, {-1,-1}
	};
	float Weights[9];

	for (int32 j = 0; j < H; ++j)
	{
		for (int32 i = 0; i < W; ++i)
		{
			const int32 idx = j * W + i;
			const float sdeg = Slope_deg[idx];
			const float adeg = Aspect_deg[idx];
			const float slopeW = Clamp01((sdeg - SMin) * SRangeInv);
			const float ang = FMath::Cos(FMath::DegreesToRadians(adeg - WindFromDeg));
			const float w = Clamp01(slopeW * FMath::Max(0.0f, ang));
			const float swe = InOut_SWE[idx];
			const float l = swe * WindwardLossFactor * w;
			Loss[idx] = l;

			// Distribute loss to neighbors downwind using directional weights
			float wsum = 0.0f;
			for (int32 k = 0; k < 9; ++k)
			{
				const int32 di = Offs[k][0];
				const int32 dj = Offs[k][1];
				const float dot = di * wx + dj * wy;
				const float wk = Clamp01(0.5f + 0.5f * dot);
				Weights[k] = wk;
				wsum += wk;
			}
			if (wsum < KINDA_SMALL_NUMBER) continue;
			const float scale = (l * (1.0f + LeewardGainFactor * 0.0f)) / wsum; // factor reserved for shaping
			for (int32 k = 0; k < 9; ++k)
			{
				const int32 di = Offs[k][0];
				const int32 dj = Offs[k][1];
				const int32 ii = FMath::Clamp(i + di, 0, W - 1);
				const int32 jj = FMath::Clamp(j + dj, 0, H - 1);
				Gain[jj * W + ii] += Weights[k] * scale;
			}
		}
	}

	// Apply: InOut = In - Loss + Gain, clamp
	for (int32 idx = 0; idx < Num; ++idx)
	{
		float v = InOut_SWE[idx] - Loss[idx] + Gain[idx];
		InOut_SWE[idx] = FMath::Max(0.0f, v);
	}
}

#if WITH_EDITOR
void URidgeLeewardRedistributor::RunConservationTest()
{
	const int32 W = 16, H = 16, Num = W * H;
	TArray<float> SWE; SWE.Init(10.0f, Num);
	TArray<float> Slope; Slope.Init(20.0f, Num);
	TArray<float> Aspect; Aspect.Init(0.0f, Num);
	const float SumBefore = Algo::Accumulate(SWE, 0.0f);
	Apply_Implementation(FIntPoint(W, H), SWE, Slope, Aspect, FVector2D(90.0f, 0.0f), 5.0f);
	const float SumAfter = Algo::Accumulate(SWE, 0.0f);
	const float Rel = (SumBefore > 0.0f) ? FMath::Abs(SumAfter - SumBefore) / SumBefore : 0.0f;
	UE_LOG(LogTemp, Display, TEXT("Redistribution conservation rel error = %g"), Rel);
}
#endif


