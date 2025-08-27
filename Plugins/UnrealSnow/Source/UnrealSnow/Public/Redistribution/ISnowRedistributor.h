#pragma once

#include "CoreMinimal.h"

// Keep this header-only compatible and lightweight. No UObject/UINTERFACE involvement.

// Place the new API in a nested namespace to avoid clashes with the legacy
// Simulation/ISnowRedistributor UObject interface that may coexist in this module.
namespace UnrealSnow
{
namespace Redistribution
{
	struct UNREALSNOW_API FSnowGrid
	{
		int32 W = 0;
		int32 H = 0;
		float CellSizeM = 1.0f;          // meters per texel (assume square cells)
		TArray<float> GroundH_M;         // terrain height [m], size W*H
		TArray<float> SnowDepth_M;       // snow depth [m], size W*H
	};

	struct UNREALSNOW_API FRedistributionParams
	{
		FVector2D WindDir = FVector2D(1.0f, 0.0f); // unit vector in XY
		float WindSpeedMS = 0.0f;                   // m/s
		float DtHours = 0.0f;                       // delta time in hours for this step
		int32 Iterations = 1;                        // small # iterations per tick (1–8)
		float AngleOfReposeDeg = 37.0f;             // ~37°
		float Diffusivity = 0.05f;                  // slope creep coefficient (m^2/h)
		float CapacityK = 0.005f;                   // wind transport scale
		float CapacityP = 2.0f;                     // exponent ~2–3
		float UstarThresh = 6.0f;                   // threshold m/s for transport
		float MaxErodeDepositRate = 0.02f;          // clamp per-iteration (m)
	};

	class UNREALSNOW_API ISnowRedistributor
	{
	public:
		virtual ~ISnowRedistributor() {}
		virtual void Step(FSnowGrid& Grid, const FRedistributionParams& P) = 0;
	};

} // namespace Redistribution
} // namespace UnrealSnow


