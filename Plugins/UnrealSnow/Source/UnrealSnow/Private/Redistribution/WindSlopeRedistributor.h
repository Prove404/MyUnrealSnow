#pragma once

#include "CoreMinimal.h"
#include "Redistribution/ISnowRedistributor.h"

namespace UnrealSnow
{
namespace Redistribution
{
	/**
	 * CPU wind + slope-based redistributor. Mass-conserving with 5-point stencil.
	 * Reuses internal buffers across calls to avoid per-tick allocations.
	 */
	class FWindSlopeRedistributor final : public ISnowRedistributor
	{
	public:
		FWindSlopeRedistributor() = default;
		virtual ~FWindSlopeRedistributor() override = default;

		virtual void Step(FSnowGrid& Grid, const FRedistributionParams& P) override;

	private:
		// Temporary working buffers (size W*H)
		TArray<float> WorkA; // general purpose (e.g., next depth)
		TArray<float> WorkB; // general purpose (e.g., flux accum)
		TArray<FVector2f> GradZ; // terrain gradient
		TArray<FVector2f> GradH; // snow gradient

		void EnsureSize(int32 NumCells);
		FORCEINLINE static int32 Idx(int32 x, int32 y, int32 W) { return y * W + x; }
		void ComputeGradients(const FSnowGrid& G);
		void SlopeCreep(FSnowGrid& G, const FRedistributionParams& P, float TanPhi);
		void WindAdvection(FSnowGrid& G, const FRedistributionParams& P, float TanPhi);
	};

} // namespace Redistribution
} // namespace UnrealSnow


