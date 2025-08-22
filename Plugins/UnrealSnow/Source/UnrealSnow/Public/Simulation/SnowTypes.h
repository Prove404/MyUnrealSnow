#pragma once
#include "CoreMinimal.h"
#include "SnowTypes.generated.h"

USTRUCT(BlueprintType)
struct FSnowState
{
	GENERATED_BODY()

public:
	// Depth in meters
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="Snow")
	float SnowDepth_m = 0.0f;

	// Water equivalent in millimeters
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="Snow")
	float SnowWaterEquivalent_mm = 0.0f;

	// Bulk density kg m^-3
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="Snow")
	float Density_kg_m3 = 100.0f;

	// (Optional) age in days to drive albedo decay
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="Snow")
	float Age_days = 0.0f;
};

USTRUCT(BlueprintType)
struct FTerrainSample
{
	GENERATED_BODY()

public:
	// Degrees [0..90]
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="Terrain")
	float Slope_deg = 0.0f;

	// Aspect degrees [0..360), 0=N, 90=E
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="Terrain")
	float Aspect_deg = 0.0f;
};
