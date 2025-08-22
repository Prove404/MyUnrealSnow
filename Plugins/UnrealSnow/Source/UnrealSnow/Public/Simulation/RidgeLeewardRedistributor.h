#pragma once

#include "CoreMinimal.h"
#include "Simulation/ISnowRedistributor.h"
#include "UObject/Object.h"
#include "RidgeLeewardRedistributor.generated.h"

/** Mass-conserving ridge/leeward redistributor using wind and slope/aspect. */
UCLASS(BlueprintType, EditInlineNew, DefaultToInstanced)
class UNREALSNOW_API URidgeLeewardRedistributor : public UObject, public ISnowRedistributor
{
	GENERATED_BODY()

public:
	// Fraction of cell SWE that can be transported off windward faces [0..1]
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Redistribution")
	float WindwardLossFactor = 0.2f;

	// Shapes downwind spread (0 = broad, 1 = focused). Mass is conserved regardless.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Redistribution")
	float LeewardGainFactor = 0.5f;

	// Slope range [deg] used to activate transport
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Redistribution")
	float MinSlopeDeg = 5.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Redistribution")
	float MaxSlopeDeg = 45.0f;

public:
	// ISnowRedistributor
	virtual void Apply_Implementation(FIntPoint GridSize,
		TArray<float>& InOut_SWE,
		const TArray<float>& Slope_deg,
		const TArray<float>& Aspect_deg,
		FVector2D WindDirDeg,
		float WindSpeed) override;

#if WITH_EDITOR
	UFUNCTION(CallInEditor, Category="UnrealSnow|Tests")
	void RunConservationTest();
#endif
};


