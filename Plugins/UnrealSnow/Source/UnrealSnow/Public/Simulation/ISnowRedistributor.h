#pragma once
#include "CoreMinimal.h"
#include "UObject/Interface.h"
#include "ISnowRedistributor.generated.h"

UINTERFACE(BlueprintType)
class USnowRedistributor : public UInterface
{
    GENERATED_BODY()
};

class ISnowRedistributor
{
    GENERATED_BODY()

public:
    UFUNCTION(BlueprintNativeEvent, BlueprintCallable, Category="UnrealSnow|Redistribution")
    void Apply(FIntPoint GridSize,
               UPARAM(ref) TArray<float>& InOut_SWE,
               const TArray<float>& Slope_deg,
               const TArray<float>& Aspect_deg,
               FVector2D WindDirDeg,
               float WindSpeed);
};
