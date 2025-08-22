#pragma once

#include "CoreMinimal.h"
#include "Engine/DeveloperSettings.h"
#include "SnowSettings.generated.h"

UCLASS(Config=Game, DefaultConfig, meta=(DisplayName="Unreal Snow Settings"))
class UNREALSNOW_API USnowSettings : public UDeveloperSettings
{
	GENERATED_BODY()

public:
	// Toggle debug features for UnrealSnow systems
	UPROPERTY(EditAnywhere, Config, Category="UnrealSnow")
	bool bEnableDebug = false;

#if WITH_EDITOR
	virtual FName GetCategoryName() const override { return TEXT("UnrealSnow"); }
#endif
};


