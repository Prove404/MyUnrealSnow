#pragma once

#include "CoreMinimal.h"
#include "UObject/Object.h"
#include "Engine/TextureRenderTarget2D.h"
#include "SWROutExporter.generated.h"

UCLASS(BlueprintType, Blueprintable, EditInlineNew)
class UNREALSNOW_API USWROutExporter : public UObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Export")
	UTextureRenderTarget2D* TargetRT = nullptr;

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Export")
	int32 ExportEveryNSteps = 60; // hourly if tick ~1min

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Export")
	FDirectoryPath OutputDir;

	// Geo transform (world-file): pixel size (m), origin (x0,y0), rotation=0
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Export")
	FVector2D PixelSizeMeters = FVector2D(1.0, 1.0);

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Export")
	FVector2D OriginMeters = FVector2D(0.0, 0.0);

public:
	UFUNCTION(CallInEditor, Category="UnrealSnow|Export")
	void ExportNow();

	void TickExport(int32 StepIndex, const FDateTime& SimTime);

private:
	void EnqueueWrite(const FString& BasePath, UTextureRenderTarget2D* RT);
};


