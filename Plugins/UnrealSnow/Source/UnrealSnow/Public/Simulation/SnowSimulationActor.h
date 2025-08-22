#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/StaticMeshComponent.h"
// forward include allowed for clarity
#include "Weather/StochasticWeatherProvider.h"
class UTextureRenderTarget2D;
#include "SnowSimulationActor.generated.h"

UCLASS()
class UNREALSNOW_API ASnowSimulationActor : public AActor
{
	GENERATED_BODY()

protected:
	// Simulation tick timestep [hours]
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Simulation")
	float TimeStepHours = 1.0f;

	// Simulation grid size. If zero, will be inferred from landscape/render target elsewhere.
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Simulation")
	FIntPoint GridSize = FIntPoint(0, 0);

	// Explicit grid controls
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Grid")
	int32 GridSizeXY = 128;  // NxN for now

	UPROPERTY(VisibleAnywhere, Category="UnrealSnow|Grid")
	int32 NumCells = 0;

	// Simulation time scale and start
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Time")
	float HoursPerSecond = 1.0f;  // sim time scale

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Time")
	FDateTime StartTime = FDateTime(2024,1,1,0,0,0);

	// Debug snowfall override
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Debug")
	bool bForceConstantSnow = true;

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Debug", meta=(EditCondition="bForceConstantSnow"))
	float DebugSnow_mmph = 1.0f; // 1 mm/hour

	// Debug visualization material
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Debug")
	UMaterialInterface* DebugMaterial = nullptr;

	UPROPERTY(Transient)
	UMaterialInstanceDynamic* DebugMID = nullptr;

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Debug")
	FName DebugTextureParam = "SnowDepthTex";

	UPROPERTY(VisibleAnywhere, Category="UnrealSnow|Debug")
	UStaticMeshComponent* DebugPlane = nullptr;

	// Pick a Blueprint class for the provider in the Details panel (no inline expansion).
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category="UnrealSnow|Weather")
	TSoftClassPtr<UStochasticWeatherProvider> WeatherProviderClass;

	// Runtime instance created from WeatherProviderClass. Not editable.
	UPROPERTY(Transient)
	UStochasticWeatherProvider* WeatherProvider = nullptr;

	// Redistributor used to apply wind-driven transport effects.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category="UnrealSnow|Redistribution")
	UObject* Redistributor = nullptr;

	// Output snow depth as a render target for visualization.
	UPROPERTY(VisibleAnywhere, Category="UnrealSnow|Outputs")
	UTextureRenderTarget2D* SnowDepthRT = nullptr;

	// Toggle GPU compute path (SnowCS) instead of CPU accumulation.
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Simulation")
	bool bUseGPU = false;

	// SWRout albedo parameters
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Radiation")
	float AlbedoFresh = 0.9f;

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Radiation")
	float AlbedoOld = 0.5f;

	UPROPERTY(EditAnywhere, Category="UnrealSnow|Radiation")
	float AlbedoDecayDays = 7.0f;

	// Output SWRout render target
	UPROPERTY(VisibleAnywhere, Category="UnrealSnow|Outputs")
	UTextureRenderTarget2D* SWROutRT = nullptr;

	// Optional exporter object
	UPROPERTY(EditAnywhere, Category="UnrealSnow|Outputs")
	class USWROutExporter* SWROutExporter = nullptr;

public:
	ASnowSimulationActor();
	virtual void BeginPlay() override;
	virtual void Tick(float DeltaSeconds) override;

	// Editor/UI helpers
	UFUNCTION(BlueprintCallable, Category="UnrealSnow|Simulation")
	void SetSimRunning(bool bRun);

	UFUNCTION(BlueprintCallable, Category="UnrealSnow|Simulation")
	void SimulateStep(float DeltaSeconds);

	UFUNCTION(BlueprintCallable, Category="UnrealSnow|Simulation")
	FString GetCurrentSimTimeISO() const { return CurrentSimTime.ToString(TEXT("yyyy-MM-dd HH:mm:ss")); }

	UFUNCTION(BlueprintCallable, Category="UnrealSnow|Simulation")
	float GetMeanSWE() const { return LastMeanSWE; }

	UFUNCTION(BlueprintCallable, Category="UnrealSnow|Simulation")
	float GetMeanSWRout() const { return LastMeanSWRout; }

private:
	// State fields for a minimal CPU simulation
	TArray<float> SnowWaterEq_mm; // SWE field [mm]
	TArray<float> Slope_deg;      // Optional input for redistribution (zeros for now)
	TArray<float> Aspect_deg;     // Optional input for redistribution (zeros for now)

	// New grid fields
	TArray<float> SWE;         // mm
	TArray<float> Depth01;     // normalized [0..1] for debug draw

	UPROPERTY(VisibleAnywhere, Category="UnrealSnow|Debug")
	float DebugMeanSWE = 0.f;

	UPROPERTY(VisibleAnywhere, Category="UnrealSnow|Time")
	FDateTime CurrentSimTime;
	int32 StepIndex = 0;
	float LastMeanSWE = 0.0f;
	float LastMeanSWRout = 0.0f;

	void InitGrid();
	void EnsureResources();
	void RunStep(float DeltaSeconds);
};


