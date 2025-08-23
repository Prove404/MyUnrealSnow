#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Components/StaticMeshComponent.h"
// forward include allowed for clarity
#include "Weather/StochasticWeatherProvider.h"
#pragma region ForwardDecls
class UTextureRenderTarget2D;
class AVirtualHeightfieldMesh;
class UMaterialInstanceDynamic;
class UVirtualHeightfieldMeshComponent;
#pragma endregion
#include "ProceduralMeshComponent.h"
#include "Kismet/GameplayStatics.h"
#include "Components/PrimitiveComponent.h"
#include "Materials/MaterialInstanceDynamic.h"
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

	// Debug helpers
	UPROPERTY(EditAnywhere, Category="Snow|Debug")
	float DebugDisplacementBoost = 1.f;

	// Optionally hide any legacy snow plane mesh if attached to this actor
	UPROPERTY(EditAnywhere, Category="Snow|Debug")
	bool bHideLegacySnowPlane = true;

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

	// --- Snow displacement (runtime heightfield) ---
	UPROPERTY(EditAnywhere, Category="Snow|Displacement")
	float SnowDensity = 300.f; // kg m^-3 (use to convert SWE->depth)

	UPROPERTY(EditAnywhere, Category="Snow|Displacement")
	bool bDepthFromSWE = true; // depth = SWE * (1000 / SnowDensity)

	UPROPERTY(EditAnywhere, Category="Snow|Displacement")
	FName SnowWriterActorTag = "SnowRVTWriter"; // actor in level with this tag

	UPROPERTY(EditAnywhere, Category="Snow|Viz")
	FName SnowDepthParam = TEXT("SnowDepthTex"); // material parameter name

	UPROPERTY(EditAnywhere, Category="Snow|Viz")
	FName SnowDepthScaleParam = TEXT("SnowDepthToWorld");

	UPROPERTY(EditAnywhere, Category="Snow|Viz")
	int32 DepthResX = 512;

	UPROPERTY(EditAnywhere, Category="Snow|Viz")
	int32 DepthResY = 512;

	UPROPERTY(Transient)
	UTexture2D* SnowDepthTex = nullptr; // R16F meters

	TArray<float> SnowDepthData;
	FUpdateTextureRegion2D SnowDepthRegion;

	// Writer MID cache
	UPROPERTY(Transient)
	UMaterialInstanceDynamic* SnowWriterMID = nullptr;

	// --- Procedural snow surface mesh ---
	UPROPERTY(VisibleAnywhere, Category="Snow|Mesh")
	UProceduralMeshComponent* SnowMesh = nullptr;

	UPROPERTY(EditAnywhere, Category="Snow|Mesh")
	bool bUseProceduralMesh = false;

	UPROPERTY(EditAnywhere, Category="Snow|Mesh") int32 GridX = 128;
	UPROPERTY(EditAnywhere, Category="Snow|Mesh") int32 GridY = 128;
	UPROPERTY(EditAnywhere, Category="Snow|Mesh") float WorldSizeX = 5000.f;
	UPROPERTY(EditAnywhere, Category="Snow|Mesh") float WorldSizeY = 5000.f;

	// Physics/displacement
	UPROPERTY(EditAnywhere, Category="Snow|Physics") float DisplacementScale = 1.f;
	UPROPERTY(EditAnywhere, Category="Snow|Physics") bool bSmoothDepth = true;
	UPROPERTY(EditAnywhere, Category="Snow|Physics") int32 SmoothKernelRadius = 1;

	// SIMPLE INTERNAL SWE FALLBACK (meters water equivalent)
	TArray<float> SWE_Grid;
	TArray<float> SWE_Temp;

	// ===== Procedural snow mesh buffers =====
	UPROPERTY() TArray<FVector> BaseVerts;     // Base terrain vertices (sampled once)
	UPROPERTY() TArray<FVector> DynVerts;      // Deformed vertices (Base + snow depth)
	UPROPERTY() TArray<FVector> Normals;
	UPROPERTY() TArray<FVector2D> UVs;
	UPROPERTY() TArray<FProcMeshTangent> Tangents;
	UPROPERTY() TArray<int32> Indices;

	// Utility
	int32 VertIndex(int32 ix, int32 iy) const { return iy*(GridX+1) + ix; }
	FVector2D GridToLocalXY(int32 ix, int32 iy) const;

public:
	ASnowSimulationActor();
	virtual void BeginPlay() override;
	virtual void OnConstruction(const FTransform& Xform) override;
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

	// Optional: expose for Blueprints if you want manual hookups
	UFUNCTION(BlueprintCallable, Category="Snow|Viz")
	UTexture2D* GetSnowDepthTexture() const { return SnowDepthTex; }

	void InitSnowDepthTexture();
	void UpdateSnowDepthTexture();
	void WireSnowDepthToVHM()
	{
		UWorld* W = GetWorld(); 
		if (!W || !SnowDepthTex) return;

		TArray<AActor*> Actors;
		UGameplayStatics::GetAllActorsOfClass(W, AActor::StaticClass(), Actors);

		int32 Wired = 0;
		for (AActor* A : Actors)
		{
			if (!A) continue;

			TArray<UPrimitiveComponent*> Prims;
			A->GetComponents<UPrimitiveComponent>(Prims);
			for (UPrimitiveComponent* P : Prims)
			{
				if (!P) continue;

				const FString ClassName = P->GetClass()->GetName(); // "VirtualHeightfieldMeshComponent"
				if (!ClassName.Contains(TEXT("VirtualHeightfieldMeshComponent")))
					continue;

				if (UMaterialInstanceDynamic* MID = P->CreateAndSetMaterialInstanceDynamic(0))
				{
					MID->SetTextureParameterValue(SnowDepthParam, SnowDepthTex);
					MID->SetScalarParameterValue(SnowDepthScaleParam, 500.f);
					P->SetMaterial(0, MID);
					++Wired;
				}
			}
		}
		UE_LOG(LogTemp, Warning, TEXT("[Snow] Wired VHM comps: %d"), Wired);
	}

	// Procedural mesh API
	UFUNCTION(BlueprintCallable, Category="Snow|Mesh")
	void RebuildSnowMesh();

	UFUNCTION(BlueprintCallable, Category="Snow|Mesh")
	float GetSWE(int32 ix, int32 iy) const;

	// Safety: ensure mesh buffers are ready before updating
	bool HasValidMeshBuffers() const { return DynVerts.Num() == (GridX+1)*(GridY+1); }

	// Console exec for quick mesh/bounds debug
	UFUNCTION(Exec)
	void SnowMeshInfo();

	// Console exec helpers
	UFUNCTION(Exec)
	void SnowMove(float X, float Y, float Z);           // set actor world location (cm)
	UFUNCTION(Exec)
	void SnowHere(float SizeMeters=500.0f); // move actor under camera & set size
	UFUNCTION(Exec)
	void SnowSize(float SizeMeters);        // just resize patch
	UFUNCTION(Exec)
	void SnowLogCenter();                   // print sample trace info

	// Debug options
	UPROPERTY(EditAnywhere, Category="Snow|Debug")
	bool bAutoSnapAtBeginPlay = true; // snap once at PIE

	// SWE grid external setter and clamped getter
	UFUNCTION(BlueprintCallable, Category="Snow|SWE")
	void SetSWEAt(int32 X, int32 Y, float SWE_m);
	float GetSWE_Clamped(int32 ix, int32 iy) const;

	// Toggle internal precip growth hack
	UPROPERTY(EditAnywhere, Category="Snow|SWE")
	bool bUseInternalPrecipHack = false;

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

	// SWE helpers
	void EnsureSWEStorage();
	void SmoothSWE();
};


