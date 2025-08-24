#include "Simulation/SnowSimulationActor.h"
#include "UnrealSnowLog.h"

#include "SnowRTUtils.h"
#include "Engine/TextureRenderTarget2D.h"
#include "Weather/StochasticWeatherProvider.h"
#include "UObject/SoftObjectPtr.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Components/SceneComponent.h"
#include "Components/PrimitiveComponent.h"
#include "Materials/Material.h"
#include "Materials/MaterialInterface.h"
#include "DrawDebugHelpers.h"
#include "Engine/CollisionProfile.h"
#include "Simulation/ISnowRedistributor.h"
#include "Simulation/RidgeLeewardRedistributor.h"
#include "Kismet/GameplayStatics.h"
#include "Components/MeshComponent.h"
#include "Engine/StaticMesh.h"
#include "Engine/Texture2D.h"
#include "Engine/StaticMeshActor.h"
#include "HAL/PlatformFilemanager.h"
#include "Misc/Paths.h"
#include "Engine/World.h"
#include "Engine/EngineTypes.h"
#include "TimerManager.h"
#include "Math/Float16.h"
#include "LandscapeProxy.h"
#include "EngineUtils.h"
#include "VirtualHeightfieldMeshComponent.h"

#ifndef SNOW_WITH_LEGACY_PROCMESH
#define SNOW_WITH_LEGACY_PROCMESH 0
#endif

// Forward decl from SnowCS.cpp
void RunSnowCS_RenderTarget(class UTextureRenderTarget2D* RT, class UTextureRenderTarget2D* SWETexture, float InvMaxDepth);

#if !SNOW_WITH_LEGACY_PROCMESH
// Stubs to satisfy any accidental references in cooked BPs
void ASnowSimulationActor::BuildProceduralSnowPlane(){}
void ASnowSimulationActor::UpdateProceduralMesh(){}
#endif

namespace
{
	static constexpr int32 kDefaultGridSize = 256;

	static float ComputeWindSpeed(const FWeatherSample& S)
	{
		return FMath::Sqrt(S.U10_ms * S.U10_ms + S.V10_ms * S.V10_ms);
	}

	static float ComputeWindDirFromDeg(const FWeatherSample& S)
	{
		// Meteorological 'from' direction in degrees, clockwise from north
		const float DirRad = FMath::Atan2(-S.U10_ms, -S.V10_ms);
		float DirDeg = FMath::RadiansToDegrees(DirRad);
		if (DirDeg < 0.0f) DirDeg += 360.0f;
		return DirDeg;
	}
}

// Diagnostics helper to introspect actor components and child VHMs
static void DumpActorComponents(AActor* A, const TCHAR* Tag)
{
    if (!A) { UE_LOG(LogTemp, Warning, TEXT("[Snow] %s: <null actor>"), Tag); return; }
    UE_LOG(LogTemp, Warning, TEXT("[Snow] %s: Actor=%s Class=%s"), Tag, *A->GetName(), *A->GetClass()->GetName());
    TArray<UActorComponent*> Cs; A->GetComponents(Cs);
    UE_LOG(LogTemp, Warning, TEXT("[Snow]   Direct components: %d"), Cs.Num());
    for (UActorComponent* C : Cs) { if (C) UE_LOG(LogTemp, Warning, TEXT("      - %s"), *C->GetClass()->GetName()); }
    TInlineComponentArray<UVirtualHeightfieldMeshComponent*> VHMs(A, /*bIncludeFromChildActors*/ true);
    UE_LOG(LogTemp, Warning, TEXT("[Snow]   Found VHM comps (incl. children): %d"), VHMs.Num());
    for (auto* V : VHMs) { if (V) UE_LOG(LogTemp, Warning, TEXT("      * VHM: %s"), *V->GetName()); }
}

static void LogSWEStats(const TArray<float>& SWE, int32 GX, int32 GY, float SnowDensity, float DisplacementScale)
{
	if (SWE.Num() == 0) return;
	float sum=0, maxSWE=0; for (float v : SWE){ sum+=v; maxSWE = FMath::Max(maxSWE, v); }
	const float mean = sum / SWE.Num();                   // m w.e.
	const float rhoWater = 1000.f;
	const float maxDepth_m = maxSWE * (rhoWater/FMath::Max(1.f, SnowDensity));
	const float maxDisp_cm = maxDepth_m * 100.f * DisplacementScale;

	UE_LOG(LogTemp, Log, TEXT("[Snow] SWE mean=%.4f mwe, max=%.4f mwe -> maxDepth=%.3f m, disp=%.1f cm (scale=%.2f)"),
		mean, maxSWE, maxDepth_m, maxDisp_cm, DisplacementScale);
}

// Helper to disable and hide a component safely (used for legacy procedural mesh cleanup)
static void DisableAndHideComponent(UActorComponent* Comp)
{
    if (!Comp) return;
    Comp->SetComponentTickEnabled(false);
    if (UPrimitiveComponent* Prim = Cast<UPrimitiveComponent>(Comp))
    {
        Prim->SetHiddenInGame(true);
        Prim->SetVisibility(false, true);
        Prim->SetCollisionEnabled(ECollisionEnabled::NoCollision);
    }
}

// Remove any ProceduralMeshComponent instances without directly depending on its header
static void RemoveLegacyProcMeshComponents(AActor* Owner)
{
    if (!Owner) return;
    TArray<UActorComponent*> AllComps;
    Owner->GetComponents(AllComps);
    for (UActorComponent* C : AllComps)
    {
        if (!C) continue;
        if (C->GetClass()->GetName().Contains(TEXT("ProceduralMeshComponent")))
        {
            DisableAndHideComponent(C);
            C->DestroyComponent(true);
        }
    }
}

ASnowSimulationActor::ASnowSimulationActor()
{
	PrimaryActorTick.bCanEverTick = true;

	Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	SetRootComponent(Root);
}

void ASnowSimulationActor::PostLoad()
{
	Super::PostLoad();
	RemoveLegacyProcMeshComponents(this);
}

void ASnowSimulationActor::OnConstruction(const FTransform& Xform)
{
	Super::OnConstruction(Xform);
}

void ASnowSimulationActor::BeginPlay()
{
	Super::BeginPlay();

	// Safety: remove any legacy procedural mesh components
	RemoveLegacyProcMeshComponents(this);

	if (bAutoSnapAtBeginPlay)
	{
		FTimerHandle H;
		GetWorldTimerManager().SetTimer(H, [this]()
		{
			SnowHere(500.f);
		}, 0.05f, false);
	}
	InitGrid();
	if (!WeatherProvider && WeatherProviderClass.IsValid())
	{
		UClass* Cls = WeatherProviderClass.LoadSynchronous();
		if (Cls)
		{
			WeatherProvider = NewObject<UStochasticWeatherProvider>(this, Cls, NAME_None, RF_Transactional);
		}
	}
	if (!Redistributor)
	{
		Redistributor = NewObject<URidgeLeewardRedistributor>(this, URidgeLeewardRedistributor::StaticClass(), NAME_None, RF_Transactional);
	}
	CurrentSimTime = StartTime;
	UE_LOG(LogUnrealSnow, Display, TEXT("Grid init: %d x %d (%d cells)"), GridSizeXY, GridSizeXY, NumCells);
	EnsureResources();

	// Debug plane hookup
	if (DebugMaterial && SnowDepthRT)
	{
		DebugMID = UMaterialInstanceDynamic::Create(DebugMaterial, this);
		DebugMID->SetTextureParameterValue(DebugTextureParam, SnowDepthRT);
		if (DebugPlane)
		{
			DebugPlane->SetMaterial(0, DebugMID);
		}
	}

	// Initialize runtime depth texture and bind to writer
	if (!SnowDepthTex) InitSnowDepthTexture();
	UpdateSnowDepthTexture();
	UE_LOG(LogTemp, Warning, TEXT("[Snow] DepthTex init: %dx%d"), DepthResX, DepthResY);
	WireSnowDepthToVHM();
	{
		TArray<AActor*> TaggedActors;
		UGameplayStatics::GetAllActorsWithTag(GetWorld(), SnowWriterActorTag, TaggedActors);
		if (TaggedActors.Num() > 0)
		{
			AActor* Writer = TaggedActors[0];
			UMeshComponent* MeshComp = Writer ? Writer->FindComponentByClass<UMeshComponent>() : nullptr;
			if (MeshComp)
			{
				UMaterialInterface* BaseMat = MeshComp->GetMaterial(0);
				if (BaseMat)
				{
					SnowWriterMID = UMaterialInstanceDynamic::Create(BaseMat, this);
					MeshComp->SetMaterial(0, SnowWriterMID);
					if (SnowWriterMID && SnowDepthTex)
					{
						// Bind the dynamic R16F depth texture
						SnowWriterMID->SetTextureParameterValue(SnowDepthParam, SnowDepthTex);
						// Verify texture parameter binding
						{
							UTexture* BoundTex = nullptr;
							const bool bGot = SnowWriterMID->GetTextureParameterValue(FName("SnowDepthTex"), BoundTex);
							if (bGot && BoundTex == SnowDepthTex)
							{
								UE_LOG(LogTemp, Warning, TEXT("[Snow] Texture param successfully bound."));
							}
							else
							{
								UE_LOG(LogTemp, Error, TEXT("[Snow] Texture param binding FAILED."));
							}
						}
						// Optional: also push the scalar params so you don’t set them by hand
						SnowWriterMID->SetScalarParameterValue(FName("DisplacementScale"), DisplacementScale);
						SnowWriterMID->SetScalarParameterValue(FName("DepthRangeMeters"), DepthRangeMeters);
						const float DomainSizeCm = WorldSizeX;
						const FVector2D OriginXY(GetActorLocation().X - DomainSizeCm * 0.5f,
							GetActorLocation().Y - DomainSizeCm * 0.5f);
						SnowWriterMID->SetVectorParameterValue(FName("SnowOriginWS"),
							FLinearColor(OriginXY.X, OriginXY.Y, 0, 0));
						SnowWriterMID->SetScalarParameterValue(FName("SnowInvSize"), (DomainSizeCm != 0.f) ? (1.0f / DomainSizeCm) : 0.f);
					}
				}
			}
		}
	}

	// Procedural mesh path removed
	EnsureSWEStorage();

	// Hide legacy primitives except SnowMesh so only the procedural surface is visible
	{
		TArray<UPrimitiveComponent*> PrimComps;
		GetComponents<UPrimitiveComponent>(PrimComps);
		for (UPrimitiveComponent* P : PrimComps)
		{
			if (P)
			{
				P->SetHiddenInGame(true);
				P->SetVisibility(false, true);
			}
		}
	}
}

void ASnowSimulationActor::InitGrid()
{
	const int32 N = FMath::Max(1, GridSizeXY);
	NumCells = N * N;
	SWE.SetNumZeroed(NumCells);
	Depth01.SetNumZeroed(NumCells);

	if (!SnowDepthRT || SnowDepthRT->SizeX != N || SnowDepthRT->SizeY != N)
	{
		if (!SnowDepthRT)
		{
			SnowDepthRT = NewObject<UTextureRenderTarget2D>(this);
			SnowDepthRT->bAutoGenerateMips = false;
			SnowDepthRT->ClearColor = FLinearColor::Black;
		}
		SnowDepthRT->RenderTargetFormat = RTF_R32f;  // single-channel float
		SnowDepthRT->InitAutoFormat(N, N);
		SnowDepthRT->UpdateResourceImmediate(true);
	}
}

void ASnowSimulationActor::EnsureResources()
{
	int32 W = GridSize.X;
	int32 H = GridSize.Y;
	if (W <= 0 || H <= 0)
	{
		W = kDefaultGridSize;
		H = kDefaultGridSize;
		GridSize = FIntPoint(W, H);
	}

	const int32 Num = W * H;
	if (SnowWaterEq_mm.Num() != Num)
	{
		SnowWaterEq_mm.SetNumZeroed(Num);
	}
	if (Slope_deg.Num() != Num)
	{
		Slope_deg.SetNumZeroed(Num);
	}
	if (Aspect_deg.Num() != Num)
	{
		Aspect_deg.SetNumZeroed(Num);
	}

	USnowRTUtils::InitFloatRT(SnowDepthRT, W, PF_R32_FLOAT);
	USnowRTUtils::InitFloatRT(SWROutRT, W, PF_R32_FLOAT);
}

void ASnowSimulationActor::EnsureSWEStorage()
{
	const int32 N = GridX * GridY;
	if (SWE_Grid.Num() != N) { SWE_Grid.SetNumZeroed(N); }
	if (SWE_Temp.Num() != N) { SWE_Temp.SetNumZeroed(N); }
}

void ASnowSimulationActor::SetSWEAt(int32 X, int32 Y, float SWE_m)
{
	X = FMath::Clamp(X, 0, GridX - 1);
	Y = FMath::Clamp(Y, 0, GridY - 1);
	EnsureSWEStorage();
	SWE_Grid[Y * GridX + X] = FMath::Max(0.f, SWE_m);
}

float ASnowSimulationActor::GetSWE_Clamped(int32 ix, int32 iy) const
{
	ix = FMath::Clamp(ix, 0, GridX - 1);
	iy = FMath::Clamp(iy, 0, GridY - 1);
	if (SWE_Grid.Num() == GridX * GridY)
	{
		return SWE_Grid[iy * GridX + ix];
	}
	return 0.f;
}

void ASnowSimulationActor::SmoothSWE()
{
	if (!bSmoothDepth || SmoothKernelRadius <= 0) return;
	const int32 R = SmoothKernelRadius;
	const int32 N = GridX * GridY;
	if (SWE_Temp.Num() != N) SWE_Temp.SetNumZeroed(N);

	for (int32 y = 0; y < GridY; ++y)
	{
		for (int32 x = 0; x < GridX; ++x)
		{
			float acc = 0.f; int32 cnt = 0;
			for (int32 dy = -R; dy <= R; ++dy)
			for (int32 dx = -R; dx <= R; ++dx)
			{
				const int32 xx = FMath::Clamp(x + dx, 0, GridX - 1);
				const int32 yy = FMath::Clamp(y + dy, 0, GridY - 1);
				acc += SWE_Grid[yy * GridX + xx]; ++cnt;
			}
			SWE_Temp[y * GridX + x] = acc / (float)cnt;
		}
	}
	SWE_Grid = SWE_Temp;
}

void ASnowSimulationActor::Tick(float DeltaSeconds)
{
	Super::Tick(DeltaSeconds);
	if (NumCells == 0)
	{
		UE_LOG(LogUnrealSnow, Error, TEXT("NumCells is 0; grid not initialized."));
		return;
	}

	const double dHours = HoursPerSecond * DeltaSeconds;
	CurrentSimTime += FTimespan::FromSeconds(dHours * 3600.0);

	float Precip_mmph = 0.f;
	if (WeatherProvider)
	{
		FWeatherSample S{};
		const bool bOk = WeatherProvider->GetForTime_Implementation(CurrentSimTime, S);
		if (bOk) Precip_mmph = S.Precip_mmph;
	}
	if (bForceConstantSnow) { Precip_mmph += DebugSnow_mmph; }

	for (int32 i = 0; i < NumCells; ++i)
	{
		SWE[i] += Precip_mmph * (float)dHours;
	}
	// Propagate to procedural mesh SWE grid (convert mm -> meters)
	{
		const int32 N = GridSizeXY;
		if (N > 1)
		{
			for (int32 i = 0; i < NumCells; ++i)
			{
				const int32 ix = i % N;
				const int32 iy = i / N;
				const int32 mx = (GridX > 1) ? FMath::Clamp(FMath::RoundToInt(ix * (GridX - 1) / float(N - 1)), 0, GridX - 1) : 0;
				const int32 my = (GridY > 1) ? FMath::Clamp(FMath::RoundToInt(iy * (GridY - 1) / float(N - 1)), 0, GridY - 1) : 0;
				const float swe_m = SWE[i] * 0.001f;
				SetSWEAt(mx, my, swe_m);
			}
		}
	}

	const float rho = 100.f;            // kg/m^3 fake bulk density
	for (int32 i = 0; i < NumCells; ++i)
	{
		const float depth_m = SWE[i] / rho;        // SWE(mm)/rho
		Depth01[i] = FMath::Clamp(depth_m / 1.0f, 0.f, 1.f); // 1 m full white for debug
	}

	USnowRTUtils::WriteArrayToRT(SnowDepthRT, Depth01);
	if (DebugMID && SnowDepthRT)
	{
		DebugMID->SetTextureParameterValue(DebugTextureParam, SnowDepthRT);
	}

	// Update runtime height texture and ensure writer binding is alive
	UpdateSnowDepthTexture();
	if (!SnowWriterMID)
	{
		TArray<AActor*> TaggedActors;
		UGameplayStatics::GetAllActorsWithTag(GetWorld(), SnowWriterActorTag, TaggedActors);
		if (TaggedActors.Num() > 0)
		{
			AActor* Writer = TaggedActors[0];
			UMeshComponent* MeshComp = Writer ? Writer->FindComponentByClass<UMeshComponent>() : nullptr;
			if (MeshComp)
			{
				UMaterialInterface* BaseMat = MeshComp->GetMaterial(0);
				if (BaseMat)
				{
					SnowWriterMID = UMaterialInstanceDynamic::Create(BaseMat, this);
					MeshComp->SetMaterial(0, SnowWriterMID);
				}
			}
		}
	}
	if (SnowWriterMID && SnowDepthTex)
	{
		SnowWriterMID->SetTextureParameterValue(SnowDepthParam, SnowDepthTex);
	}

	// --- TEMPORARY demo growth (disable by default) ---
	if (bUseInternalPrecipHack)
	{
		EnsureSWEStorage();
		const float PrecipRate_mps = 0.00001f; // 0.01 mm/s placeholder
		for (int32 i = 0; i < SWE_Grid.Num(); ++i)
			SWE_Grid[i] = FMath::Max(0.f, SWE_Grid[i] + PrecipRate_mps * DeltaSeconds);
	}
	// --- END TEMP ---

	SmoothSWE();

	// Assert depth values sanity to catch the km bug early
	float MaxDepthMeters = 0.f;
	{
		const float rhoWater_local = 1000.f;
		if (SWE_Grid.Num() == GridX * GridY)
		{
			for (int32 i = 0; i < SWE_Grid.Num(); ++i)
			{
				const float depth_m = SWE_Grid[i] * (rhoWater_local / FMath::Max(1.f, SnowDensity));
				MaxDepthMeters = FMath::Max(MaxDepthMeters, depth_m);
			}
		}
	}
	ensureAlwaysMsgf(MaxDepthMeters < 50.f, TEXT("[Snow] Depth sanity check failed: %.3fm (expected < 50m). Check SWE->depth scale."), MaxDepthMeters);

	LogSWEStats(SWE_Grid, GridX, GridY, SnowDensity, DisplacementScale);

#if SNOW_WITH_LEGACY_PROCMESH
	if (bEnableLegacyProceduralMesh)
	{
		if (!HasValidMeshBuffers())
		{
			UE_LOG(LogTemp, Warning, TEXT("[Snow] Mesh buffers not ready; procedural path removed"));
			if (!HasValidMeshBuffers()) return;
		}

		const float rhoWater = 1000.f;
		const int32 VX = GridX + 1;
		const int32 VY = GridY + 1;
		if (DynVerts.Num() == VX * VY)
		{
			for (int32 iy = 0; iy < VY; ++iy)
			{
				for (int32 ix = 0; ix < VX; ++ix)
				{
					const int32 vi = VertIndex(ix, iy);
					FVector v = BaseVerts.IsValidIndex(vi) ? BaseVerts[vi] : DynVerts[vi];
					const int32 cx = FMath::Clamp(ix, 0, GridX - 1);
					const int32 cy = FMath::Clamp(iy, 0, GridY - 1);
					const float SWE_m = GetSWE_Clamped(cx, cy);
					const float depth_m = SWE_m * (rhoWater / FMath::Max(1.f, SnowDensity));
					v.Z = BaseVerts[vi].Z + depth_m * 100.f * DisplacementScale * DebugDisplacementBoost;
					DynVerts[vi] = v;
				}
			}
			if (SnowMesh)
			{
				if (DynVerts.Num() != (GridX+1)*(GridY+1))
				{
					UE_LOG(LogTemp, Warning, TEXT("[Snow] DynVerts size %d != expected %d; procedural path removed"),
						DynVerts.Num(), (GridX+1)*(GridY+1));
					if (DynVerts.Num() != (GridX+1)*(GridY+1)) return;
				}
				SnowMesh->UpdateMeshSection_LinearColor(0, DynVerts, Normals, UVs, TArray<FLinearColor>(), Tangents);
			}
		}
	}
#endif
}

void ASnowSimulationActor::RunStep(float DeltaSeconds)
{
    if (!WeatherProvider)
    {
        return;
    }

    EnsureResources();

    FWeatherSample Sample;
    bool bHasWeather = false;
    static thread_local bool bInProviderCall = false;
    if (!bInProviderCall && WeatherProvider)
    {
        bInProviderCall = true;
        bHasWeather = WeatherProvider->GetForTime_Implementation(CurrentSimTime, Sample);
        bInProviderCall = false;
    }
    if (!bHasWeather)
    {
        return;
    }

    const int32 W = GridSize.X;
    const int32 H = GridSize.Y;
    const int32 Num = W * H;

	// 1) Accumulation with dt from tick (in hours)
	const float DtHours = FMath::Max(DeltaSeconds, 0.0f) / 3600.0f;
	float Precip_mmph = FMath::Max(Sample.Precip_mmph, 0.0f);
	if (bForceConstantSnow) { Precip_mmph += DebugSnow_mmph; }
	const float dSWE = Precip_mmph * DtHours; // [mm]
	// Accumulate SWE on CPU as source for both paths
	for (int32 i = 0; i < Num; ++i)
	{
		SnowWaterEq_mm[i] += dSWE;
	}
	// Propagate to procedural mesh SWE grid (convert mm -> meters)
	{
		if (W > 1 && H > 1)
		{
			for (int32 i = 0; i < Num; ++i)
			{
				const int32 ix = i % W;
				const int32 iy = i / W;
				const int32 mx = (GridX > 1) ? FMath::Clamp(FMath::RoundToInt(ix * (GridX - 1) / float(W - 1)), 0, GridX - 1) : 0;
				const int32 my = (GridY > 1) ? FMath::Clamp(FMath::RoundToInt(iy * (GridY - 1) / float(H - 1)), 0, GridY - 1) : 0;
				const float swe_m = SnowWaterEq_mm[i] * 0.001f;
				SetSWEAt(mx, my, swe_m);
			}
		}
	}

	// 2) Optional redistribution
	if (Redistributor)
	{
		const float WindSpeed = ComputeWindSpeed(Sample);
		const float WindDirFromDeg = ComputeWindDirFromDeg(Sample);
		ISnowRedistributor::Execute_Apply(
			Redistributor,
			GridSize,
			SnowWaterEq_mm,
			Slope_deg,
			Aspect_deg,
			FVector2D(WindDirFromDeg, 0.0f),
			WindSpeed);
	}

	// 3) Convert SWE to snow depth [m] using temporary constant density 100 kg m^-3
	TArray<float> Depth_m;
	Depth_m.SetNumUninitialized(Num);
	for (int32 i = 0; i < Num; ++i)
	{
		// SWE [mm] -> [m water] = mm / 1000; then divide by density ratio 100/1000 => /100
		Depth_m[i] = SnowWaterEq_mm[i] / 100.0f;
	}

	// 4) Normalize to [0,1] for visualization
	float MaxDepth = 0.0f;
	for (int32 i = 0; i < Num; ++i)
	{
		MaxDepth = FMath::Max(MaxDepth, Depth_m[i]);
	}
	TArray<float> Normalized;
	Normalized.SetNumUninitialized(Num);
	float InvMaxDepth = 1.0f;
	if (MaxDepth > KINDA_SMALL_NUMBER)
	{
		InvMaxDepth = 1.0f / MaxDepth;
		for (int32 i = 0; i < Num; ++i)
		{
			Normalized[i] = FMath::Clamp(Depth_m[i] * InvMaxDepth, 0.0f, 1.0f);
		}
	}
	else
	{
		for (int32 i = 0; i < Num; ++i)
		{
			Normalized[i] = 0.0f;
		}
	}

	// 5) CPU write or GPU dispatch (SWE->Depth normalized)
	if (!bUseGPU)
	{
		USnowRTUtils::WriteArrayToRT(SnowDepthRT, Normalized);
	}
	else
	{
		// Upload SWE to a temporary RT and run compute to convert and normalize
		static UTextureRenderTarget2D* SWERT = nullptr;
		if (!SWERT || SWERT->SizeX != W || SWERT->SizeY != H)
		{
			USnowRTUtils::InitFloatRT(SWERT, W, PF_R32_FLOAT);
		}
		USnowRTUtils::WriteArrayToRT(SWERT, SnowWaterEq_mm);
		RunSnowCS_RenderTarget(SnowDepthRT, SWERT, InvMaxDepth);
	}

	if (DebugMID && SnowDepthRT)
	{
		DebugMID->SetTextureParameterValue(DebugTextureParam, SnowDepthRT);
	}

	// 5b) Compute SWRout = alpha * SWdown_eff (placeholder SWdown_eff=1)
	static TArray<float> SnowAgeDays;
	if (SnowAgeDays.Num() != Num)
	{
		SnowAgeDays.SetNumZeroed(Num);
	}
	const float DtDays = FMath::Max(DeltaSeconds, 0.0f) / 86400.0f;
	for (int32 i = 0; i < Num; ++i)
	{
		SnowAgeDays[i] += DtDays;
		if (dSWE > 0.0f)
		{
			SnowAgeDays[i] = 0.0f;
		}
	}
	const float Decay = (AlbedoDecayDays > KINDA_SMALL_NUMBER) ? (1.0f / AlbedoDecayDays) : 1.0f;
	TArray<float> SWRout;
	SWRout.SetNumUninitialized(Num);
	for (int32 i = 0; i < Num; ++i)
	{
		const float t = SnowAgeDays[i];
		const float alpha = AlbedoOld + (AlbedoFresh - AlbedoOld) * FMath::Exp(-t * Decay);
		const float SWdown_eff = 1.0f; // placeholder (SWdown * shadow)
		SWRout[i] = alpha * SWdown_eff;
	}
	USnowRTUtils::WriteArrayToRT(SWROutRT, SWRout);

	// Metrics
	{
		double sumSWE = 0.0;
		for (int32 i = 0; i < Num; ++i) sumSWE += SnowWaterEq_mm[i];
		LastMeanSWE = (Num > 0) ? (float)(sumSWE / Num) : 0.0f;
		// SWRout mean equals alpha mean (since SWdown_eff=1)
		double sumAlpha = 0.0;
		for (int32 i = 0; i < Num; ++i) sumAlpha += Normalized[i];
		LastMeanSWRout = (Num > 0) ? (float)(sumAlpha / Num) : 0.0f;
	}

	// Advance step index
	++StepIndex;
}

void ASnowSimulationActor::SetSimRunning(bool bRun)
{
	SetActorTickEnabled(bRun);
}

void ASnowSimulationActor::SimulateStep(float DeltaSeconds)
{
	RunStep(DeltaSeconds);
}

void ASnowSimulationActor::InitSnowDepthTexture()
{
	SnowDepthTex = UTexture2D::CreateTransient(DepthResX, DepthResY, PF_R16F);
	SnowDepthTex->SRGB = false;
	SnowDepthTex->CompressionSettings = TC_HDR;
	SnowDepthTex->Filter = TF_Nearest;
	SnowDepthTex->AddressX = TA_Clamp; SnowDepthTex->AddressY = TA_Clamp;
	SnowDepthTex->UpdateResource();

	// Zero-initialize first mip once after creation
	if (SnowDepthTex && SnowDepthTex->GetPlatformData() && SnowDepthTex->GetPlatformData()->Mips.Num() > 0)
	{
		FTexture2DMipMap& Mip = SnowDepthTex->GetPlatformData()->Mips[0];
		void* Data = Mip.BulkData.Lock(LOCK_READ_WRITE);
		FMemory::Memset(Data, 0, Mip.BulkData.GetBulkDataSize());
		Mip.BulkData.Unlock();
		SnowDepthTex->UpdateResource();
	}

	SnowDepthRegion = FUpdateTextureRegion2D(0, 0, 0, 0, DepthResX, DepthResY);
}

void ASnowSimulationActor::UpdateSnowDepthTexture()
{
	if (!SnowDepthTex) return;

	if (SnowDepthTex->GetSizeX() != DepthResX || SnowDepthTex->GetSizeY() != DepthResY)
	{
		InitSnowDepthTexture();
	}

	const int32 W = DepthResX, H = DepthResY;
	const float rho_w = 1000.f; // kg/m^3
	const float invRhoSnow = (SnowDensity > 1.f) ? (1.0f / SnowDensity) : (1.0f / 300.f);

	float MaxDepthMetersThisTick = 0.f;

	if (!SnowDepthTex->GetPlatformData() || SnowDepthTex->GetPlatformData()->Mips.Num() == 0)
	{
		return;
	}

	FTexture2DMipMap& Mip = SnowDepthTex->GetPlatformData()->Mips[0];
	uint8* Raw = static_cast<uint8*>(Mip.BulkData.Lock(LOCK_READ_WRITE));
	FFloat16* Out = reinterpret_cast<FFloat16*>(Raw);
	const int32 N = W * H;

	for (int32 y = 0; y < H; ++y)
	{
		const float ty = (float)y / (H - 1);
		for (int32 x = 0; x < W; ++x)
		{
			const float tx = (float)x / (W - 1);
			const int32 gx = FMath::Clamp(FMath::RoundToInt(tx * (GridX - 1)), 0, GridX - 1);
			const int32 gy = FMath::Clamp(FMath::RoundToInt(ty * (GridY - 1)), 0, GridY - 1);

			const int32 i = y * W + x;
			const float SWE_mwe = GetSWE_Clamped(gx, gy);
			float depth_m = SWE_mwe * rho_w * invRhoSnow; // meters of snow
			depth_m = FMath::Clamp(depth_m, 0.f, MaxClampMeters);
			MaxDepthMetersThisTick = FMath::Max(MaxDepthMetersThisTick, depth_m);

			if (bNormalizeDepthTexture)
			{
				const float norm = FMath::Clamp(depth_m / FMath::Max(DepthRangeMeters, 0.01f), 0.f, 1.f);
				Out[i] = FFloat16(norm);
			}
			else
			{
				Out[i] = FFloat16(depth_m); // store meters directly
			}
		}
	}

	Mip.BulkData.Unlock();
	SnowDepthTex->UpdateResource();

	ensureAlwaysMsgf(MaxDepthMetersThisTick <= MaxClampMeters + 0.01f,
		TEXT("[Snow] Depth sanity failed: %.3fm > MaxClampMeters=%.3fm. Check SWE->depth conversion."),
		MaxDepthMetersThisTick, MaxClampMeters);

	// Improved logging
	if (bNormalizeDepthTexture)
	{
		const float normMax = FMath::Clamp(MaxDepthMetersThisTick / FMath::Max(DepthRangeMeters, 0.01f), 0.f, 1.f);
		UE_LOG(LogTemp, Warning, TEXT("[Snow] DepthTex max=%.3fm (norm=%.3f of %.2fm range)"),
			MaxDepthMetersThisTick, normMax, DepthRangeMeters);
	}
	else
	{
		UE_LOG(LogTemp, Warning, TEXT("[Snow] DepthTex max=%.3fm (raw meters in R16F)"),
			MaxDepthMetersThisTick);
	}
}

void ASnowSimulationActor::WireSnowDepthToVHM()
{
	if (!SnowDepthTex) { UE_LOG(LogTemp, Error, TEXT("[Snow] Wire: DepthTex is null")); return; }

	auto BindOne = [&](UVirtualHeightfieldMeshComponent* VHM)->int32
	{
		if (!VHM) return 0;
		
		// --- compute domain (unchanged) ---
		float DomainSizeCm = WorldSizeX;
		FVector2D OriginXY(GetActorLocation().X - DomainSizeCm*0.5f,
			GetActorLocation().Y - DomainSizeCm*0.5f);
		{
			TArray<AActor*> Landscapes;
			UGameplayStatics::GetAllActorsOfClass(GetWorld(), ALandscapeProxy::StaticClass(), Landscapes);
			float Best = TNumericLimits<float>::Max(); ALandscapeProxy* LBest = nullptr;
			for (AActor* LA : Landscapes)
			{
				const float D2 = FVector::DistSquared(LA->GetActorLocation(), GetActorLocation());
				if (D2 < Best) { Best = D2; LBest = Cast<ALandscapeProxy>(LA); }
			}
			if (LBest)
			{
				FVector O, Ext; LBest->GetActorBounds(false, O, Ext);
				DomainSizeCm = FMath::Max(Ext.X*2.f, Ext.Y*2.f);
				OriginXY = FVector2D(O.X - Ext.X, O.Y - Ext.Y);
			}
		}

		// ---- ensure slot 0 exists and has a base material ----
		int32 NumSlots = VHM->GetNumMaterials();
		UE_LOG(LogTemp, Warning, TEXT("[Snow] VHM %s: NumSlots=%d"), *VHM->GetName(), NumSlots);

		UPrimitiveComponent* PrimComp = static_cast<UPrimitiveComponent*>(VHM);
		UMaterialInterface* BaseMtl = (NumSlots > 0) ? PrimComp->GetMaterial(0) : nullptr;
		if (!BaseMtl)
		{
			UMaterialInterface* Fallback = DefaultVHMMaterial.LoadSynchronous();
			if (!Fallback)
			{
				UE_LOG(LogTemp, Error, TEXT("[Snow] VHM %s has no material; set DefaultVHMMaterial in SnowSimulationActor!"),
					*VHM->GetName());
				return 0;
			}
			PrimComp->SetMaterial(0, Fallback);
			BaseMtl = PrimComp->GetMaterial(0);
			UE_LOG(LogTemp, Warning, TEXT("[Snow] Set default material %s on %s"), *Fallback->GetName(), *VHM->GetName());
		}

		// ---- create MID and bind params ----
		UMaterialInstanceDynamic* MID = Cast<UMaterialInstanceDynamic>(BaseMtl);
		if (!MID)
		{
			MID = UMaterialInstanceDynamic::Create(BaseMtl, this);
			PrimComp->SetMaterial(0, MID);
			UE_LOG(LogTemp, Warning, TEXT("[Snow] Created MID on %s slot 0 from %s"), *VHM->GetName(), *BaseMtl->GetName());
		}

		MID->SetTextureParameterValue("SnowDepthTex", SnowDepthTex);
		MID->SetScalarParameterValue("DepthRangeMeters", DepthRangeMeters);
		MID->SetScalarParameterValue("DisplacementScale", DisplacementScale);
		MID->SetVectorParameterValue("SnowOriginWS", FLinearColor(OriginXY.X, OriginXY.Y, 0, 0));
		MID->SetScalarParameterValue("SnowInvSize", 1.0f / FMath::Max(DomainSizeCm, 1.f));

		UTexture* Bound = nullptr;
		MID->GetTextureParameterValue(FMaterialParameterInfo("SnowDepthTex"), Bound);
		UE_LOG(LogTemp, Warning, TEXT("[Snow] Bound=%s Size=%.0fcm Origin=(%.0f,%.0f) InvSize=%.8f"),
			Bound ? *Bound->GetName() : TEXT("NONE"),
			DomainSizeCm, OriginXY.X, OriginXY.Y, 1.0f/FMath::Max(DomainSizeCm,1.f));

		return 1;
	};

	// ----- Try the explicit actor first (with correct component enumeration)
	int32 Total = 0;
	if (TargetVHMActor.IsValid())
	{
		AActor* A = TargetVHMActor.Get();
		UE_LOG(LogTemp, Warning, TEXT("[Snow] TargetVHMActor: %s (%s)"), *A->GetName(), *A->GetClass()->GetName());
		// Correct API: collect VHM comps including child actors
		TArray<UVirtualHeightfieldMeshComponent*> VHMOnActor;
		A->GetComponents<UVirtualHeightfieldMeshComponent>(VHMOnActor, /*bIncludeFromChildActors*/ true);
		UE_LOG(LogTemp, Warning, TEXT("[Snow]   VHM comps on actor (incl. children): %d"), VHMOnActor.Num());
		for (auto* C : VHMOnActor) { Total += BindOne(C); }
		if (Total > 0) { UE_LOG(LogTemp, Warning, TEXT("[Snow] Wired via TargetVHMActor (%d comps)"), Total); }
	}

	// ----- Global fallback if still zero
	if (Total == 0)
	{
		for (TActorIterator<AActor> It(GetWorld()); It; ++It)
		{
			TInlineComponentArray<UVirtualHeightfieldMeshComponent*> VHMOnActor(*It, /*bIncludeFromChildActors*/ true);
			for (auto* C : VHMOnActor) { Total += BindOne(C); }
		}
		UE_LOG(LogTemp, Warning, TEXT("[Snow] Wired VHM comps (global): %d"), Total);
	}

	// ----- Retry once if still nothing (VHM may init after us)
	if (Total == 0 && bRetryWireIfMissing)
	{
		UE_LOG(LogTemp, Warning, TEXT("[Snow] No VHM wired; retrying in %.2fs"), WireRetryDelay);
		FTimerHandle H; GetWorldTimerManager().SetTimer(H, this, &ASnowSimulationActor::WireSnowDepthToVHM, WireRetryDelay, false);
	}
}

FVector2D ASnowSimulationActor::GridToLocalXY(int32 ix, int32 iy) const
{
	const float fx = (static_cast<float>(ix) / static_cast<float>(GridX) - 0.5f) * WorldSizeX;
	const float fy = (static_cast<float>(iy) / static_cast<float>(GridY) - 0.5f) * WorldSizeY;
	return FVector2D(fx, fy);
}

float ASnowSimulationActor::GetSWE(int32 ix, int32 iy) const
{
	const int32 Wsrc0 = GridSize.X > 0 ? GridSize.X : GridSizeXY;
	const int32 Hsrc0 = GridSize.Y > 0 ? GridSize.Y : GridSizeXY;
	const int32 Wsrc = Wsrc0;
	const int32 Hsrc = Hsrc0;
	const bool bUseFull = SnowWaterEq_mm.Num() == Wsrc * Hsrc;
	const TArray<float>& Src = bUseFull ? SnowWaterEq_mm : SWE;
	if (Wsrc <= 0 || Hsrc <= 0 || Src.Num() != Wsrc * Hsrc)
	{
		return 0.0f;
	}
	const int32 sx = (Wsrc > 1) ? FMath::Clamp(FMath::RoundToInt(ix * (Wsrc - 1) / float(GridX)), 0, Wsrc - 1) : 0;
	const int32 sy = (Hsrc > 1) ? FMath::Clamp(FMath::RoundToInt(iy * (Hsrc - 1) / float(GridY)), 0, Hsrc - 1) : 0;
	const int32 sIdx = sy * Wsrc + sx;
	const float swe_mm = (sIdx >= 0 && sIdx < Src.Num()) ? Src[sIdx] : 0.0f;
	return swe_mm * 0.001f; // meters water equivalent
}

// Legacy rebuild guarded by compile-time switch
#if SNOW_WITH_LEGACY_PROCMESH
void ASnowSimulationActor::RebuildSnowMesh()
{
	UE_LOG(LogTemp, Warning, TEXT("[Snow] RebuildSnowMesh() start Grid=%dx%d World=%.0fx%.0f cm"),
		GridX, GridY, WorldSizeX, WorldSizeY);

	if (GridX < 2 || GridY < 2 || WorldSizeX <= 0 || WorldSizeY <= 0)
	{
		UE_LOG(LogTemp, Error, TEXT("[Snow] Invalid grid/size. GridX,Y>=2 and WorldSize>0 required."));
		return;
	}

	if (!LegacyProcMesh)
	{
		return;
	}

	const int32 VX = GridX + 1;
	const int32 VY = GridY + 1;
	BaseVerts.SetNum(VX * VY);
	DynVerts.SetNum(VX * VY);
	Normals.SetNum(VX * VY);
	UVs.SetNum(VX * VY);
	Tangents.SetNum(VX * VY);
	Indices.Reset();
	Indices.Reserve(GridX * GridY * 6);

	UWorld* W = GetWorld();
	if (!W) { UE_LOG(LogTemp, Error, TEXT("[Snow] World==nullptr")); return; }
	FCollisionQueryParams Q(FName(TEXT("SnowBaseTrace")), /*bTraceComplex=*/false, this);
	Q.AddIgnoredActor(this);
	if (LegacyProcMesh) { Q.AddIgnoredComponent(LegacyProcMesh); }
	Q.bFindInitialOverlaps = false;
	FCollisionObjectQueryParams Obj;
	Obj.AddObjectTypesToQuery(ECC_WorldStatic);
	Obj.AddObjectTypesToQuery(ECC_WorldDynamic);

	int32 hitCount = 0;

	const FVector Center = GetActorLocation();

	for (int32 iy = 0; iy < VY; ++iy)
	{
		for (int32 ix = 0; ix < VX; ++ix)
		{
			const FVector2D xy = GridToLocalXY(ix, iy);
			const float WX = Center.X + xy.X;
			const float WY = Center.Y + xy.Y;
			const FVector StartWorld(WX, WY, Center.Z + 1000000.f);
			const FVector EndWorld  (WX, WY, Center.Z - 1000000.f);
			FHitResult Hit;
			float zLocal = 0.f; FVector nLocal = FVector::UpVector;
			if (W && W->LineTraceSingleByObjectType(Hit, StartWorld, EndWorld, Obj, Q) && Hit.bBlockingHit)
			{
				zLocal = Hit.ImpactPoint.Z - Center.Z;
				nLocal = Hit.ImpactNormal.GetSafeNormal();
				++hitCount;

				// Draw a few corner markers in world space
				if ((ix==0||ix==VX-1) && (iy==0||iy==VY-1))
					DrawDebugSphere(W, Hit.ImpactPoint + FVector(0,0,5), 25.f, 12, FColor::Cyan, false, 10.f, 0, 2.f);
			}
			const int32 vi = VertIndex(ix, iy);
			BaseVerts[vi] = FVector(xy.X, xy.Y, zLocal);
			DynVerts[vi]  = BaseVerts[vi];
			Normals[vi]   = nLocal;
			UVs[vi]       = FVector2D((float)ix / GridX, (float)iy / GridY);
			Tangents[vi]  = FProcMeshTangent(1, 0, 0);
		}
	}

	for (int32 iy = 0; iy < GridY; ++iy)
	{
		for (int32 ix = 0; ix < GridX; ++ix)
		{
			const int32 v00 = VertIndex(ix, iy);
			const int32 v10 = VertIndex(ix + 1, iy);
			const int32 v01 = VertIndex(ix, iy + 1);
			const int32 v11 = VertIndex(ix + 1, iy + 1);
			Indices.Add(v00); Indices.Add(v11); Indices.Add(v10);
			Indices.Add(v00); Indices.Add(v01); Indices.Add(v11);
		}
	}

	// Draw 4 corner world-space spheres to locate the patch footprint
	{
		auto CornerWS = [&](float sx, float sy)
		{
			const FVector2D xy = FVector2D((sx-0.5f)*WorldSizeX, (sy-0.5f)*WorldSizeY);
			const FVector S = GetActorTransform().TransformPosition(FVector(xy.X, xy.Y,  200000.f));
			const FVector E = GetActorTransform().TransformPosition(FVector(xy.X, xy.Y, -200000.f));
			FHitResult H;
			if (GetWorld()->LineTraceSingleByObjectType(H, S, E,
				FCollisionObjectQueryParams(ECC_WorldStatic),
				FCollisionQueryParams(FName("SnowCorner"), true, this)) && H.bBlockingHit)
			{
				DrawDebugSphere(GetWorld(), H.ImpactPoint + FVector(0,0,15), 40.f, 12, FColor::Green, false, 10.f, 0, 2.f);
			}
		};
		CornerWS(0,0); CornerWS(1,0); CornerWS(0,1); CornerWS(1,1);
	}

	LegacyProcMesh->ClearAllMeshSections();
	LegacyProcMesh->CreateMeshSection_LinearColor(0, DynVerts, Indices, Normals, UVs, TArray<FLinearColor>(), Tangents, true);
	{
		UMaterialInterface* DefaultMat = UMaterial::GetDefaultMaterial(MD_Surface);
		if (DefaultMat)
		{
			LegacyProcMesh->SetMaterial(0, DefaultMat);
		}
		LegacyProcMesh->SetRelativeLocation(FVector::ZeroVector);
		LegacyProcMesh->SetRelativeRotation(FRotator::ZeroRotator);
		LegacyProcMesh->SetRelativeScale3D(FVector(1));
		LegacyProcMesh->MarkRenderStateDirty();
	}
	LegacyProcMesh->ContainsPhysicsTriMeshData(true);

	const FBoxSphereBounds B = LegacyProcMesh->Bounds;
	UE_LOG(LogTemp, Warning, TEXT("[Snow] Rebuild done. Hit=%d  Verts=%d  Tris=%d  Bounds Origin=%s Extent=%s"),
		hitCount, DynVerts.Num(), Indices.Num()/3, *B.Origin.ToString(), *B.BoxExtent.ToString());

	UE_LOG(LogTemp, Warning, TEXT("[Snow] Rebuild: Hit=%d of %d (%.1f%%). Actor at %s"),
		hitCount, (GridX+1)*(GridY+1), 100.f*hitCount/float((GridX+1)*(GridY+1)), *GetActorLocation().ToString());

	DrawDebugBox(W, B.Origin, B.BoxExtent, FColor::Yellow, false, 10.f, 0, 2.f);
}
#else
void ASnowSimulationActor::RebuildSnowMesh(){}
#endif


#if SNOW_WITH_LEGACY_PROCMESH
void ASnowSimulationActor::SnowMeshInfo()
{
	if (!LegacyProcMesh)
	{
		UE_LOG(LogTemp, Warning, TEXT("[Snow] LegacyProcMesh missing"));
		return;
	}
	const FBoxSphereBounds B = LegacyProcMesh->Bounds;
	UE_LOG(LogTemp, Log, TEXT("[Snow] Sections:%d  Verts:%d  Bounds O:%s Ext:%s  WorldLoc:%s"),
		LegacyProcMesh->GetNumSections(), DynVerts.Num(),
		*B.Origin.ToString(), *B.BoxExtent.ToString(),
		*GetActorLocation().ToString());

	// Visualize bounds
	DrawDebugBox(GetWorld(), B.Origin, B.BoxExtent, FColor::Cyan, false, 5.f, 0, 2.f);
}
#else
void ASnowSimulationActor::SnowMeshInfo()
{
	// Legacy path disabled; no-op
}
#endif


void ASnowSimulationActor::SnowMove(float X, float Y, float Z)
{
	SetActorLocation(FVector(X, Y, Z));
	UE_LOG(LogTemp, Warning, TEXT("[Snow] SnowMove -> %s"), *GetActorLocation().ToString());
	// Legacy procedural mesh removed
}


void ASnowSimulationActor::SnowSize(float SizeMeters)
{
	const float S = FMath::Max(1.f, SizeMeters) * 100.f; // m -> cm
	WorldSizeX = S;
	WorldSizeY = S;
	UE_LOG(LogTemp, Warning, TEXT("[Snow] SnowSize set to %.0f cm"), S);
	// Legacy procedural mesh removed
}

void ASnowSimulationActor::SnowHere(float SizeMeters)
{
	UWorld* W = GetWorld(); if(!W) { UE_LOG(LogTemp, Error, TEXT("[Snow] No world")); return; }
	APlayerController* PC = UGameplayStatics::GetPlayerController(W, 0); if(!PC) { UE_LOG(LogTemp, Error, TEXT("[Snow] No PC")); return; }

	FVector CamLoc; FRotator CamRot; PC->GetPlayerViewPoint(CamLoc, CamRot);
	const FVector Start = CamLoc;
	const FVector End   = CamLoc + CamRot.Vector() * 1000000.f;

	FHitResult Hit;
	FCollisionQueryParams Q(FName("SnowHere"), true, this);
	FCollisionObjectQueryParams Obj; Obj.AddObjectTypesToQuery(ECC_WorldStatic); Obj.AddObjectTypesToQuery(ECC_WorldDynamic);

	if (W->LineTraceSingleByObjectType(Hit, Start, End, Obj, Q) && Hit.bBlockingHit)
	{
		SetActorLocation(Hit.ImpactPoint + FVector(0,0,10));
		UE_LOG(LogTemp, Warning, TEXT("[Snow] SnowHere at %s on %s"),
			*Hit.ImpactPoint.ToString(), *GetNameSafe(Hit.GetActor()));
		SnowSize(SizeMeters);
	}
	else
	{
		UE_LOG(LogTemp, Error, TEXT("[Snow] SnowHere: camera ray found no ground; move closer and retry."));
	}
}

void ASnowSimulationActor::SnowLogCenter()
{
	const FVector2D xy(0,0);
	const FVector S = GetActorTransform().TransformPosition(FVector(xy.X, xy.Y,  200000.f));
	const FVector E = GetActorTransform().TransformPosition(FVector(xy.X, xy.Y, -200000.f));
	UE_LOG(LogTemp, Warning, TEXT("[Snow] Center trace Start=%s End=%s Loc=%s Size=%.0f x %.0f (cm)"),
		*S.ToString(), *E.ToString(), *GetActorLocation().ToString(), WorldSizeX, WorldSizeY);
}

