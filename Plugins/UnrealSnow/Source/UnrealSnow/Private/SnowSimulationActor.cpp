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

// Forward decl from SnowCS.cpp
void RunSnowCS_RenderTarget(class UTextureRenderTarget2D* RT, class UTextureRenderTarget2D* SWETexture, float InvMaxDepth);

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

ASnowSimulationActor::ASnowSimulationActor()
{
	PrimaryActorTick.bCanEverTick = true;

	if (!RootComponent)
	{
		USceneComponent* Root = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
		SetRootComponent(Root);
	}
	SnowMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("SnowMesh"));
	SnowMesh->SetupAttachment(RootComponent);
	SnowMesh->SetMobility(EComponentMobility::Movable);
	SnowMesh->bUseAsyncCooking = true;
	SnowMesh->SetVisibility(true, true);
	SnowMesh->SetHiddenInGame(false);
	SnowMesh->bRenderInMainPass = true;
	SnowMesh->SetCastShadow(false);
	SnowMesh->SetBoundsScale(3.0f);

	// Default to hidden/disabled (VHM primary path)
	if (SnowMesh)
	{
		SnowMesh->SetHiddenInGame(true);
		SnowMesh->SetVisibility(false, true);
		SnowMesh->SetCollisionEnabled(ECollisionEnabled::NoCollision);
	}
}

void ASnowSimulationActor::OnConstruction(const FTransform& Xform)
{
	Super::OnConstruction(Xform);
	RebuildSnowMesh();
}

void ASnowSimulationActor::BeginPlay()
{
	Super::BeginPlay();

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
						SnowWriterMID->SetTextureParameterValue(SnowDepthParam, SnowDepthTex);
					}
				}
			}
		}
	}

	// Build the procedural snow mesh surface
	RebuildSnowMesh();
	if (bHideLegacySnowPlane)
	{
		// If your actor has an older plane component, hide it here, e.g.:
		// if (SnowPlaneComponent) SnowPlaneComponent->SetVisibility(false, true);
		// (Leave as-is if you don’t have such a component.)
	}
	EnsureSWEStorage();

	// Hide legacy primitives except SnowMesh so only the procedural surface is visible
	{
		TArray<UPrimitiveComponent*> PrimComps;
		GetComponents<UPrimitiveComponent>(PrimComps);
		for (UPrimitiveComponent* P : PrimComps)
		{
			if (P && P != SnowMesh)
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

	LogSWEStats(SWE_Grid, GridX, GridY, SnowDensity, DisplacementScale);

	if (!HasValidMeshBuffers())
	{
		UE_LOG(LogTemp, Warning, TEXT("[Snow] Mesh buffers not ready; calling RebuildSnowMesh()"));
		RebuildSnowMesh();
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
				UE_LOG(LogTemp, Warning, TEXT("[Snow] DynVerts size %d != expected %d; rebuilding"),
					DynVerts.Num(), (GridX+1)*(GridY+1));
				RebuildSnowMesh();
				if (DynVerts.Num() != (GridX+1)*(GridY+1)) return;
			}
			SnowMesh->UpdateMeshSection_LinearColor(0, DynVerts, Normals, UVs, TArray<FLinearColor>(), Tangents);
		}
	}
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
	SnowDepthTex->Filter = TF_Bilinear;
	SnowDepthTex->AddressX = TA_Clamp; SnowDepthTex->AddressY = TA_Clamp;
	SnowDepthTex->UpdateResource();

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
	static TArray<uint16> Buf;
	Buf.SetNumUninitialized(W * H);

	const float rhoWater = 1000.f;
	for (int32 y = 0; y < H; ++y)
	{
		const float ty = (float)y / (H - 1);
		for (int32 x = 0; x < W; ++x)
		{
			const float tx = (float)x / (W - 1);
			const int32 gx = FMath::Clamp(FMath::RoundToInt(tx * (GridX - 1)), 0, GridX - 1);
			const int32 gy = FMath::Clamp(FMath::RoundToInt(ty * (GridY - 1)), 0, GridY - 1);

			const float SWE_m = GetSWE_Clamped(gx, gy);
			const float depth_m = SWE_m * (rhoWater / FMath::Max(1.f, SnowDensity));

			Buf[y * W + x] = FFloat16(depth_m).Encoded;
		}
	}

	float maxDepth = 0.f;
	for (int32 i=0; i < W*H; ++i) { maxDepth = FMath::Max(maxDepth, FFloat16(Buf[i]).GetFloat()); }
	UE_LOG(LogTemp, Warning, TEXT("[Snow] DepthTex max=%.3fm (pre-scale)"), maxDepth);

	const uint32 SrcPitch = W * sizeof(uint16);
	const uint32 SrcBpp   = sizeof(uint16);
	SnowDepthTex->UpdateTextureRegions(0, 1, &SnowDepthRegion, SrcPitch, SrcBpp,
		reinterpret_cast<uint8*>(Buf.GetData()));
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

void ASnowSimulationActor::RebuildSnowMesh()
{
	UE_LOG(LogTemp, Warning, TEXT("[Snow] RebuildSnowMesh() start Grid=%dx%d World=%.0fx%.0f cm"),
		GridX, GridY, WorldSizeX, WorldSizeY);

	if (GridX < 2 || GridY < 2 || WorldSizeX <= 0 || WorldSizeY <= 0)
	{
		UE_LOG(LogTemp, Error, TEXT("[Snow] Invalid grid/size. GridX,Y>=2 and WorldSize>0 required."));
		return;
	}

	if (!SnowMesh)
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
	if (SnowMesh) { Q.AddIgnoredComponent(SnowMesh); }
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

	SnowMesh->ClearAllMeshSections();
	SnowMesh->CreateMeshSection_LinearColor(0, DynVerts, Indices, Normals, UVs, TArray<FLinearColor>(), Tangents, true);
	{
		UMaterialInterface* DefaultMat = UMaterial::GetDefaultMaterial(MD_Surface);
		if (DefaultMat)
		{
			SnowMesh->SetMaterial(0, DefaultMat);
		}
		SnowMesh->SetRelativeLocation(FVector::ZeroVector);
		SnowMesh->SetRelativeRotation(FRotator::ZeroRotator);
		SnowMesh->SetRelativeScale3D(FVector(1));
		SnowMesh->MarkRenderStateDirty();
	}
	SnowMesh->ContainsPhysicsTriMeshData(true);

	const FBoxSphereBounds B = SnowMesh->Bounds;
	UE_LOG(LogTemp, Warning, TEXT("[Snow] Rebuild done. Hit=%d  Verts=%d  Tris=%d  Bounds Origin=%s Extent=%s"),
		hitCount, DynVerts.Num(), Indices.Num()/3, *B.Origin.ToString(), *B.BoxExtent.ToString());

	UE_LOG(LogTemp, Warning, TEXT("[Snow] Rebuild: Hit=%d of %d (%.1f%%). Actor at %s"),
		hitCount, (GridX+1)*(GridY+1), 100.f*hitCount/float((GridX+1)*(GridY+1)), *GetActorLocation().ToString());

	DrawDebugBox(W, B.Origin, B.BoxExtent, FColor::Yellow, false, 10.f, 0, 2.f);
}


void ASnowSimulationActor::SnowMeshInfo()
{
	if (!SnowMesh)
	{
		UE_LOG(LogTemp, Warning, TEXT("[Snow] SnowMesh missing"));
		return;
	}
	const FBoxSphereBounds B = SnowMesh->Bounds;
	UE_LOG(LogTemp, Log, TEXT("[Snow] Sections:%d  Verts:%d  Bounds O:%s Ext:%s  WorldLoc:%s"),
		SnowMesh->GetNumSections(), DynVerts.Num(),
		*B.Origin.ToString(), *B.BoxExtent.ToString(),
		*GetActorLocation().ToString());

	// Visualize bounds
	DrawDebugBox(GetWorld(), B.Origin, B.BoxExtent, FColor::Cyan, false, 5.f, 0, 2.f);
}


void ASnowSimulationActor::SnowMove(float X, float Y, float Z)
{
	SetActorLocation(FVector(X, Y, Z));
	UE_LOG(LogTemp, Warning, TEXT("[Snow] SnowMove -> %s"), *GetActorLocation().ToString());
	RebuildSnowMesh();
}


void ASnowSimulationActor::SnowSize(float SizeMeters)
{
	const float S = FMath::Max(1.f, SizeMeters) * 100.f; // m -> cm
	WorldSizeX = S;
	WorldSizeY = S;
	UE_LOG(LogTemp, Warning, TEXT("[Snow] SnowSize set to %.0f cm"), S);
	RebuildSnowMesh();
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

