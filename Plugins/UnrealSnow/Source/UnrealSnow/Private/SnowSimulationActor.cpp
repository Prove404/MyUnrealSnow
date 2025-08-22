#include "Simulation/SnowSimulationActor.h"
#include "UnrealSnowLog.h"

#include "SnowRTUtils.h"
#include "Engine/TextureRenderTarget2D.h"
#include "Weather/StochasticWeatherProvider.h"
#include "UObject/SoftObjectPtr.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Simulation/ISnowRedistributor.h"
#include "Simulation/RidgeLeewardRedistributor.h"

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

ASnowSimulationActor::ASnowSimulationActor()
{
	PrimaryActorTick.bCanEverTick = true;
	PrimaryActorTick.bStartWithTickEnabled = true;

	DebugPlane = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("DebugPlane"));
	SetRootComponent(DebugPlane);
	static ConstructorHelpers::FObjectFinder<UStaticMesh> PlaneMesh(TEXT("/Engine/BasicShapes/Plane"));
	if (PlaneMesh.Succeeded())
	{
		DebugPlane->SetStaticMesh(PlaneMesh.Object);
		DebugPlane->SetRelativeScale3D(FVector(10.f, 10.f, 1.f));
	}
}

void ASnowSimulationActor::BeginPlay()
{
	Super::BeginPlay();
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

	if (DebugMaterial && SnowDepthRT)
	{
		DebugMID = UMaterialInstanceDynamic::Create(DebugMaterial, this);
		DebugMID->SetTextureParameterValue(DebugTextureParam, SnowDepthRT);
		if (DebugPlane)
		{
			DebugPlane->SetMaterial(0, DebugMID);
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

	static float OverlayAccum = 0.f;
	OverlayAccum += DeltaSeconds;
	if (OverlayAccum > 1.f)
	{
		OverlayAccum = 0.f;
		float MinV = FLT_MAX, MaxV = -FLT_MAX;
		for (float v : Depth01) { MinV = FMath::Min(MinV, v); MaxV = FMath::Max(MaxV, v); }
		if (GEngine)
		{
			const FString Msg = FString::Printf(TEXT("Depth01 min=%.4f max=%.4f SWEmean=%.3f"), MinV, MaxV, DebugMeanSWE);
			GEngine->AddOnScreenDebugMessage((uint64)this, 1.1f, FColor::Cyan, Msg);
		}
	}

	static float Acc = 0.f;
	Acc += DeltaSeconds;
	if (Acc > 1.f)
	{
		Acc = 0.f;
		double Sum = 0.0;
		for (float v : SWE) Sum += v;
		DebugMeanSWE = (NumCells > 0) ? float(Sum / NumCells) : 0.f;
		UE_LOG(LogUnrealSnow, Display, TEXT("t=%s Precip=%.3f mm/h SWEmean=%.3f"),
			*CurrentSimTime.ToString(), Precip_mmph, DebugMeanSWE);
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


