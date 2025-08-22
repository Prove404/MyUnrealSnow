#include "Simulation/SWROutExporter.h"

#include "Engine/TextureRenderTarget2D.h"
#include "Async/Async.h"
#include "ImageUtils.h"
#include "Containers/Array.h"
#include "Misc/FileHelper.h"
#include "HAL/PlatformFilemanager.h"
#include "RHI.h"

void USWROutExporter::ExportNow()
{
	if (!TargetRT) return;
	const FString Dir = OutputDir.Path;
	IPlatformFile& PF = FPlatformFileManager::Get().GetPlatformFile();
	if (!PF.DirectoryExists(*Dir))
	{
		PF.CreateDirectoryTree(*Dir);
	}
	const FString Base = FPaths::Combine(Dir, TEXT("SWRout_Manual"));
	EnqueueWrite(Base, TargetRT);
}

void USWROutExporter::TickExport(int32 StepIndex, const FDateTime& SimTime)
{
	if (!TargetRT || ExportEveryNSteps <= 0) return;
	if ((StepIndex % ExportEveryNSteps) != 0) return;
	const FString Dir = OutputDir.Path;
	IPlatformFile& PF = FPlatformFileManager::Get().GetPlatformFile();
	if (!PF.DirectoryExists(*Dir))
	{
		PF.CreateDirectoryTree(*Dir);
	}
	const FString Stamp = SimTime.ToString(TEXT("yyyyMMdd_HHmm"));
	const FString Base = FPaths::Combine(Dir, FString::Printf(TEXT("SWRout_%s"), *Stamp));
	EnqueueWrite(Base, TargetRT);
}

void USWROutExporter::EnqueueWrite(const FString& BasePath, UTextureRenderTarget2D* RT)
{
	FTextureRenderTargetResource* Res = RT->GameThread_GetRenderTargetResource();
	// Readback on render thread, then write files on background thread
	TArray<FFloat16Color> SurfaceData;
	FIntRect Rect(0, 0, RT->SizeX, RT->SizeY);
	FReadSurfaceDataFlags Flags(RCM_UNorm);
	Res->ReadFloat16Pixels(SurfaceData, Flags, Rect);

	const int32 W = RT->SizeX;
	const int32 H = RT->SizeY;

	Async(EAsyncExecution::ThreadPool, [SurfaceData = MoveTemp(SurfaceData), W, H, BasePath, this]()
	{
		// Convert to grayscale colors for PNG compression
		TArray<FColor> ColorData;
		ColorData.Reserve(W * H);
		for (int32 i = 0; i < W * H; ++i)
		{
			float g = SurfaceData[i].R.GetFloat();
			g = FMath::Clamp(g, 0.0f, 1.0f);
			uint8 b = (uint8)FMath::RoundToInt(g * 255.0f);
			ColorData.Add(FColor(b, b, b, 255));
		}

		// Write PNG
		TArrayView64<const FColor> ColorView(ColorData.GetData(), ColorData.Num());
		TArray64<uint8> PNG;
		FImageUtils::PNGCompressImageArray(W, H, ColorView, PNG);
		FFileHelper::SaveArrayToFile(PNG, *(BasePath + TEXT(".png")));

		// Write world file (.pgw) with pixel size and origin
		const FString PGW = FString::Printf(TEXT("%f\n0\n0\n-%f\n%f\n%f\n"),
			PixelSizeMeters.X, PixelSizeMeters.Y, OriginMeters.X, OriginMeters.Y + H * PixelSizeMeters.Y);
		FFileHelper::SaveStringToFile(PGW, *(BasePath + TEXT(".pgw")));

		// Write projection (.prj)
		const FString PRJ = TEXT("GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433]]");
		FFileHelper::SaveStringToFile(PRJ, *(BasePath + TEXT(".prj")));
	});
}


