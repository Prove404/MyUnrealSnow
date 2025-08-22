#include "SnowRTUtils.h"

#include "Engine/TextureRenderTarget2D.h"
#include "RenderUtils.h"
#include "RHI.h"
#include "RHICommandList.h"
#include "RHIResources.h"
#include "RenderResource.h"

bool USnowRTUtils::InitFloatRT(UTextureRenderTarget2D*& RT, int32 Size, EPixelFormat PF)
{
	if (Size <= 0)
	{
		return false;
	}

	if (!RT)
	{
		RT = NewObject<UTextureRenderTarget2D>(GetTransientPackage());
		RT->AddToRoot();
	}

	RT->ClearColor = FLinearColor::Black;
	RT->bAutoGenerateMips = false;
	RT->AddressX = TA_Clamp;
	RT->AddressY = TA_Clamp;
	RT->InitCustomFormat(Size, Size, PF, false);
	RT->UpdateResourceImmediate(true);

	// One-time debug: write a horizontal gradient on first init for PF_R32_FLOAT targets
	static bool bDidDebugGradient = false;
	if (!bDidDebugGradient && PF == PF_R32_FLOAT)
	{
		bDidDebugGradient = true;
		TArray<float> Gradient;
		Gradient.SetNumUninitialized(Size * Size);
		for (int32 y = 0; y < Size; ++y)
		{
			for (int32 x = 0; x < Size; ++x)
			{
				const float v = (Size > 1) ? (static_cast<float>(x) / static_cast<float>(Size - 1)) : 0.0f;
				Gradient[y * Size + x] = v;
			}
		}
		USnowRTUtils::WriteArrayToRT(RT, Gradient);
	}
	return true;
}

void USnowRTUtils::WriteArrayToRT(UTextureRenderTarget2D* RT, const TArray<float>& Data)
{
	if (!RT)
	{
		return;
	}

	const int32 Width = RT->SizeX;
	const int32 Height = RT->SizeY;
	if (Width <= 0 || Height <= 0)
	{
		return;
	}

	const int32 Num = Width * Height;
	if (Data.Num() != Num)
	{
		return;
	}

	FTextureRenderTargetResource* Resource = RT->GameThread_GetRenderTargetResource();

	// Copy data so we can safely pass it to render thread
	TArray<float> DataCopy = Data;
	const EPixelFormat PF = RT->GetFormat();

	ENQUEUE_RENDER_COMMAND(WriteSnowRT)
	([
		Resource,
		DataCopy = MoveTemp(DataCopy),
		Width,
		Height,
		PF
	](FRHICommandListImmediate& RHICmdList)
	{
		FRHITexture* TextureRHI = Resource->GetRenderTargetTexture();
		if (TextureRHI == nullptr)
		{
			return;
		}

		if (PF != PF_R32_FLOAT)
		{
			// Only PF_R32_FLOAT is supported by this writer
			return;
		}

		const uint32 SrcPitchBytes = static_cast<uint32>(Width * sizeof(float));
		const uint8* SrcBytes = reinterpret_cast<const uint8*>(DataCopy.GetData());
		FUpdateTextureRegion2D Region(0, 0, 0, 0, Width, Height);
		RHICmdList.UpdateTexture2D(TextureRHI, 0, Region, SrcPitchBytes, SrcBytes);
	});
}


