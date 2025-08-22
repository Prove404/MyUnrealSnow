#pragma once

#include "CoreMinimal.h"
#include "Engine/TextureRenderTarget2D.h"

/** Internal RT utilities for UnrealSnow (CPU side). */
class USnowRTUtils final
{
public:
	/** Create or reinitialize a float render target with square Size. */
	static bool InitFloatRT(UTextureRenderTarget2D*& RT, int32 Size, EPixelFormat PF = PF_R32_FLOAT);

	/** Write a single-channel float array (width*height) to the render target. */
	static void WriteArrayToRT(UTextureRenderTarget2D* RT, const TArray<float>& Data);
};


