#include "CoreMinimal.h"
#include "Engine/TextureRenderTarget2D.h"
#include "GlobalShader.h"
#include "ShaderParameterStruct.h"
#include "RenderGraphBuilder.h"
#include "RenderGraphUtils.h"
#include "RHICommandList.h"
#include "RHIResources.h"

// Compute shader wrapper
class FSnowCS : public FGlobalShader
{
public:
    DECLARE_GLOBAL_SHADER(FSnowCS);
    SHADER_USE_PARAMETER_STRUCT(FSnowCS, FGlobalShader);

    BEGIN_SHADER_PARAMETER_STRUCT(FParameters, )
        SHADER_PARAMETER(FIntPoint, OutputSize)
        SHADER_PARAMETER(float, InvMaxDepth)
        SHADER_PARAMETER_RDG_TEXTURE_SRV(Texture2D<float>, InSWE)
        SHADER_PARAMETER_SAMPLER(SamplerState, InSWESampler)
        SHADER_PARAMETER_RDG_TEXTURE_UAV(RWTexture2D<float>, OutSnowDepth)
    END_SHADER_PARAMETER_STRUCT()

    static bool ShouldCompilePermutation(const FGlobalShaderPermutationParameters& Parameters)
    {
        return true;
    }
};

IMPLEMENT_GLOBAL_SHADER(FSnowCS, "/UnrealSnow/SnowCS.usf", "MainCS", SF_Compute);

// Helper to dispatch the compute shader
static void DispatchSnowCS(FRHICommandListImmediate& RHICmdList, FRDGBuilder& GraphBuilder,
    FRDGTextureRef InSWE, FRDGTextureRef TargetTexture, float InvMaxDepth)
{
    FSnowCS::FParameters* PassParameters = GraphBuilder.AllocParameters<FSnowCS::FParameters>();
    PassParameters->OutputSize = FIntPoint(TargetTexture->Desc.Extent.X, TargetTexture->Desc.Extent.Y);
    PassParameters->InvMaxDepth = InvMaxDepth;
    PassParameters->InSWE = GraphBuilder.CreateSRV(FRDGTextureSRVDesc::Create(InSWE));
    PassParameters->InSWESampler = TStaticSamplerState<SF_Bilinear, AM_Clamp, AM_Clamp, AM_Clamp>::GetRHI();
    PassParameters->OutSnowDepth = GraphBuilder.CreateUAV(FRDGTextureUAVDesc(TargetTexture));

    TShaderMapRef<FSnowCS> ComputeShader(GetGlobalShaderMap(GMaxRHIFeatureLevel));

    const FIntVector GroupCounts(
        FMath::DivideAndRoundUp(PassParameters->OutputSize.X, 8),
        FMath::DivideAndRoundUp(PassParameters->OutputSize.Y, 8),
        1);

    FComputeShaderUtils::AddPass(
        GraphBuilder,
        RDG_EVENT_NAME("UnrealSnow::SnowCS"),
        ComputeShader,
        PassParameters,
        GroupCounts);
}

// Public entry point callable from game thread to run the compute shader on a render target
void RunSnowCS_RenderTarget(class UTextureRenderTarget2D* RT, class UTextureRenderTarget2D* SWETexture, float InvMaxDepth)
{
    if (!RT)
    {
        return;
    }
    FTextureRenderTargetResource* Resource = RT->GameThread_GetRenderTargetResource();
    FTextureRenderTargetResource* SWEResource = SWETexture ? SWETexture->GameThread_GetRenderTargetResource() : nullptr;
    ENQUEUE_RENDER_COMMAND(RunSnowCS)
    ([Resource, SWEResource, InvMaxDepth](FRHICommandListImmediate& RHICmdList)
    {
        FRDGBuilder GraphBuilder(RHICmdList);
        FRDGTextureRef RDGTex = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(Resource->GetRenderTargetTexture(), TEXT("SnowDepthRT")));
        FRDGTextureRef RDGSWE = GraphBuilder.RegisterExternalTexture(CreateRenderTarget(SWEResource->GetRenderTargetTexture(), TEXT("SWETexture")));
        DispatchSnowCS(RHICmdList, GraphBuilder, RDGSWE, RDGTex, InvMaxDepth);
        GraphBuilder.Execute();
    });
}


