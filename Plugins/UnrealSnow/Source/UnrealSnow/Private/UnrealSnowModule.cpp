#include "UnrealSnowModule.h"
#include "Interfaces/IPluginManager.h"
#include "Misc/Paths.h"
#include "ShaderCore.h"
// Ensure actor class is linked in and registered
#include "Simulation/SnowSimulationActor.h"
#include "UnrealSnowLog.h"

 

void FUnrealSnowModule::StartupModule()
{
	const TSharedPtr<IPlugin> Plugin = IPluginManager::Get().FindPlugin(TEXT("UnrealSnow"));
	if (Plugin.IsValid())
	{
		const FString PluginBaseDir = Plugin->GetBaseDir();
		const FString ShaderDir = FPaths::Combine(PluginBaseDir, TEXT("Shaders"));
		AddShaderSourceDirectoryMapping(TEXT("/UnrealSnow"), ShaderDir);
		UE_LOG(LogTemp, Log, TEXT("UnrealSnow shader dir mapped: %s"), *ShaderDir);
	}

	// NOTE(UE5): removed runtime material graph mutation; use a preauthored material instance in editor/actor.
}

void FUnrealSnowModule::ShutdownModule()
{
}

IMPLEMENT_MODULE(FUnrealSnowModule, UnrealSnow)


