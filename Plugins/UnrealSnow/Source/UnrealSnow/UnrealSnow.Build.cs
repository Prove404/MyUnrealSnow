using UnrealBuildTool;

public class UnrealSnow : ModuleRules
{
	public UnrealSnow(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;

		PublicDependencyModuleNames.AddRange(new string[]
		{
			"Core", "CoreUObject", "Engine", "InputCore", "RenderCore", "RHI", "Projects", "DeveloperSettings", "ProceduralMeshComponent"
		});

		PrivateDependencyModuleNames.AddRange(new string[]
		{
			"Slate", "SlateCore", "Renderer"
		});

		PublicIncludePaths.AddRange(new[] {
			System.IO.Path.Combine(ModuleDirectory, "Public"),
			System.IO.Path.Combine(ModuleDirectory, "Public/Simulation"),
			System.IO.Path.Combine(ModuleDirectory, "Public/Weather")
		});

		if (Target.bBuildEditor)
		{
			PrivateDependencyModuleNames.AddRange(new string[]
			{
				"UnrealEd", "AssetTools", "AssetRegistry"
			});
		}
	}
}

