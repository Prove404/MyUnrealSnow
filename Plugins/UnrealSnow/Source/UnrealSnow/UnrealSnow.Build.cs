using UnrealBuildTool;

public class UnrealSnow : ModuleRules
{
	public UnrealSnow(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;


		PublicDependencyModuleNames.AddRange(new string[]
		{
			"Core", "CoreUObject", "Engine", "InputCore", "Landscape", "DeveloperSettings"
		});


		PrivateDependencyModuleNames.AddRange(new string[]
		{
			"Projects", "RenderCore", "RHI", "VirtualHeightfieldMesh"
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

