using UnrealBuildTool;

public class UnrealSnow : ModuleRules
{
	public UnrealSnow(ReadOnlyTargetRules Target) : base(Target)
	{
		PCHUsage = PCHUsageMode.UseExplicitOrSharedPCHs;


		PublicDependencyModuleNames.AddRange(new string[]
		{
			"Core",
			"CoreUObject",
			"Engine",
			"RenderCore",
			"RHI",
			"Projects",
			"VirtualHeightfieldMesh",   // AVirtualHeightfieldMesh + UVirtualHeightfieldMeshComponent
			"Landscape",
			"DeveloperSettings"
		});


		PrivateDependencyModuleNames.AddRange(new string[]
		{
		});

		PublicIncludePaths.AddRange(new[] {
			System.IO.Path.Combine(ModuleDirectory, "Public"),
			System.IO.Path.Combine(ModuleDirectory, "Public/Simulation"),
			System.IO.Path.Combine(ModuleDirectory, "Public/Weather"),
			System.IO.Path.Combine(ModuleDirectory, "Public/Redistribution")
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

