#pragma once

#include "Modules/ModuleManager.h"

class FUnrealSnowModule : public IModuleInterface
{
public:
	virtual void StartupModule() override;
	virtual void ShutdownModule() override;
};


