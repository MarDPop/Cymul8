#pragma once

#include "../Vehicle.h"
#include "Components.h"

#include <vector>

class Vehicle_Component : public virtual Vehicle
{
	std::vector<Component> _components;

public:

};