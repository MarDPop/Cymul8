#pragma once

#include "../../physics/Body.h"
#include "../../physics/Action.h"
#include "Resource.h"


class Component 
{
	Eigen::Vector3d _body_location;

	double _mass;

	MomentOfInertia<MOMENT_CONSTANTS::FULL> _inertia;

	BodyAction _action;

public:

	virtual void update(const Vehicle_Component& vehicle, double time) = 0;

	const BodyAction& get_action() const
	{
		return _action;
	}

	const BodyAction& get_action() const
	{
		return _action;
	}

};