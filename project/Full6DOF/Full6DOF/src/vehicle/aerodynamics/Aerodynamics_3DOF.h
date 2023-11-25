#include "Aerodynamics.h"

class Aerodynamics_3DOF_None : public virtual Aerodynamics<Eigen::Vector3d>
{

};

class Aerodynamics_3DOF_Const_Drag_Coefficient : public virtual Aerodynamics<Eigen::Vector3d>
{
    const double _CDA;

public:

    Aerodynamics_3DOF_Const_Drag_Coefficient(double __CD,
        double __ref_area) :
        _CDA(-__CD*__ref_area) {}

    void update(const AeroData& aeroData,
        double t) override
    {
        this->_action = (_CDA*aeroData.dynamic_pressure)*aeroData.air_velocity_ecef_unit;
    }

};