#pragma once

#include <string>

class Resource
{

    std::string _name;

    double _mass_per_unit;

    double _cost_per_unit;

public:

    inline Resource(std::string name, 
                    double mass_per_unit = 0, 
                    double cost_per_unit = 0) :
                    _name(name),
                    _mass_per_unit(mass_per_unit),
                    _cost_per_unit(cost_per_unit) {}

    inline std::string get_name() const
    {
        return this->_name;
    }

    inline double get_mass_per_unit() const
    {
        return this->_mass_per_unit;
    }

    inline double get_cost_per_unit() const
    {
        return this->_cost_per_unit;
    }


};


class ResourceQuantity
{
protected:

    const Resource& _resource;

    double _resource_quanitity_full;

    double _resource_quantity;
    
public:

    ResourceQuantity(const Resource& _resource, double resource_quanitity_full);
    virtual ~ResourceQuantity();


};