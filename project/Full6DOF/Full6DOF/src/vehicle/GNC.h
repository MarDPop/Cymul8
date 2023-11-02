#pragma once

template<class V>
class Navigation
{
    const V& _vehicle;

    double* const _estimated_state;

public:

    Navigation(const V& vehicle) :
        _vehicle(vehicle)
    {}

    virtual ~Navigation
    {
        delete[] _estimated_state;
    }

    virtual void estimate_state(const double time)
    {
        _vehicle.get_state(_estimated_state);
    }

    const double* get_estimated_state() const
    {
        return _estimated_state;
    }
};

template<class V>
class Guidance
{
    const V& _vehicle;

    double* const _desired_state;

public:

    const unsigned nVehicleStates;

    Guidance(const V& vehicle) :
        _vehicle(vehicle)
    {}

    virtual ~Guidance
    {
        delete[] _desired_state;
    }

    virtual void set_estimated_state(const double* estimated, 
                                    const double time)
    {
        memcpy(_desired_state, esimated, _nVehicleStates * sizeof(double));
    }

    const double* get_desired_state() const
    {
        return _desired_state;
    }
};

class Control
{

public:

    const unsigned nControlStates;

    Control(unsigned __nControlStates = 0u) :
        nControlStates(__nControlStates) {}

    virtual void set_desired_state(const double* desired_state, 
                                    const double* control_states,
                                    const double time,
                                    double* dx_control)
    {

    }
};

template<class V, class G, class N, class C>
class GuidanceNavigationControl
{
    static_assert(std::is_base_of<Guidance, G>::value, "G not derived from BaseClass");
    static_assert(std::is_base_of<Navigation, N>::value, "N not derived from BaseClass");
    static_assert(std::is_base_of<Control, C>::value, "C not derived from BaseClass");

    double _delay = 0;

public:

    G guidance;

    N navigation;

    C control;

    GuidanceNavigationControl(const V& __vehicle) :
                guidance(__vehicle),
                navigation(__vehicle)
    {}

    void set_delay(double delay)
    {
        _delay = delay;
    }

    void update(const double* control_states, const double time, double* dx_control)
    {
        double time_delay = time - _delay;
        _navigation.estimate_state(time_delay);
        _guidance.set_estimated_state(time_delay, _navigation.get_estimated_state());
        _control.estimate_state(_guidance.get_desired_state(), control_states, time_delay, dx_control);
    }

    
};
