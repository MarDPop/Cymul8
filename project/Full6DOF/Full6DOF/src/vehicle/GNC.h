#pragma once

#include "../physics/Body.h"
#include <memory>

struct GNC
{
    virtual unsigned get_control_states() const = 0;

    virtual void update(const double* control_states, const double time, double* dx_control) = 0;
};

template<class S>
class Navigation
{
protected:

    S _estimated_state;

public:

    virtual void update(const double time) {}

    const S& get_estimated_state() const
    {
        return _estimated_state;
    }
};

template<class S>
class Guidance
{
protected:

    S _desired_state;

public:

    virtual void update(const S& estimated,
        const double time) {}

    const S& get_desired_state() const
    {
        return _desired_state;
    }
};

template<class S>
class Control
{

public:

    const unsigned N_CONTROL_STATES;

    Control(unsigned __nControlStates = 0u) :
        N_CONTROL_STATES(__nControlStates) {}

    virtual void update(const S& desired_state,
        const S& estimated_state,
        const double* control_states,
        const double time,
        double* dx_control) {}
};

template<class S, class G, class N, class C>
class GuidanceNavigationControl_T : public virtual GNC
{
    static_assert(std::is_base_of<Guidance<S>, G>::value, "G not derived from Guidance");
    static_assert(std::is_base_of<Navigation<S>, N>::value, "N not derived from Navigation");
    static_assert(std::is_base_of<Control<S>, C>::value, "C not derived from Control");

    double _delay = 0;

    G _guidance; // check to see if this devirtualizes update()

    N _navigation;

    C _control;

public:

    const G& get_guidance() const
    {
        *_guidance;
    }

    const N& get_navigation() const
    {
        *_navigation;
    }

    const C& get_control() const
    {
        *_control;
    }

    void set_delay(double delay)
    {
        _delay = delay;
    }

    unsigned get_control_states() const override
    {
        return C.N_CONTROL_STATES;
    }

    void update(const double* control_states, const double time, double* dx_control) override
    {
        double time_delay = time - _delay;

        _navigation.update(time_delay);

        _guidance.update(   time_delay, 
                            _navigation.get_estimated_state());

        _control.update(    _navigation.get_estimated_state(),
                            _guidance.get_desired_state(),
                            control_states, 
                            time_delay, 
                            dx_control);
    }    
};

template<class S>
class GuidanceNavigationControl : public virtual GNC
{

    double _delay = 0;

    std::unique_ptr<Guidance<S>> _guidance = std::make_unique<Guidance<S>>();

    std::unique_ptr<Navigation<S>> _navigation = std::make_unique<Guidance<S>>();

    std::unique_ptr<Control<S>> _control = std::make_unique<Guidance<S>>();

public:

    const Guidance<S>& get_guidance() const
    {
        *_guidance;
    }

    void set_guidance(std::unique_ptr<Guidance<S>> guidance)
    {
        _guidance = std::move(guidance);
    }

    const Navigation<S>& get_navigation() const
    {
        *_navigation;
    }

    void set_navigation(std::unique_ptr<Navigation<S>> navigation)
    {
        _navigation = std::move(navigation);
    }

    const Control<S>& get_control() const
    {
        *_control;
    }
    
    void set_control(std::unique_ptr<Control<S>> control)
    {
        _control = std::move(control);
    }

    void set_delay(double delay)
    {
        _delay = delay;
    }

    unsigned get_control_states() const override
    {
        return _control->N_CONTROL_STATES;
    }

    void update(const double* control_states, const double time, double* dx_control) override
    {
        double time_delay = time - _delay;

        _navigation->update(time_delay);

        _guidance->update(_navigation->get_estimated_state(),
            time_delay);

        _control->update(_navigation->get_estimated_state(),
            _guidance->get_desired_state(),
            control_states,
            time_delay,
            dx_control);
    }
};
