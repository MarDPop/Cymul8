#pragma once

#include "../physics/Body.h"
#include <memory>

template<class B>
struct GNC
{
    static_assert(std::is_base_of<Body_Base, B>::value, "GNC needs to be a base of Body");

    virtual unsigned get_control_states() = 0;

    virtual void update(const double* control_states, const double time, double* dx_control) = 0;
};

template<class B>
class Navigation
{
protected:

    B _estimated_state;

public:

    virtual void update(const double time) {}

    const B& get_estimated_state() const
    {
        return _estimated_state;
    }
};

template<class B>
class Guidance
{
protected:

    B _desired_state;

public:

    virtual void update(const B& estimated,
        const double time) {}

    const B& get_desired_state() const
    {
        return _desired_state;
    }
};

template<class B>
class Control
{

public:

    const unsigned N_CONTROL_STATES;

    Control(unsigned __nControlStates = 0u) :
        N_CONTROL_STATES(__nControlStates) {}

    virtual void update(const B& desired_state,
        const B& estimated_state,
        const double* control_states,
        const double time,
        double* dx_control) {}
};

template<class B, class G, class N, class C>
class GuidanceNavigationControl_T : public virtual GNC<B>
{
    static_assert(std::is_base_of<Guidance<B>, G>::value, "G not derived from Guidance");
    static_assert(std::is_base_of<Navigation<B>, N>::value, "N not derived from Navigation");
    static_assert(std::is_base_of<Control<B>, C>::value, "C not derived from Control");

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

template<class B>
class GuidanceNavigationControl : public virtual GNC<B>
{

    double _delay = 0;

    std::unique_ptr<Guidance<B>> _guidance = std::make_unique<Guidance<B>>();

    std::unique_ptr<Navigation<B>> _navigation = std::make_unique<Guidance<B>>();

    std::unique_ptr<Control<B>> _control = std::make_unique<Guidance<B>>();

public:

    const Guidance<B>& get_guidance() const
    {
        *_guidance;
    }

    void set_guidance(std::unique_ptr<Guidance<B>> guidance)
    {
        _guidance = std::move(guidance);
    }

    const Navigation<B>& get_navigation() const
    {
        *_navigation;
    }

    void set_navigation(std::unique_ptr<Navigation<B>> navigation)
    {
        _navigation = std::move(navigation);
    }

    const Control<B>& get_control() const
    {
        *_control;
    }
    
    void set_control(std::unique_ptr<Control<B>> control)
    {
        _control = std::move(control);
    }

    void set_delay(double delay)
    {
        _delay = delay;
    }


    void update(const double* control_states, const double time, double* dx_control) override
    {
        double time_delay = time - _delay;

        _navigation->update(time_delay);

        _guidance->update(time_delay,
            _navigation->get_estimated_state());

        _control->update(_navigation.get_estimated_state(),
            _guidance->get_desired_state(),
            control_states,
            time_delay,
            dx_control);
    }
};
