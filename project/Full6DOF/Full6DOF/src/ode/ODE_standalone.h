#pragma once

#include "../util/Vector.h"
#include <algorithm>

template<class T, typename Float>
class ode
{

    static_assert(std::is_arithmetic<T>::value, "Must use a arithmetic type");

public:

    enum class STEP
    {
        EULER,
        HUEN,
        RK4,
        RK23
    };

    struct step_options
    {
        friend class ode;

        Float min_timestep;
        Float max_timestep;
        Float relative_error;
        Float* const inv_absolute_error;

    public:
        
        step_options(unsigned N) :
            inv_absolute_error(new Float[N]) {}

        ~step_options()
        {
            delete[] absolute_error;
        }
    };

    struct run_options
    {
        Float start_time;
        Float end_time;
        Float record_interval;
        Float initial_time_step;
        fixed_vector<Float> initial_state;
        STEP step;

        run_options(unsigned N) :
            initial_state(N)
        {}
    };

    class recording
    {
        friend class ode;

        std::vector<Float> _times;

        std::vector<fixed_vector<Float>> _states;

        const unsigned _N;

        void put(Float time, Float* states)
        {
            _times.push_back(time);
            _states.emplace_back(states, _N);
        }

    public:

        void reserve(unsigned N)
        {
            _times.reserve(n);
            _states.reserve(n);
        }

        recording(unsigned n) : _N(n) {}

        recording(const run_options& options) :
            _N(options.initial_state.size())
        {
            unsigned nRecordings = static_cast<unsigned>(
                (options.end_time - options.start_time) / options.record_interval
                ) + 2u;
            this->reserve(nRecordings);
        }

        const std::vector<Float>& get_times()
        {
            return this->_times;
        }

        const std::vector<fixed_vector<Float>>& get_states()
        {
            return this->_states;
        }

        fixed_vector<Float> get(Float time)
        {
            auto it = std::lower_bound(_times.begin(), _times.end());
            if (*it == times.end())
            {
                return fixed_vector<Float>(_states.end(), _N);
            }
            auto idx = std::distance(it, _times.begin());
            Float delta = (time - *it) / (*(it + 1) - *it);

            fixed_vector<Float> lo(_states[idx], _N);
            auto& hi = _states[idx + 1];
            for (unsigned i = 0; i < _N; i++)
            {
                lo[i] += (hi[i] - lo[i]) * delta;
            }
            return lo;
        }
    };

private:

    Float* const _state; 

    Float* const _state_tmp;

    Float* _state_rate = nullptr;

    Float _time; // could be faster if time in state, but probably insignificant for anything but smallest of odes

    Float _dt;

    const unsigned _N;

    const unsigned _state_bytes;

    static void euler_step(ode& __ode)
    {
        __ode.dynamics(__ode._state, _ode._time, __ode._state_rate);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += __ode._state_rate[i] * __ode._dt;
        }
        __ode._time += __ode._dt;
    }

    static void huen_step(ode& __ode)
    {
        const Float dt_half = __ode._dt * static_cast<Float>(0.5);

        __ode.dynamics(__ode._state, _ode._time, __ode._state_rate);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += __ode._state_rate[i] * __ode._dt;
        }
        __ode._time += __ode._dt;

        __ode.dynamics(__ode._state, _ode._time, __ode._state_rate + _N);
        
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += (__ode._state_rate[i + _N] - __ode._state_rate[i]) * dt_half;
        }
    }

    static void rk4_step(ode& __ode)
    {
        const Float dt_half = __ode._dt * static_cast<Float>(0.5);
        Float* x0 = __ode._state + _N;
        Float* k2 = __ode._state_rate + _N;
        Float* k3 = k2 + _N;
        Float* k4 = k3 + _N;

        memcpy(x0, __ode._state, __ode._state_bytes);

        __ode.dynamics(__ode._state, _ode._time, __ode._state_rate);
        
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += __ode._state_rate[i] * dt_half;
        }
        __ode._time += dt_half;

        __ode.dynamics(__ode._state, _ode._time, k2);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] = x0[i] + k2[i] * dt_half;
        }

        __ode.dynamics(__ode._state, _ode._time, k3);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] = x0[i] + k3[i] * __ode._dt;
        }

        __ode._time += dt_half;
        __ode.dynamics(__ode._state, _ode._time, k4);
        const Float dt_6 = __ode._dt * static_cast<Float>(0.166666666666666666667);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] = x0[i] + (__ode._state_rate[i] + k4[i] + 2.0 * (k2[i] + k3[i])) * dt_6;
        }

    }

    static void euler_huen_step(ode& __ode)
    {
        const Float dt_half = __ode._dt * static_cast<Float>(0.5);
        Float* x0 = __ode._state + _N;
        Float t0 = __ode._time;

        __ode.dynamics(__ode._state, _ode._time, __ode._state_rate);

        constexpr unsigned MAX_ITER = 10u;
        for (unsigned _ = 0; _ < MAX_ITER; _++)
        {
            for (unsigned i = 0; i < _N; i++)
            {
                __ode._state[i] = x0[i] + __ode._state_rate[i] * __ode._dt;
            }
            __ode._time = t0 + __ode._dt;

            __ode.dynamics(__ode._state, _ode._time, __ode._state_rate + _N);

            Float max_err = 0.0;
            for (unsigned i = 0; i < _N; i++)
            {
                Float delta = (__ode._state_rate[i + _N] - __ode._state_rate[i]) * dt_half;
                __ode._state[i] += delta;
                max_err = std::max(max_err, fabs(delta * __ode.options.inv_absolute_error[i]));
            }

            constexpr Float dt_factor = 0.9;
            __ode._dt *= dt_factor * std::min(std::max(1.0 / max_err, 0.2), 2.0);

            if (max_err < 1.0)
            {
                break;
            }
        }
    }

    static void rk23_step(ode& __ode)
    {
        Float* x0 = __ode._state + _N;
        Float* k2 = __ode._state_rate + _N;
        Float* k3 = k2 + _N;
        Float* k4 = k3 + _N;
        Float* z = x0 + _N;

        const Float t0 = __ode._time;
        memcpy(x0, __ode._state, __ode._state_bytes);
        memcpy(__ode._state_rate, k4, __ode._state_bytes); // FSAL

        constexpr unsigned MAX_ITER = 10u;
        for (unsigned _ = 0; _ < MAX_ITER; _++)
        {
            Float dt = __ode._dt * static_cast<Float>(0.5);
            for (unsigned i = 0; i < _N; i++)
            {
                __ode._state[i] = x0[i] + __ode._state_rate[i] * dt;
            }
            _ode._time = t0 + dt;

            __ode.dynamics(__ode._state, _ode._time, k2);
            dt = __ode._dt * static_cast<Float>(0.75);
            for (unsigned i = 0; i < _N; i++)
            {
                __ode._state[i] = x0[i] + k2[i] * dt;
            }
            _ode._time = t0 + dt;

            __ode.dynamics(__ode._state, _ode._time, k3);

            _ode._time = t0 + __ode._dt;

            dt = __ode._dt * static_cast<Float>(0.333333333333333333333);
            for (unsigned i = 0; i < _N; i++)
            {
                __ode._state[i] = x0[i] + ((__ode._state_rate[i] + k2[i] + k3[i]) + 
                    (k3[i] - __ode._state_rate[i]) * static_cast<Float>(0.333333333333333333333)) * dt;
            }

            __ode.dynamics(__ode._state, _ode._time, k4);

            Float dt_1 = __ode._dt * static_cast<Float>(0.29166666666666666666666666666);
            Float dt_2 = __ode._dt * static_cast<Float>(0.25);
            Float dt_4 = __ode._dt * static_cast<Float>(0.125);
            for (unsigned i = 0; i < _N; i++)
            {
                z[i] = x0[i] + (__ode._state_rate[i] * dt_1 + k2[i] * dt_2 + k3[i] * dt + k4[i] * dt_4);
            }

            dt_1 = 0.0;
            for (unsigned i = 0; i < _N; i++)
            {
                dt_2 = fabs( (z[i] - __ode._state[i])*__ode.options.inv_absolute_error[i] );
                dt_1 = std::max(dt_1, dt_2);
            }

            constexpr Float dt_factor = 0.9;
            __ode._dt *= dt_factor * std::min(std::max( 1.0/dt_2, 0.2), 2.0);

            if (dt_1 < 1.0)
            {
                break;
            }
        }
    }

public:

    step_options options;

    T dynamics;

    ode() : 
        _state(new Float[_N]),
        _state_tmp(new Float[_N]),
        _N(T.get_num_states()),
        _state_bytes(_N*sizeof(Float)),
        options(_N)
    {}

    ~ode()
    {
        delete[] _state;
        delete[] _state_rate;
    }

    inline const double* get_state()
    {
        return _state;
    }

    inline const double* get_state_rate()
    {
        return _state_rate;
    }

    /**
    * @brief run and record
    */
    inline recording run(run_options& options)
    {
        // Initialize recording
        recording record(options);

        // Initial state and time
        memcpy(_state, options.initial_state.data(), _state_bytes);
        _time = options.start_time;
        _dt = options.initial_time_step;

        record.emplace_back(_state, _N);

        void (*step)(ode&) = &ode::euler_step;
        switch (options.step)
        {
        case HUEN:
            step = &ode::huen_step;
            this->_state_rate = new Float[_N];
            break;
        case RK4:
            step = &ode::rk4_step;
            this->_state_rate = new Float[_N*5];
            break;
        case EULER_HUEN:
            step = &ode::huen_step;
            this->_state_rate = new Float[_N*2];
            break;
        case RK23:
            step = &ode::rk23_step;
            this->_state_rate = new Float[_N*4];
            // need to initialize K1 for FSAL
            dynamics(_state, _time, _state_rate + 4*N);
            break;
        }

        while (_time < options.end_time && T.valid())
        {
            step(*this);
            if (_time > time_record)
            {
                record.emplace_back(_state, _N);
            }
        }

        record.emplace_back(_state, _N);

        delete[] this->_k;

        return record;
    }
};
