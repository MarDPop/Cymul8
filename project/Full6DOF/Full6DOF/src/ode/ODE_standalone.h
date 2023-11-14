#pragma once

#define _GNU_SOURCE
#define _IOSC11_SOURCE 

#include "../util/Vector.h"
#include <algorithm>
#include <cstdlib>
#include <new>

template<class T, typename Float>
class ode
{

    static_assert(std::is_arithmetic<Float>::value, "Must use a arithmetic type");

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
        Float min_timestep;
        Float max_timestep;
        Float relative_error;
        Float* const inv_absolute_error;

        step_options() = delete;
        step_options(unsigned N) :
            min_timestep(1.0e-6),
            max_timestep(1.0),
            relative_error(1e-6),
            inv_absolute_error(aligned_alloc(32,N*sizeof(Float))) {}

        ~step_options()
        {
            free(inv_absolute_error);
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

        void reserve(unsigned n)
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
            if (*it == _times.end())
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

    const unsigned _N;

    const unsigned _state_bytes;

    T& dynamics;

    step_options options;

private:

    alignas(32) Float* const _state;

    alignas(32) Float* _state_tmp = nullptr;

    alignas(32) Float* _state_rate = nullptr;

    Float _time; // could be faster if time in state, but probably insignificant for anything but smallest of odes

    Float _dt;

public:

    void euler_step()
    {
        this->dynamics(this->_state, this->_time, this->_state_rate);
        for (unsigned i = 0; i < _N; i++)
        {
            this->_state[i] += this->_state_rate[i] * this->_dt;
        }
        this->_time += this->_dt;
    }

    void huen_step()
    {
        const Float dt_half = this->_dt * static_cast<Float>(0.5);

        this->dynamics(this->_state, this->_time, this->_state_rate);
        for (unsigned i = 0; i < _N; i++)
        {
            this->_state[i] += this->_state_rate[i] * this->_dt;
        }
        this->_time += this->_dt;

        this->dynamics(this->_state, this->_time, this->_state_rate + _N);

        for (unsigned i = 0; i < _N; i++)
        {
            this->_state[i] += (this->_state_rate[i + _N] - this->_state_rate[i]) * dt_half;
        }
    }

    void rk4_step()
    {
        const Float dt_half = this->_dt * static_cast<Float>(0.5);
        Float* k2 = this->_state_rate + _N;
        Float* k3 = k2 + _N;
        Float* k4 = k3 + _N;

        memcpy(this->_state_tmp, this->_state, this->_state_bytes);

        this->dynamics(this->_state, this->_time, this->_state_rate);

        for (unsigned i = 0; i < _N; i++)
        {
            this->_state[i] += this->_state_rate[i]*dt_half;
        }
        this->_time += dt_half;

        this->dynamics(this->_state, this->_time, k2);
        for (unsigned i = 0; i < _N; i++)
        {
            this->_state[i] = this->_state_tmp[i] + k2[i]*dt_half;
        }

        this->dynamics(this->_state, this->_time, k3);
        for (unsigned i = 0; i < _N; i++)
        {
            this->_state[i] = this->_state_tmp[i] + k3[i]*this->_dt;
        }

        this->_time += dt_half;
        this->dynamics(this->_state, this->_time, k4);
        const Float dt_6 = this->_dt * static_cast<Float>(0.166666666666666666667);
        for (unsigned i = 0; i < _N; i++)
        {
            this->_state[i] = this->_state_tmp[i] + (this->_state_rate[i] + k4[i] + 2.0*(k2[i] + k3[i]))*dt_6;
        }
    }

    void euler_huen_step()
    {
        const Float dt_half = this->_dt * static_cast<Float>(0.5);
        Float t0 = this->_time;

        memcpy(this->_state_tmp, this->_state, this->_state_bytes);
        this->dynamics(this->_state, this->_time, this->_state_rate);

        constexpr unsigned MAX_ITER = 10u;
        for (unsigned _ = 0; _ < MAX_ITER; _++)
        {
            for (unsigned i = 0; i < _N; i++)
            {
                this->_state[i] = this->_state_tmp[i] + this->_state_rate[i]*this->_dt;
            }
            this->_time = t0 + this->_dt;

            this->dynamics(this->_state, this->_time, this->_state_rate + _N);

            Float max_err = 0.0;
            for (unsigned i = 0; i < _N; i++)
            {
                Float delta = (this->_state_rate[i + _N] - this->_state_rate[i])*dt_half;
                this->_state[i] += delta;
                max_err = std::max(max_err, fabs(delta * this->options.inv_absolute_error[i]));
            }

            constexpr Float dt_factor = 0.9;
            this->_dt *= dt_factor * std::min(std::max(1.0 / max_err, 0.2), 2.0);

            if (max_err < 1.0)
            {
                break;
            }
        }
    }

    void rk23_step()
    {
        Float* k2 = this->_state_rate + _N;
        Float* k3 = k2 + _N;
        Float* k4 = k3 + _N;
        Float* z = this->_state_tmp + _N;

        const Float t0 = this->_time;
        memcpy(this->_state_tmp, this->_state, this->_state_bytes);
        memcpy(this->_state_rate, k4, this->_state_bytes); // FSAL

        constexpr unsigned MAX_ITER = 10u;
        for (unsigned _ = 0; _ < MAX_ITER; _++)
        {
            Float dt_tmp = this->_dt*static_cast<Float>(0.5);
            for (unsigned i = 0; i < _N; i++)
            {
                this->_state[i] = this->_state_tmp[i] + this->_state_rate[i]*dt_tmp;
            }
            this->_time = t0 + dt_tmp;

            this->dynamics(this->_state, this->_time, k2);
            dt_tmp = this->_dt*static_cast<Float>(0.75);
            for (unsigned i = 0; i < _N; i++)
            {
                this->_state[i] = _state_tmp[i] + k2[i]*dt_tmp;
            }
            this->_time = t0 + dt_tmp;

            this->dynamics(this->_state, this->_time, k3);

            this->_time = t0 + this->_dt;

            dt_tmp = this->_dt*static_cast<Float>(0.333333333333333333333);
            for (unsigned i = 0; i < _N; i++)
            {
                this->_state[i] = _state_tmp[i] + ((this->_state_rate[i] + k2[i] + k3[i]) +
                    (k3[i] - this->_state_rate[i])*static_cast<Float>(0.333333333333333333333))*dt_tmp;
            }

            this->dynamics(this->_state, this->_time, k4);

            Float dt_1 = this->_dt*static_cast<Float>(0.29166666666666666666666666666);
            Float dt_2 = this->_dt*static_cast<Float>(0.25);
            Float dt_4 = this->_dt*static_cast<Float>(0.125);
            for (unsigned i = 0; i < _N; i++)
            {
                z[i] = this->_state_tmp[i] + (this->_state_rate[i]*dt_1 + k2[i]*dt_2 + k3[i]*dt_tmp + k4[i]*dt_4);
            }

            dt_1 = 0.0;
            for (unsigned i = 0; i < _N; i++)
            {
                dt_2 = fabs((z[i] - this->_state[i])*this->options.inv_absolute_error[i]);
                dt_1 = std::max(dt_1, dt_2);
            }

            constexpr Float dt_factor = 0.9;
            this->_dt *= dt_factor*std::min(std::max(1.0/dt_2, 0.2), 2.0);

            if (dt_1 < 1.0)
            {
                break;
            }
        }
    }

    ode(T& __dynamics) :
        _N(__dynamics.get_num_states()),
        _state_bytes(_N*sizeof(Float)),
        dynamics(__dynamics),
        options(_N),
        _state(new Float[_N]),
        _state_tmp(new Float[_N]),
        _state_rate(new Float[_N]),
        _time(0.0),
        _dt(1.0)
    {}

    ~ode()
    {
        delete[](_state, std::align_val_t(32));
        delete[](_state_tmp, std::align_val_t(32));
        delete[](_state_rate, std::align_val_t(32));
    }

    const double* get_state()
    {
        return _state;
    }

    const double* get_state_rate()
    {
        return _state_rate;
    }

    /**
    * @brief run and record
    */
    recording run(run_options& options)
    {
        // Initialize recording
        recording record(options);

        // Initial state and time
        memcpy(_state, options.initial_state.data(), _state_bytes);
        _time = options.start_time;
        _dt = options.initial_time_step;

        record.emplace_back(_state, _N);

        delete[] (_state_tmp, std::align_val_t(32));
        delete[] (_state_rate, std::align_val_t(32));

        void (ode::*step)(void) = nullptr;
        switch (options.step)
        {
        case STEP::EULER:
            step = &ode::euler_step;
            _state_rate = new(std::align_val_t(32)) Float[_N];
            break;
        case STEP::HUEN:
            step = &ode::huen_step;
            _state_rate = new(std::align_val_t(32)) Float[_N];
            break;
        case STEP::RK4:
            step = &ode::rk4_step;
            _state_rate = new(std::align_val_t(32)) Float[_N*5];
            break;
        case STEP::EULER_HUEN:
            step = &ode::huen_step;
            _state_rate = new(std::align_val_t(32)) Float[_N*2];
            break;
        case STEP::RK23:
            step = &ode::rk23_step;
            _state_rate = new(std::align_val_t(32)) Float[_N*4];
            _state_tmp = new(std::align_val_t(32)) Float[_N*2];
            // need to initialize K1 for FSAL
            dynamics(_state, _time, _state_rate + 3*_N);
            break;
        }

        Float time_record = 0;
        while (_time < options.end_time && T.valid())
        {
            this->step();
            if (_time > time_record)
            {
                record.emplace_back(_state, _N);
                time_record += options.record_interval;
            }
        }

        record.emplace_back(_state, _N);

        return record;
    }
};
