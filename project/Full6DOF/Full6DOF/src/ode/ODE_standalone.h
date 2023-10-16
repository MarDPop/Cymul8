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
        Float min_timestep;
        Float max_timestep;
        Float relative_error;
        Float* const absolute_error;

        step_options(unsigned N) :
            absolute_error(new Float[N]) {}
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
        Float* initial_state;
        STEP step;

        run_options(unsigned N) :
            initial_state(new Float[N]) {}

        ~run_options()
        {
            delete[] initial_state;
        }
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

        recording(unsigned n) : _N(n) {}

        void reserve(unsigned N)
        {
            _times.reserve(n);
            _states.reserve(n);
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

    Float* const _state_rate;

    Float* _k = nullptr;

    Float _time;

    Float _dt;

    const unsigned _N;

    const unsigned _state_bytes;

    static void euler_step(ode& __ode)
    {
        __ode.dynamics(__ode._state, __ode._state_rate);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += __ode._state_rate[i] * __ode._dt;
        }
        __ode._time += __ode._dt;
    }

    static void huen_step(ode& __ode)
    {
        const Float dt_half = __ode._state_rate[_N] * static_cast<Float>(0.5);

        __ode.dynamics(__ode._state, __ode._state_rate);
        memcpy(__ode._k, __ode._state_rate, __ode._state_bytes);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += __ode._state_rate[i] * __ode._dt;
        }
        __ode._time += __ode._dt;

        __ode.dynamics(__ode._state, __ode._state_rate);
        
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += (__ode.get_state_rate[i] - __ode._k[i]) * dt_half;
        }
    }

    static void rk4_step(ode& __ode)
    {
        const Float dt_half = __ode._state_rate[_N] * static_cast<Float>(0.5);

        memcpy(__ode._k, __ode._state, __ode._state_bytes);

        __ode.dynamics(__ode._state, __ode._state_rate);
        
        Float* k1 = __ode._k + _N;
        memcpy(k1, __ode._state_rate, __ode._state_bytes);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] += __ode._state_rate[i] * dt_half;
        }
        __ode._time += dt_half;

        __ode.dynamics(__ode._state, __ode._state_rate);
        Float* k2 = k1 + _N;
        memcpy(k2, __ode._state_rate, __ode._state_bytes);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] = dx0[i] + k2[i] * dt_half;
        }

        __ode.dynamics(__ode._state, __ode._state_rate);
        Float* k3 = k2 + _N;
        memcpy(k3, __ode._state_rate, __ode._state_bytes);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] = dx0[i] + k3[i] * __ode._dt;
        }

        __ode._state[_N] += dt_half;
        __ode.dynamics(__ode._state, __ode._state_rate);
        Float* k4 = k3 + _N;
        memcpy(k4, __ode._state_rate, __ode._state_bytes);

        const Float dt_6 = __ode._dt * static_cast<Float>(0.166666666666666666667);
        for (unsigned i = 0; i < _N; i++)
        {
            __ode._state[i] = dx0[i] + (k1[i] + k4[i] + 2.0 * (k2[i] + k3[i])) * dt_6;
        }

    }

    static void rk23_step(ode& __ode)
    {
        __ode.dynamics(__ode._state, __ode._state_rate);
        const Float t0 = __ode._state[_N];
        array_wrapper<Float> x0(__ode._state, _N);
        array_wrapper<Float> k1(__ode._state_rate, _N);
        array_wrapper<Float> k2(_N);
        array_wrapper<Float> k3(_N);
        array_wrapper<Float> k4(_N);

        const auto bytes = _N * sizeof(Float);

        constexpr unsigned MAX_ITER = 20u;
        for (unsigned _ = 0; _ < MAX_ITER; _++)
        {
            const Float dt_1 = __ode._state_rate[_N] * static_cast<Float>(0.5);
            for (unsigned i = 0; i < _N; i++)
            {
                __ode._state[i] = dx0[i] + k1[i] * dt_1;
            }
            _ode._state[_N] = t0 + dt1;

            __ode.dynamics(__ode._state, __ode._state_rate);
            memcpy(k2.data(), ode._state_rate, bytes);

            const Float dt_2 = __ode._state_rate[_N] * static_cast<Float>(0.75);
            for (unsigned i = 0; i < _N; i++)
            {
                __ode._state[i] = dx0[i] + k2[i] * dt_2;
            }
            _ode._state[_N] = t0 + dt_2;

            __ode.dynamics(__ode._state, __ode._state_rate);
            memcpy(k3.data(), ode._state_rate, bytes);

            for (unsigned i = 0; i < _N; i++)
            {
                __ode._state[i] = dx0[i] + k2[i] * __ode._state_rate[_N];
            }
            _ode._state[_N] = t0 + __ode._state_rate[_N];
        }
    }

public:

    step_options options;

    T dynamics;

    ode() : _N(T.get_num_states()),
        _state_bytes(_N*sizeof(Float)),
        _state(new Float[_N]),
        _state_rate(new Float[_N]),
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

    inline recording run(run_options& options)
    {

        void (*step)(ode&) = &ode::euler_step;
        switch (options.step)
        {
        case HUEN:
            step = &ode::huen_step;
            this->_k = new Float[_N];
            break;
        case RK4:
            step = &ode::rk4_step;
            this->_k = new Float[_N*4];
            break;
        case RK23:
            step = &ode::rk23_step;
            this->_k = new Float[_N*4];
            break;
        }

        recording record;
        unsigned nRecordings = static_cast<unsigned>(
            (options.end_time - options.start_time) / options.record_interval
            ) + 2u;
        record.times.reserve(nRecordings);
        record.states.reserve(nRecordings);

        memcpy(_state, options.initial_state, _state_bytes);
        _time = options.start_time;
        _dt = options.initial_time_step;

        record.emplace_back(_state, _N);

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
