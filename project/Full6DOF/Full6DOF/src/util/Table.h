#pragma once

#include <exception>

class Table
{

public:

};

template<typename Float>
class LinearTable
{
	Float* _x;

	Float* _v;

	Float* _dv;

	unsigned _length;

public:

	LinearTable(unsigned n) : _length(n)
	{
		_x = new Float[n];
		_v = new Float[n];
		_dv = new Float[n];
	}

	LinearTable(const Float* const x, const Float* const v, unsigned n)
	{
		_x = new Float[n];
		_v = new Float[n];
		_dv = new Float[n];
		auto bytes = n*sizeof(Float);
		memcpy(_x, x, bytes);
		memcpy(_v, v, bytes);
		for (auto i = 1u; i < n; i++)
		{
			_dv[i - 1] = (_v[i] - _v[i - 1]) / (_x[i] - _x[i - 1]);
		}
	}

	~LinearTable()
	{
		delete[] _x;
		delete[] _v;
		delete[] _dx;
	}

	/**
	* Don't use this 
	*/
	void add(const Float x, const Float v)
	{
		auto i = _length++;
		auto bytes = _length * sizeof(Float);
		Float* _tmp_x = _x;
		Float* _tmp_v = _v;
		Float* _tmp_dv = _dv;
		_x = realloc(_x, bytes);
		_v = realloc(_v, bytes);
		_dv = realloc(_dv, bytes);
		if (!dv)
		{
			_x = _tmp_x;
			_v = _tmp_v;
			_dv = _tmp_dv;
			_length--;
			throw std::runtime_error("Couldn't reserve memory");
		}

		_x[i] = x;
		_v[i] = v;
		_dv[i - 1] = (_v[i] - _v[i - 1]) / (_x[i] - _x[i - 1]);
	}

	void set(unsigned i, const Float x, const Float v)
	{
		assert(i < _length);
		_x[i] = x;
		_v[i] = v;
		if (i != _length - 1)
		{
			_dv[i] = (_v[i + 1] - _v[i]) / (_x[i + 1] - _x[i]);
		}
	}

	Float get(Float x)
	{
		auto it = std::lower_bound(_x, _x + _length, x);
		auto delta = x - *it;
		auto idx = it - _x;
		return _v[idx] + delta*_dv[idx];
	}

};

template<typename Float>
class CubicTable
{

	Float* _x;

	Float* _v;

	Float* _dv;

public:

};

template<typename Float, unsigned N>
class LinearTable
{
	Float* _x;

	Float* _v;

	Float* _dv;

	unsigned _length;

	constexpr ROW_BYTES = N * sizeof(Float);

public:

	LinearTable(unsigned n) : _length(n)
	{
		_x = new Float[n];
		_v = new Float[n*N];
		_dv = new Float[n*N];
	}

	LinearTable(const Float* const x, const Float* const v, unsigned n)
	{
		_x = new Float[n];
		_v = new Float[n*N];
		_dv = new Float[n*N];

		memcpy(_x, x, n*sizeof(Float));
		memcpy(_v, v, n*ROW_BYTES);
		
		auto dv_ = _dv;
		auto v_ = _v;
		for (auto i = 1u; i < n; i++)
		{
			auto den = static_cast<Float>(1.0) / (_x[i] - _x[i - 1]);
			for (auto j = 0u; j < N; j++)
			{
				dv_[j] = (v_[N + j] - v_[j])*den;
				dv_ += N;
				v_ += N;
			}
			
		}
	}

	~LinearTable()
	{
		delete[] _x;
		delete[] _v;
		delete[] _dx;
	}

	/**
	* Don't use this
	*/
	void add(Float x, Float* v)
	{
		auto i = _length++;
		auto bytes = _length * sizeof(Float);
		Float* _tmp_x = _x;
		Float* _tmp_v = _v;
		Float* _tmp_dv = _dv;
		_x = realloc(_x, _length*sizeof(Float));
		_v = realloc(_v, _length*N*sizeof(Float));
		_dv = realloc(_dv, _length*N*sizeof(Float));
		if (!dv)
		{
			_x = _tmp_x;
			_v = _tmp_v;
			_dv = _tmp_dv;
			_length--;
			throw std::runtime_error("Couldn't reserve memory");
		}

		_x[i] = x;
		memcpy(_v + i * N, v, ROW_BYTES);
		_dv[i - 1] = (_v[i] - _v[i - 1]) / (_x[i] - _x[i - 1]);
	}

	void set(unsigned i, const Float x, const Float* const v)
	{
		assert(i < _length);
		_x[i] = x;
		auto v_ = = _v + i * N;
		memcpy(v_, v, ROW_BYTES);
		if (i != _length - 1)
		{
			auto dv_ = _dv + i*N;
			
			auto den = static_cast<Float>(1.0) / (_x[i + 1] - _x[i]);
			for (auto j = 0u; j < N; j++)
			{
				dv_[j] = (v_[N + j] - v_[j]) * den;
			}
		}
	}

	void get(Float x, Float* v)
	{
		auto it = std::lower_bound(_x, _x + _length, x);
		auto delta = x - *it;
		auto idx = (it - _x)*N;
		auto v_ = = _v + idx;
		auto dv_ = = _dv + idx;
		for (auto j = 0u; j < N; j++)
		{
			v[j] = v_[j] + delta * dv_[j];
		} 
	}

};