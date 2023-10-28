#pragma once

#include <exception>
#include <algorithm>
#include <vector>
#include <array>

template<typename Float, unsigned NCOLS>
class Table
{
	std::vector<Float> _x;

	std::vector<std::array<Float, NCOLS>> _v;

	std::array<Float, NCOLS> interpolate(Float x) const
	{
		auto it = std::lower_bound(_x.begin(), _x.end(), x);
		auto idx = std::distance(_x.begin(), it);
		std::array<Float, NCOLS> p1;
		if (interp > 1)
		{
			p1 = _v[idx + interp];
		}
		else
		{
			auto dx = *(it + 1) - *it;
			auto delta = x - *it;
			if (interp == NEAREST)
			{
				p1 = _v[idx + (static_cast<Float>(2.0)*delta > dx)];
			}
			else
			{
				auto delta /= dx;
				p1 = _v[idx];
				const auto& p2 = _v[idx + 1];
				switch (interp)
				{
				case LINEAR:
					for (auto i = 0u; i < NCOLS; i++)
					{
						p1[i] += (p2[i] - p1[i])*delta;
					}
					break;
				case CUBIC:
					const auto& p0 = _v[idx - (idx > 0)];
					const auto& p3 = _v[idx + 1 + (idx < (_x.size() - 2))];
					const auto halfx = 0.5*delta;
					for (auto j = 0u; j < NCOLS; j++)
					{
						p1[j] += halfx*(p2[j] - p0[j] +
							delta*(2.0*p0[j] - 5.0*p1[j] + 4.0*p2[j] - p3[j] + 
								delta*(3.0*(p1[j] - p2[j]) + p3[j] - p0[j])));
					}
					break;
				}
				
			}
			
		}
		return p1;
	}

	std::array<Float, NCOLS> extrapolate(Float x, Float dx) const
	{
		
		unsigned idx = (dx > 0) * (_x.size() - 2);
		dx /= (_x[idx + 1] - _x[idx]);
		std::array<Float, NCOLS> p1 = _v[idx];
		const auto& p2 = _v[idx + 1];
		for (auto i = 0u; i < NCOLS; i++)
		{
			p1[i] += (p2[i] - p1[i])*dx;
		}
		return p1;
	}

public:

	enum INTERPOLATION
	{
		FLOOR = 0,
		CEIL,
		NEAREST,
		LINEAR,
		CUBIC
	};

	enum EXTRAPOLATION
	{
		HOLD_LAST_VALUE = 0,
		LINEAR
	};

	Table() {}

	void set(unsigned idx, const Float& x, const std::array<Float, NCOLS>& v)
	{
		if (idx > _x.size())
		{
			_x.resize(idx);
			_v.resize(idx);
		}
		_x[idx] = x;
		_v[idx] = v;
	}

	void add(const Float& x, const std::array<Float, NCOLS>& v)
	{
		_x.push_back(x);
		_v.push_back(v);
	}

	std::array<Float, NCOLS> get(Float x, 
									INTERPOLATION interp = LINEAR,
									EXTRAPOLATION extrap = HOLD_LAST_VALUE) const
	{
		if (x < _x[0])
		{
			if (extrap == HOLD_LAST_VALUE)
			{
				return _v[0];
			}
			return this->extrapolate(x, x - _x[0]);
		}
		if (x > _x.back())
		{
			if (extrap == HOLD_LAST_VALUE)
			{
				return _v.back();
			}
			return this->extrapolate(x, x - _x.back());
		}
		return this->interpolate(x);
	}
};

template<typename Float, unsigned NROWS, unsigned NCOLS>
class LinearTable
{
	Float _x[NROWS];

	Float _v[NROWS][NCOLS];

	Float _dv[NROWS][NCOLS];

	constexpr std::size_t ROW_BYTES = NCOLS * sizeof(Float);

	void finish()
	{
		auto* v = _v;
		auto* dv = _dv;
		for (auto i = 1u; i < NROWS; i++)
		{
			auto dx = 1.0 / (_x[i] - _x[i - 1]);
			for (auto j = 0u; j < NCOLS; j++)
			{
				dv[j] = (v[j + NCOLS] - v[j]) * dx;
			}
		}
	}

public:

	void set(const Float* x, const Float* v)
	{
		memcpy(_x, x, NROWS*sizeof(Float));
		memcpy(_v, v, NROWS*ROW_BYTES);
		finish();
	}

	void set(const std::vector<Float>& x, const std::vector<std::array<Float,NCOLS>>& v)
	{
		assert(x.size() >= NROWS && v.size() >= NROWS);
		memcpy(_x, x.data(), NROWS * sizeof(Float));
		memcpy(_v, &v[0][0], NROWS * ROW_BYTES);
		finish();
	}

	void get(Float x, Float* v) const
	{
		auto it = std::lower_bound(_x, _x + _length, x);
		auto delta = x - *it;
		auto idx = it - _x;
		for (auto i = 0u; i < NCOLS; i++)
		{
			v[i] = _v[i] + delta*_dv[i];
		}
	}

};

template<typename Float>
class DynamicTable
{
	Float* const _x;

	Float* const _v;

	Float* _p;

	void (DynamicTable::*_interpolate)(Float, Float*);

	void (DynamicTable::*_finish)(void);

	constexpr unsigned NCOEF_CUBIC = 3;

public:

	const unsigned NROWS;

	const unsigned NCOLS;

	const unsigned ROW_BYTES;

	DynamicTable(	unsigned __NROWS,
					unsigned __NCOLS,
					Table::INTERPOLATION method) : 
		NROWS(__NROWS),
		NCOLS(__NCOLS),
		ROW_BYTES(__NCOLS*sizeof(Float)),
		_x(new Float[__NROWS]),
		_v(new Float[__NROWS*__NCOLS])
	{
		if (method == CUBIC && __NROWS < 4)
		{
			method = _NROWS > 1 ? LINEAR : FLOOR;
		}
		if (method == LINEAR && __NROWS < 2)
		{
			method = FLOOR;
		}

		switch (method)
		{
		case NEAREST:
			_p = new Float[__NROWS];
			_finish = &Table::finish_nearest();
			_interpolate = &Table::interpolate_nearest();
			break;
		case LINEAR:
			_p = new Float[__NROWS * __NCOLS];
			_finish = &Table::finish_linear();
			_interpolate = &Table::interpolate_nearest();
			break;
		case CUBIC:
			_p = new Float[__NROWS * __NCOLS * NCOEF_CUBIC];
			_finish = &Table::finish_cubic();
			_interpolate = &Table::interpolate_nearest();
			break;
		default:
			_p = nullptr;
			break;
		}
	}

	~DynamicTable()
	{
		delete[] _x;
		delete[] _v;
		delete[] _p;
	}

	void set(const Float* x, const Float* v)
	{
		memcpy(_x, x, NROWS * sizeof(Float));
		memcpy(_v, v, NROWS * ROW_BYTES);
		finish(_p, _x, _v);
	}

	void set(const std::vector<Float>& x, const std::vector<std::array<Float, NCOLS>>& v)
	{
		assert(x.size() >= NROWS && v.size() >= NROWS);
		memcpy(_x, x.data(), NROWS * sizeof(Float));
		memcpy(_v, &v[0][0], NROWS * ROW_BYTES);
		finish(_p, _x, _v);
	}

	void get(Float x, Float* v) const
	{
		if (x < _x[0])
		{

		}
		if (x > _x.back())
		{

		}
		_interpolate(x, v);
	}

private:

	void finish_nearest()
	{
		for (auto i = 1u; i < NROWS; i++)
		{
			_p[i - 1] = 0.5*(_x[i] + _x[i - 1]);
		}
	}

	void finish_linear()
	{
		auto* v = _v;
		auto* dv = _p;
		for (auto i = 1u; i < NROWS; i++)
		{
			auto dx = 1.0 / (_x[i] - _x[i - 1]);
			for (auto j = 0u; j < NCOLS; j++)
			{
				dv[j] = (v[j + NCOLS] - v[j]) * dx;
			}
		}
	}

	void finish_cubic()
	{
		if (NROWS < 3)
		{
			finish_linear();
		}
		auto* v = _v;
		auto* p = _p;

		const auto NEXT_ROW = 2 * NCOLS;
		const auto LAST_ROW = 3 * NCOLS;

		Float A[9] = {0,0,1,0,1,0,0,1,0};

		Float A_inv[9];
		Float y[3];
		
		auto dx0 = 1.0;
		auto dx1 = 1.0 / (_x[1] - _x[0]);
		auto dx2 = 1.0 / (_x[2] - _x[1]);

		p += ncoefs;

		for (auto i = 1u; i < NROWS - 1; i++)
		{
			auto dx = _x[1] - _x[0];
			auto dx_sq = dx*dx;
			auto twodx = 2*dx;
			auto dx0 = dx1;
			auto dx1 = dx2;
			auto dx2 = 1.0 / (_x[i + 2] - _x[i + 1]);

			A[0] = dx_sq;
			A[1] = dx;
			A[3] = twodx;
			for (auto j = 0u; j < NCOLS; j++)
			{
				dvdx0 = v[j + NCOLS] - v[j];
				dvdx1 = v[j + NEXT_ROW] - v[j + NCOLS];
				dvdx2 = v[j + LAST_ROW] - v[j + NEXT_ROW];

				y[0] = v[j + NCOLS];

				p += NCOEF_CUBIC;
			}
			v += NCOLS;
		}
	}

	void interpolate_nearest(Float x, Float* v)
	{
		auto it = std::lower_bound(_x, _x + NROWS, x);
		auto idx = (it - _x);
		idx += (x > _p[idx]);
		memcpy(v, _v + idx*NCOLS, ROW_BYTES);
	}

	void interpolate_linear(Float x, Float* v)
	{
		auto it = std::lower_bound(_x, _x + NROWS, x);
		auto idx = static_cast<unsigned>(it - _x)*NCOLS;
		const auto* __v = _v + idx;
		const auto* __p = _p + idx;
		const auto delta = x - *it;
		for (auto j = 0u; j < NCOLS; j++)
		{
			v[j] = __v[j] + delta * __p[j];
		}

	}

	void interpolate_cubic(Float x, Float* v)
	{
		auto it = std::lower_bound(_x, _x + NROWS, x);
		auto idx = static_cast<unsigned>(it - _x) * NCOLS;
		const auto* __v = _v + idx;
		const auto* __p = _p + idx;
		const auto delta = x - *it;
		for (auto j = 0u; j < NCOLS; j++)
		{
			v[j] = __v[j] + delta * __p[j];
		}
	}
};