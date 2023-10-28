#pragma once

#include <exception>
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include "fast_math.h"

template<typename Float>
class Table
{
public:
	const unsigned NCOLS;

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

private:

	std::vector<std::string> _column_names;

	std::vector<Float> _x;

	std::vector<std::vector<Float>> _v;

	std::vector<Float> interpolate(Float x) const
	{
		auto it = std::lower_bound(_x.begin(), _x.end(), x);
		auto idx = std::distance(_x.begin(), it);
		std::vector<Float> p1(NCOLS);
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

	std::vector<Float> extrapolate(Float x, Float dx) const
	{
		unsigned idx = (dx > 0) * (_x.size() - 2);
		dx /= (_x[idx + 1] - _x[idx]);
		std::vector<Float> p1 = _v[idx];
		const auto& p2 = _v[idx + 1];
		for (auto i = 0u; i < NCOLS; i++)
		{
			p1[i] += (p2[i] - p1[i])*dx;
		}
		return p1;
	}

public:

	Table(unsigned __NCOLS) : 
			NCOLS(__NCOLS),
			_column_names(NCOLS)
	{}

	const std::vector<std::string>& get_column_names()
	{
		return _column_names;
	}

	void set_column_name(unsigned idx, std::string name)
	{
		_column_names[idx] = name;
	}

	void set(unsigned idx, const Float& x, const std::vector<Float>& v)
	{
		assert(v.size() == NCOLS);

		if (idx > _x.size())
		{
			_x.resize(idx);
			_v.resize(idx);
		}
		_x[idx] = x;
		_v[idx] = v;
	}

	void add(const Float& x, const std::vector<Float>& v)
	{
		assert(v.size() == NCOLS);

		_x.push_back(x);
		_v.push_back(v);
	}

	std::vector<Float> get(Float x, 
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
public:

	constexpr unsigned NCOEF_CUBIC = 3;

	const unsigned NROWS;

	const unsigned NCOLS;

	const unsigned ROW_BYTES;

	const unsigned LAST_IDX;

private: 
	void (DynamicTable::* _extrapolate)(Float, Float*);

	void (DynamicTable::* _interpolate)(Float, Float*);

	void (DynamicTable::* _finish)(void);

	Float* const _x;

	Float* const _v;

	Float* _p;

	void finish_nearest()
	{
		for (auto i = 1u; i < NROWS; i++)
		{
			_p[i - 1] = 0.5 * (_x[i] + _x[i - 1]);
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
				*dv++ = (*(v + NCOLS) - *v) * dx;
				v++;
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

		Float A[4];
		Float A_inv[4];
		Float y[2];

		// First Entry
		auto half_dx0 = 1.0;
		auto half_dx1 = 0.5 / (_x[1] - _x[0]);
		auto half_dx2 = 0.5 / (_x[2] - _x[1]);

		auto dx = _x[1] - _x[0];
		auto dx_sq = dx * dx;

		A[0] = dx_sq * dx;
		A[1] = dx_sq;
		A[2] = 3.0 * dx_sq;
		A[3] = 2.0 * dx;
		inverse2x2(A, A_inv);
		for (auto j = 0u; j < NCOLS; j++)
		{
			auto dvdx1 = (v[j + NCOLS] - v[j]) * dx1;
			auto dvdx2 = (v[j + NEXT_ROW] - v[j + NCOLS]) * dx2;

			p[2] = dvdx1 * 2.0;

			y[0] = (_v[j + NCOLS] - p[2] * dx - _v[j];
			y[1] = (dvdx2 + dvdx1) - p[2];

			mult2x2(A_inv, y, p);

			p += NCOEF_CUBIC;
		}
		v += NCOLS;
		// Middle entries
		for (auto i = 1u; i < NROWS - 1; i++)
		{
			dx = _x[i + 1] - _x[i];
			dx_sq = dx * dx;
			half_dx0 = half_dx1;
			half_dx1 = half_dx2;
			half_dx2 = 0.5 / (_x[i + 2] - _x[i + 1]);

			A[0] = dx_sq * dx;
			A[1] = dx_sq;
			A[2] = 3.0 * dx_sq;
			A[3] = 2.0 * dx;
			inverse2x2(A, A_inv);
			for (auto j = 0u; j < NCOLS; j++)
			{
				auto dvdx0 = (v[j] - v[j - NCOLS]) * dx0;
				auto dvdx1 = (v[j + NCOLS] - v[j]) * dx1;
				auto dvdx2 = (v[j + NEXT_ROW] - v[j + NCOLS]) * dx2;

				p[2] = dvdx1 + dvdx0; // c or dvdx at point x1

				y[0] = (_v[j + NCOLS] - p[2] * dx - _v[j];
				y[1] = (dvdx2 + dvdx1) - p[2];

				mult2x2(A_inv, y, p);

				p += NCOEF_CUBIC;
			}
			v += NCOLS;
		}

		// Last entry
		dx = _x.back() - _x[_size() - 2];
		dx_sq = dx * dx;

		A[0] = dx_sq * dx;
		A[1] = dx_sq;
		A[2] = 3.0 * dx_sq;
		A[3] = 2.0 * dx;
		inverse2x2(A, A_inv);
		for (auto j = 0u; j < NCOLS; j++)
		{
			auto dvdx0 = (v[j] - v[j - NCOLS]) * dx0;
			auto dvdx1 = (v[j + NCOLS] - v[j]) * dx1;

			p[2] = dvdx1 + dvdx0;

			y[0] = (_v[j + NCOLS] - p[2] * dx - _v[j];
			y[1] = dvdx1 * 2.0 - p[2];

			mult2x2(A_inv, y, p);

			p += NCOEF_CUBIC;
		}
	}

	void interpolate_nearest(Float x, Float* v)
	{
		auto it = std::lower_bound(_x, _x + NROWS, x);
		auto idx = (it - _x);
		idx += (x > _p[idx]);
		memcpy(v, _v + idx * NCOLS, ROW_BYTES);
	}

	void interpolate_linear(Float x, Float* v)
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

	void interpolate_cubic(Float x, Float* v)
	{
		auto it = std::lower_bound(_x, _x + NROWS, x);
		auto idx = static_cast<unsigned>(it - _x) * NCOLS;
		const auto* __v = _v + idx;
		const auto* __p = _p + idx;
		const auto dx = x - *it;
		// could be faster with simd and precomputing dx^2 dx^3
		for (auto j = 0u; j < NCOLS; j++)
		{
			v[j] = __v[j] + dx * (__p[2] + dx * (__p[1] + dx * __p[0]);
			__p += NCOEF_CUBIC;
		}
	}

	void extrapolate_hold_last_value(Float dx, Float* v)
	{
		memcpy(v, _v + LAST_IDX * (dx > 0), ROW_BYTES);
	}

	void extrapolate_linear(Float dx, Float* v)
	{
		Float* row;
		Float dx_ratio;
		if (dx > 0)
		{
			row = _v + (LAST_IDX - NCOLS);
			dx_ratio = dx / (_x[NROWS - 1] - _x[NROWS - 2]);
		}
		else
		{
			row = _v;
			dx_ratio = dx / (_x[1] - _x[0]);
		}
		for (auto i = 0u; i < NCOLS; i++)
		{
			v[i] = row[i] + (row[i + NCOLS] - row[i])*dx_ratio;
		}
	}

public:

	DynamicTable(unsigned __NROWS,
				unsigned __NCOLS,
				Table::INTERPOLATION interp = LINEAR,
				Table::EXTRAPOLATION extrap = HOLD_LAST_VALUE) :
		_x(new Float[__NROWS]),
		_v(new Float[__NROWS*__NCOLS]),
		LAST_IDX((__NROWS-1)*NCOLS),
		NROWS(__NROWS),
		NCOLS(__NCOLS),
		ROW_BYTES(__NCOLS * sizeof(Float))
	{
		if (interp == CUBIC && __NROWS < 4)
		{
			interp = _NROWS > 1 ? LINEAR : FLOOR;
		}
		if (interp == LINEAR && __NROWS < 2)
		{
			interp = FLOOR;
		}

		switch (interp)
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

		if (extrap == LINEAR && __NROWS < 2)
		{
			extrap = HOLD_LAST_VALUE;
		}

		if (extrap == LINEAR)
		{
			_extrapolate = &Table::extrapolate_linear();
		}
		else
		{
			_extrapolate = &Table::extrapolate_hold_last_value();
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
			_extrapolate(x - _x[0], v);
			return;
		}
		if (x > _x.back())
		{
			_extrapolate(x - _x.back(), v);
			return;
		}
		_interpolate(x, v);
	}

	void interpolate(Float x, Float* v) const
	{
		_interpolate(x, v);
	}

};