#pragma once

#include <exception>
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include "fast_math.h"

namespace Tablulate
{
	enum class INTERPOLATION
	{
		FLOOR = 0,
		CEIL,
		NEAREST,
		LINEAR,
		CUBIC
	};

	enum class EXTRAPOLATION
	{
		HOLD_LAST_VALUE = 0,
		LINEAR
	};
}

class Table
{
public:
	const unsigned NCOLS;

private:

	std::vector<std::string> _column_names;

	std::vector<double> _x;

	std::vector<std::vector<double>> _v;

	std::vector<double> interpolate(double x, Tablulate::INTERPOLATION interp) const;

	std::vector<double> extrapolate(double x, double dx) const;

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

	void set(unsigned idx, const double& x, const std::vector<double>& v)
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

	void add_back(const double& x, const std::vector<double>& v)
	{
		assert(v.size() == NCOLS);

		_x.push_back(x);
		_v.push_back(v);
	}

	std::vector<double> get(double x,
		Tablulate::INTERPOLATION interp = Tablulate::INTERPOLATION::LINEAR,
		Tablulate::EXTRAPOLATION extrap = Tablulate::EXTRAPOLATION::HOLD_LAST_VALUE) const;
};

class XYTable
{
	std::vector<double> _x;

	std::vector<std::array<double,2>> _entries;

public:

	XYTable(std::vector<double> x,
		std::vector<double> y);

	void operator=(const XYTable& table)
	{
		_x = table._x;
		_entries = table._entries;
	}

	double get(double x) const
	{
		auto it = std::lower_bound(_x.begin(),_x.end(), x);
		const auto& entry = _entries[it - _x.begin()];
		return entry[0] + (x - *it) * entry[1];
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

template<typename Float, unsigned NCOLS> 
class BasicTable
{
	std::vector<Float> _x;

	std::vector<std::array<Float, NCOLS>> _v;

	std::vector<std::array<Float, NCOLS>> _dv;

	unsigned (*_get_index)(const std::vector<Float>&, Float);

public:

	static unsigned linear_search(const std::vector<Float>& _x, Float x)
	{
		auto it = std::find_if(_x.begin(), _x.end(), (const auto & a, const auto & b)[] { a > b; });
		idx = std::distance(it, _x.begin());
	}

	static unsigned bisect_search(const std::vector<Float>& _x, Float x)
	{
		auto it = std::lower_bound(_x.begin(), _x.end(), x);
		return std::distance(it, _x.begin());
	}

	BasicTable() {}
	BasicTable(std::vector<Float> x
		std::vector<std::array<Float, NCOLS>> v) :
			_x(std::move(x)),
			_v(std::move(v))
	{
		_dv.resize(_x.size());
		for (auto i = 1u; i < _x.size(); i++)
		{
			auto& row = _dv[i - 1];
			Float dx = _x[i] - _x[i - 1];
			for (auto j = 0u; j < NCOLS; j++)
			{
				row[j] = (_v[i][j] - _v[i - 1][j])/dx;
			}
		}
	}

	~BasicTable() {}

	void clear()
	{
		_x.clear();
		_v.clear();
		_dv.clear();
	}

	void add(const Float x, const std::array<Float, NCOLS>& v)
	{
		if (_x.empty())
		{
			_x.push_back(x);
			_v.push_back(v);
			return;
		}

		std::array<Float, NCOLS> dv;

		if (x > _x.back())
		{

			Float dx = 1.0 / (x - _x.back());
			for (auto i = 0u; i < NCOLS; i++)
			{
				dv[i] = (v[i] - _v.back()[i])*dx;
			}
			_x.push_back(x);
			_v.push_back(v);
			_dv.push_back(dv);
			return;
		}

		unsigned idx = 0;
		while (idx < _x.size())
		{
			if (x < _x[idx])
			{
				break;
			}
			idx++;
		}
		Float dx1 = 1.0 / (x - _x[idx]);
		Float dx2 = 1.0 / (_x[idx + 1]  - x);
		for (auto i = 0u; i < NCOLS; i++)
		{
			dv[i] = (v[i] - _v[idx][i])*dx1;
			_dv[idx][i] = (_v[idx+1][i] - v[i])*dx2;
		}

		_x.insert(_x.begin() + idx, x);
		_v.insert(_v.begin() + idx, v);
		_dv.insert(_dv.begin() + idx, dv);

		if (_x.size() >= 64u)
		{
			_get_index = &bisect_search;
		}
		else
		{
			_get_index = &linear_search;
		}
	}

	void get(const Float x, Float* v) const
	{
		unsigned idx = _get_index(_x, x);
		
		Float dx = x - _x[idx];
		for (auto i = 0u; i < NCOLS; i++)
		{
			v[i] = _v[idx] + _dv[idx] * dx;
		}
	}

};

template<typename Float>
class DynamicTable
{
public:

	static constexpr unsigned NCOEF_CUBIC = 3;

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