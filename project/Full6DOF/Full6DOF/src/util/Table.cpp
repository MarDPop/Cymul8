#include "Table.h"

#include "Vector.h"
#include <assert.h>

std::vector<double> Table::interpolate(double x, Tabulate::INTERPOLATION interp) const
{
	auto it = std::lower_bound(_x.begin(), _x.end(), x);
	auto idx = std::distance(_x.begin(), it);
	std::vector<double> p1(NCOLS);
	if (static_cast<int>(interp) > 1)
	{
		p1 = _v[idx + static_cast<int>(interp)];
	}
	else
	{
		auto dx = *(it + 1) - *it;
		auto delta = x - *it;
		if (interp == Tabulate::INTERPOLATION::NEAREST)
		{
			p1 = _v[idx + (static_cast<double>(2.0) * delta > dx)];
		}
		else
		{
			delta /= dx;
			p1 = _v[idx];
			const auto& p2 = _v[idx + 1];
			switch (interp)
			{
			case Tabulate::INTERPOLATION::LINEAR:
				for (auto i = 0u; i < NCOLS; i++)
				{
					p1[i] += (p2[i] - p1[i]) * delta;
				}
				break;
			case Tabulate::INTERPOLATION::CUBIC:
				const auto& p0 = _v[idx - (idx > 0)];
				const auto& p3 = _v[idx + 1 + (idx < (_x.size() - 2))];
				const auto halfx = 0.5 * delta;
				for (auto j = 0u; j < NCOLS; j++)
				{
					p1[j] += halfx * (p2[j] - p0[j] +
						delta * (2.0 * p0[j] - 5.0 * p1[j] + 4.0 * p2[j] - p3[j] +
							delta * (3.0 * (p1[j] - p2[j]) + p3[j] - p0[j])));
				}
				break;
			}

		}

	}
	return p1;
}

std::vector<double> Table::extrapolate(double x, double dx) const
{
	unsigned idx = (dx > 0) * (_x.size() - 2);
	dx /= (_x[idx + 1] - _x[idx]);
	std::vector<double> p1 = _v[idx];
	const auto& p2 = _v[idx + 1];
	for (auto i = 0u; i < NCOLS; i++)
	{
		p1[i] += (p2[i] - p1[i]) * dx;
	}
	return p1;
}

std::vector<double> Table::get(double x,
	Tabulate::INTERPOLATION interp,
	Tabulate::EXTRAPOLATION extrap) const
{
	if (x < _x[0])
	{
		if (extrap == Tabulate::EXTRAPOLATION::HOLD_LAST_VALUE)
		{
			return _v[0];
		}
		return this->extrapolate(x, x - _x[0]);
	}
	if (x > _x.back())
	{
		if (extrap == Tabulate::EXTRAPOLATION::HOLD_LAST_VALUE)
		{
			return _v.back();
		}
		return this->extrapolate(x, x - _x.back());
	}
	return this->interpolate(x, interp);
}

XYTable::XYTable(std::vector<double> x,
	std::vector<double> y) :
	_x(std::move(x))
{
	assert(_x.size() == y.size());

	_entries.resize(_x.size());

	if (!std::is_sorted(_x.begin(), _x.end()))
	{
		auto idx = VectorUtil::get_ascending_order(_x);
		VectorUtil::permute_data(_x.data(), idx.data(), _x.size());
		VectorUtil::permute_data(y.data(), idx.data(), _x.size());
	}
	
	for (auto i = 0u; i < (_x.size() - 1); i++)
	{
		_entries[i][0] = y[i];
		_entries[i][1] = (y[i+1] - y[i]) /( _x[i + 1] - _x[i]);
	}
	_entries.back()[0] = y.back();
}