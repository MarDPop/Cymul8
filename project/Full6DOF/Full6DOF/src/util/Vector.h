#pragma once

#include <cstring>
#include <algorithm>

template<typename T>
class array_wrapper
{
    T* _data;

public:

    array_wrapper() : 
        _data(nullptr) {}

    array_wrapper( unsigned n ) : 
        _data(new T[n]) {}

    array_wrapper(T* data, unsigned n) :
        _data(new T[n]) 
    {
        memcpy(_data, data, n * sizeof(T));
    }

    array_wrapper(array_wrapper&& copy) :
        _data( std::move(copy._data) ) {}

    array_wrapper(array_wrapper& copy) = delete;

    const T* data() const
    {
        return _data;
    }

    T* data() 
    {
        return _data;
    }

    T& operator[](unsigned idx)
    {
        return _data[idx];
    }

    const T& operator[](unsigned idx) const
    {
        return _data[idx];
    }

    ~array_wrapper()
    {
        delete[] _data;
    }
};

template<typename T>
class fixed_vector
{
    T* const _data;

    const unsigned _size;

public:

    typedef T* iterator;
    typedef const T* const_iterator;

    fixed_vector(unsigned n) :
        _data(new T[n]),
        _size(n)
    {}

    fixed_vector(const fixed_vector& copy) :
        _data(new T[copy._size]),
        _size(copy._size)
    {
        memcpy(_data, copy._data, _size * sizeof(T));
    }

    fixed_vector(fixed_vector&& copy) noexcept :
        _data(std::move(copy._data)),
        _size(copy._size)
    {}

    ~fixed_vector()
    {
        delete[] _data;
    }

    T& operator[](unsigned idx)
    {
        return _data[idx];
    }

    const T& operator[](unsigned idx) const
    {
        return _data[idx];
    }

    void operator=(const fixed_vector& copy)
    {
        memcpy(_data, copy._data, _size * sizeof(T));
    }

    iterator begin() { return _data;}
    const_iterator begin() const { return _data; }
    iterator end() { return _data + _size; }
    const_iterator end() const { return _data + _size; }

    unsigned size() const { return _size; }

    void set(const T* values)
    {
        memcpy(_data, values, _size*sizeof(T));
    }
};

namespace VectorUtil
{

    template<typename T>
    void get_ascending_order(const T* data, unsigned* indices, const unsigned N)
    {
        for (auto i = 0u; i < N; i++)
        {
            indices[i] = i;
        }
        std::sort(indices, indices + N, [&](const auto& i, const auto& j) { return data[i] < data[j]; });
    }

    template<typename T>
    void permute_data(T* data, const unsigned* indices, const unsigned N)
    {
        T* sorted_data = new T[N];
        for (auto i = 0u; i < N; i++)
        {
            sorted_data[i] = data[indices[i]];
        }
        memset(data, sorted_data, N * sizeof(T));
        delete[] sorted_data;
    }

    template<typename T>
    std::vector<unsigned> get_ascending_order(const std::vector<T>& data)
    {
        std::vector<unsigned> indices(data.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](const auto& i, const auto& j) { return data[i] < data[j]; });
        return indices;
    }

}