#pragma once

#include <cstring>

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

    const T* data() const
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

    inline fixed_vector(unsigned n) :
        _data(new T[n]),
        _size(n)
    {}

    inline fixed_vector(const fixed_vector& copy) :
        _data(new T[copy._size]),
        _size(copy._size)
    {
        memset(_data, copy._data, _size * sizeof(T));
    }

    inline fixed_vector(fixed_vector&& copy) :
        _data(std::move(copy._data)),
        _size(copy._size)
    {
        memset(_data, copy._data, _size * sizeof(T));
    }

    inline ~fixed_vector()
    {
        delete[] _data;
    }

    inline T& operator[](unsigned idx)
    {
        return _data[idx];
    }

    inline const T& operator[](unsigned idx) const
    {
        return _data[idx];
    }

    inline void operator=(const fixed_vector& copy)
    {
        memset(_data, copy._data, _size * sizeof(T));
    }

    inline iterator begin() { return _data;}
    inline const_iterator begin() const { return _data; }
    inline iterator end() { return _data + _size; }
    inline const_iterator end() const { return _data + _size; }

    inline unsigned size() const { return _size; }

    inline void set(const T* values)
    {
        memset(_data, values, _size * sizeof(T));
    }
};