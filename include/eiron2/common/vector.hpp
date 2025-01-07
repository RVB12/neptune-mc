
#ifndef EIRON2_VECTOR_HPP
#define EIRON2_VECTOR_HPP

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <random>

namespace eiron2 {

// Forward declaration of the vectorview class
template <typename T = double>
class vectorview;

/**
 * A vector class that can be used to store vectors of arbitrary size.
 *
 * This class owns the data it points to.
 *
 * @tparam T this is a template parameter that specifies the type of the
 *           elements of the vector.
 */
template <typename T = double>
class vector {
    friend class vectorview<T>;

   public:
    typedef T val_t;

    /**
     * @brief Construct a new vector
     *
     * @param n  integer
     *           Number of elements of the vector
     */
    vector(int n) : _n(n), _data(new T[_n])
    {
#ifndef NDEBUG
        if constexpr (std::is_floating_point<T>::value) {
            for (int i = 0; i < _n; ++i) {
                _data[i] = std::nan("1");
            }
        }
#endif
    }

    /**
     * @brief Construct a new vector and initialize it with a value
     *
     * @param n     integer
     *              Number of elements of the vector
     *
     * @param value T
     *              Value to initialize the vector with
     */
    vector(int n, T value) : _n(n), _data(new T[_n])
    {
        for (int i = 0; i < _n; ++i) {
            _data[i] = value;
        }
    }

    /**
     * @brief Construct a new vector from an initializer list
     */
    vector(std::initializer_list<T> list) : _n(list.size()), _data(new T[_n])
    {
        std::copy(list.begin(), list.end(), _data.get());
    }

    // Destructor
    ~vector()
    {
        // Nothing to do here, the unique pointer will take care of the memory
    }

    /**
     * @brief Construct a new vector as a copy of another vector
     *        This does a deep copy.
     *
     * @param v vector
     *          vector to copy
     */
    vector(const vector& v) : _n(v.size()), _data(new T[_n])
    {
        for (int i = 0; i < _n; ++i) {
            _data[i] = v[i];
        }
    }

    /**
     * @brief Construct a new vector as a copy of another vector
     *        This does a deep copy.
     *
     * This is a templated version, you can safely ignore this until you
     * have learned about expression templates and more advanced C++ features.
     *
     * @param v vector
     *          vector to copy
     */
    vector(const vectorview<T>& v) : _n(v.size()), _data(new T[_n])
    {
        for (int i = 0; i < _n; ++i) {
            _data[i] = v[i];
        }
    }

    /**
     * @brief Construct a new vector object by moving the data of another
     *        vector.
     *        This does a shallow copy.
     *
     * @param v vector
     *          vector to move
     */
    vector(vector&& v) : _n(v.size()), _data(std::move(v._data)) {}

    /**
     * @brief Assign the data of another vector to this vector.
     *        This does a deep copy.
     *
     * @param v vector
     *          vector to copy
     * @return  vector&
     *          Reference to this vector (used for chaining assignments)
     */
    vector& operator=(const vector& v)
    {
        assert(v.size() == _n);
        for (int i = 0; i < _n; ++i) {
            _data[i] = v[i];
        }
        return *this;
    }

    /**
     * @brief Assign the data of another vector to this vector
     *        templated version.
     *        This does a deep copy.
     *
     * @param v vector
     *          vector to copy
     * @return  vector&
     *          Reference to this vector (used for chaining assignments)
     */
    vector& operator=(const vectorview<T>& v)
    {
        assert(v.size() == _n);
        for (int i = 0; i < _n; ++i) {
            _data[i] = v[i];
        }
        return *this;
    }

    /**
     * @brief Assign the data of another vector to this vector
     *        This does a shallow copy
     *
     * @param v vector
     *          vector to move
     * @return  vector&
     *          Reference to this vector (used for chaining assignments)
     */
    vector& operator=(vector&& v)
    {
        assert(v.size() == _n);
        _data = std::move(v._data);
        return *this;
    }

    // Returns the size of the vector
    inline int size() const { return _n; }

    /**
     * @brief Access the i-th element of the vector
     *
     * @param i integer
     *          index
     * @return  const T&
     *          Because this vector is const, we return a const reference.
     */
    inline const T& operator[](const int i) const
    {
        assert(i < _n);
        assert(i >= 0);
        return _data[i];
    }

    /**
     * @brief Access the i-th element of the vector
     *
     * @param i integer
     *          index
     * @return  T&
     */
    inline T& operator[](const int i)
    {
        assert(i < _n);
        assert(i >= 0);
        return _data[i];
    }

    /**
     * Return a subvector of the vector.
     * The subvector starts at index start and ends at index end - 1.
     * The stride is the step size between two elements of the subvector.
     *
     * @param start  first index of the subvector
     * @param end    last index of the subvector (not inclusive)
     * @param stride step size between two elements of the subvector
     *
     * @example
     * vector<float> v(10);
     * auto v2 = v.subvector(1,9,2); // v2 is a length 4 vector containing the
     * elements v[1], v[3], v[5], v[7]
     *
     * Note: this function will also be called if the vector is const, even
     * though the returned object is not const. Solving this problem would
     * require either an inefficient or a more complex design.
     */
    inline vectorview<T> subvector(int start, int end, int stride = 1) const
    {
        assert(start >= 0);
        assert(end <= _n);
        assert(start < end);
        assert(stride > 0);
        return vectorview<T>((end - start) / stride, &_data[start], stride);
    }

    /**
     * Return a pointer to the data of the vector.
     */
    inline const T* data() const { return _data.get(); }

    /**
     * Return a pointer to the data of the vector.
     */
    inline T* data() { return _data.get(); }

    /**
     * Return a pointer to the beginning of the vector.
     */
    inline T* begin() { return _data.get(); }

    /**
     * Return a pointer to the end of the vector.
     */
    inline T* end() { return _data.get() + _n; }

    /**
     * Return a const pointer to the beginning of the vector.
     */
    inline const T* cbegin() const { return _data.get(); }

    /**
     * Return a const pointer to the end of the vector.
     */
    inline const T* cend() const { return _data.get() + _n; }

   private:
    const int _n;
    std::unique_ptr<T[]> _data;
};

/**
 * A vector class that can be used to store vectors of arbitrary size.
 *
 * This class does not own the data it points to.
 *
 * @tparam T this is a template parameter that specifies the type of the
 *           elements of the vector.
 */
template <typename T>
class vectorview {
   public:
    typedef T val_t;

    // Constructor
    vectorview()
        : _n(0), _stride(1), _data(nullptr)
    {}

    vectorview(int n, T* data, int stride = 1)
        : _n(n), _stride(stride), _data(data)
    {}

    // Destructor
    ~vectorview()
    {
        // Nothing to do here, the view does not own the data
    }

    /**
     * @brief Construct a new vectorview as a copy of another vectorview
     *        This does a shallow copy.
     *
     * @param v vector
     *          vector to copy
     */
    vectorview(const vectorview& v)
        : _n(v._n), _stride(v._stride), _data(v._data)
    {}

    /**
     * @brief Construct a new vectorview to view the data of a vector
     *        This does a shallow copy.
     *
     */
    vectorview(const vector<T>& v)
        : _n(v.size()), _stride(1), _data(v._data.get())
    {}

    /**
     * @brief Construct a new vectorview object by moving the data of another
     *        vector. This does a shallow copy.
     *
     * @param v vector
     *          vector to move
     */
    vectorview(vectorview&& v)
        : _n(v._n), _stride(v._stride), _data(v._data)
    {}

    /**
     * @brief Assign the data of another vector to this vector
     *        This does a shallow copy.
     *
     * @param v vectorview
     *          vectorview to copy
     * @return  vectorview&
     *          Reference to this vector (used for chaining assignments)
     */
    vectorview& operator=(const vectorview& v)
    {
        _n = v._n;
        _stride = v._stride;
        _data = v._data;
        return *this;
    }

    /**
     * @brief Assign the data of another vector to this vector
     *        This does a shallow copy
     *
     * @param v vector
     *          vector to move
     * @return  vector&
     *          Reference to this vector (used for chaining assignments)
     */
    vectorview& operator=(vectorview&& v)
    {
        _n = v._n;
        _stride = v._stride;
        _data = v._data;
        return *this;
    }

    // Returns the size of the vector
    inline int size() const { return _n; }

    // Returns the stride of the vector
    inline int stride() const { return _stride; }

    /**
     * @brief Access the i-th element of the vector
     *
     * @param i integer
     *          index
     * @return  const T&
     *          Because this vector is const, we return a const reference.
     */
    inline const T& operator[](const int i) const
    {
        assert(i < _n);
        assert(i >= 0);
        return _data[i * _stride];
    }

    /**
     * @brief Access the i-th element of the vector
     *
     * @param i integer
     *          index
     * @return  T&
     */
    inline T& operator[](const int i)
    {
        assert(i < _n);
        assert(i >= 0);
        return _data[i * _stride];
    }

    /**
     * Return a subvector of the vector.
     * The subvector starts at index start and ends at index end - 1.
     * The stride is the step size between two elements of the subvector.
     *
     * @param start  first index of the subvector
     * @param end    last index of the subvector (not inclusive)
     * @param stride step size between two elements of the subvector
     *
     * @example
     * vector<float> v(10);
     * auto v2 = v.subvector(1,9,2); // v2 is a length 4 vector containing the
     * elements v[1], v[3], v[5], v[7]
     *
     * Note: this function will also be called if the vector is const, even
     * though the returned object is not const. Solving this problem would
     * require either an inefficient or a more complex design.
     */
    inline vectorview subvector(int start, int end, int stride = 1) const
    {
        assert(start >= 0);
        assert(end <= _n);
        assert(start < end);
        assert(stride > 0);
        return vectorview((end - start) / stride, &_data[start * _stride], _stride * stride);
    }

    /**
     * Return a pointer to the data of the vector.
     */
    inline const T* data() const { return _data; }

    /**
     * Return a pointer to the data of the vector.
     */
    inline T* data() { return _data; }

   private:
    int _n;
    int _stride;
    T* _data;
};

}  // namespace tws

#endif