//
// Copyright (c) 2008-2017 the Urho3D project.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

#pragma once

#include "../Container/Str.h"
#include "../Math/MathDefs.h"

namespace Urho3D
{

template <class T, unsigned N> class VectorN
{
    /// Normalize to unit length.
    void Normalize();
    /// Return length.
    /// \todo Return double if T is double
    float Length() const;
    /// Return squared length.
    T LengthSquared() const;
    /// Calculate dot product.
    T DotProduct(const VectorN& rhs) const;
    /// Calculate absolute dot product.
    T AbsDotProduct(const VectorN& rhs) const;
    /// Project vector onto axis.
    /// \todo Return double if T is double
    float ProjectOntoAxis(const VectorN& axis) const;
    /// Return absolute vector.
    VectorN Abs() const;
    /// Linear interpolation with another vector.
    VectorN Lerp(const VectorN& rhs, float t) const;
    /// Test for equality with another vector with epsilon.
    bool Equals(const VectorN& rhs) const;
    /// Returns the angle between this vector and another vector in degrees.
    /// \todo Return double if T is double
    float Angle(const VectorN& rhs) const;
    /// Return whether is NaN.
    bool IsNaN() const;
    /// Return normalized to unit length.
    VectorN Normalized() const;
    /// Return as string.
    String ToString() const;
    /// Return hash value for HashSet & HashMap.
    unsigned ToHash() const;

};

template <class T> class VectorN<T, 2>
{
public:
    /// Construct zero.
    VectorN() = default;
    /// Copy construct.
    VectorN(const VectorN& vector) = default;
    /// Copy construct from another vector type.
    template <class U> explicit VectorN(const VectorN<U, 2>& vector) :
        x_(static_cast<T>(vector.x_)),
        y_(static_cast<T>(vector.y_))
    {
    }
    /// Construct from coordinates.
    VectorN(T x, T y) : x_(x), y_(y) { }
    /// Construct from array.
    explicit VectorN(const T* data) : x_(data[0]), y_(data[1]) { }
    /// Assign from another vector.
    VectorN& operator =(const VectorN& rhs) = default;

    /// Return vector data.
    const T* Data() const { return &x_; }
    /// Return vector element by index.
    T& operator[](unsigned idx) { return (&x_)[idx]; }
    /// Return const vector element by index.
    T operator[](unsigned idx) const { return (&x_)[idx]; }

    /// X coordinate.
    T x_ = T{};
    /// Y coordinate.
    T y_ = T{};

    /// Zero vector.
    static const VectorN ZERO;
    /// (-1,0) vector.
    static const VectorN LEFT;
    /// (1,0) vector.
    static const VectorN RIGHT;
    /// (0,1) vector.
    static const VectorN UP;
    /// (0,-1) vector.
    static const VectorN DOWN;
    /// (1,1) vector.
    static const VectorN ONE;

public:
    // Copy-paste from VectorN<T, N>
    void Normalize();
    float Length() const;
    T LengthSquared() const;
    T DotProduct(const VectorN& rhs) const;
    T AbsDotProduct(const VectorN& rhs) const;
    float ProjectOntoAxis(const VectorN& axis) const;
    VectorN Abs() const;
    VectorN Lerp(const VectorN& rhs, float t) const;
    bool Equals(const VectorN& rhs) const;
    float Angle(const VectorN& rhs) const;
    bool IsNaN() const;
    VectorN Normalized() const;
    String ToString() const;
    unsigned ToHash() const;

};

template <class T> class VectorN<T, 3>
{
public:
    /// Construct zero.
    VectorN() = default;
    /// Copy construct.
    VectorN(const VectorN& vector) = default;
    /// Copy construct from another vector type.
    template <class U> explicit VectorN(const VectorN<U, 3>& vector) :
        x_(static_cast<T>(vector.x_)),
        y_(static_cast<T>(vector.y_)),
        z_(static_cast<T>(vector.z_))
    {
    }
    /// Construct from a two-dimensional vector and the Z coordinate.
    VectorN(const VectorN<T, 2>& vector, T z) : x_(vector.x_), y_(vector.y_), z_(z) { }
    /// Construct from a two-dimensional vector (for Urho2D).
    VectorN(const VectorN<T, 2>& vector) : x_(vector.x_), y_(vector.y_), z_(0) { }
    /// Construct from two-dimensional coordinates (for Urho2D).
    VectorN(T x, T y) : x_(x), y_(y), z_(0) { }
    /// Construct from coordinates.
    VectorN(T x, T y, T z) : x_(x), y_(y), z_(z) { }
    /// Construct from array.
    explicit VectorN(const T* data) : x_(data[0]), y_(data[1]), z_(data[2]) { }
    /// Assign from another vector.
    VectorN& operator =(const VectorN& rhs) = default;

    /// Return vector data.
    const T* Data() const { return &x_; }
    /// Return vector element by index.
    T& operator[](unsigned idx) { return (&x_)[idx]; }
    /// Return const vector element by index.
    T operator[](unsigned idx) const { return (&x_)[idx]; }

    /// Make vector orthogonal to the axis.
    VectorN Orthogonalize(const VectorN& axis) const { return axis.CrossProduct(*this).CrossProduct(axis).Normalized(); }
    /// Calculate cross product.
    VectorN CrossProduct(const VectorN& rhs) const
    {
        return VectorN(
            y_ * rhs.z_ - z_ * rhs.y_,
            z_ * rhs.x_ - x_ * rhs.z_,
            x_ * rhs.y_ - y_ * rhs.x_
        );
    }

    /// X coordinate.
    T x_ = T{};
    /// Y coordinate.
    T y_ = T{};
    /// Z coordinate.
    T z_ = T{};

    /// Zero vector.
    static const VectorN ZERO;
    /// (-1,0,0) vector.
    static const VectorN LEFT;
    /// (1,0,0) vector.
    static const VectorN RIGHT;
    /// (0,1,0) vector.
    static const VectorN UP;
    /// (0,-1,0) vector.
    static const VectorN DOWN;
    /// (0,0,1) vector.
    static const VectorN FORWARD;
    /// (0,0,-1) vector.
    static const VectorN BACK;
    /// (1,1,1) vector.
    static const VectorN ONE;

public:
    // Copy-paste from VectorN<T, N>
    void Normalize();
    float Length() const;
    T LengthSquared() const;
    T DotProduct(const VectorN& rhs) const;
    T AbsDotProduct(const VectorN& rhs) const;
    float ProjectOntoAxis(const VectorN& axis) const;
    VectorN Abs() const;
    VectorN Lerp(const VectorN& rhs, float t) const;
    bool Equals(const VectorN& rhs) const;
    float Angle(const VectorN& rhs) const;
    bool IsNaN() const;
    VectorN Normalized() const;
    String ToString() const;
    unsigned ToHash() const;

};

template <class T> class VectorN<T, 4>
{
public:
    /// Construct zero.
    VectorN() = default;
    /// Copy construct.
    VectorN(const VectorN& vector) = default;
    /// Copy construct from another vector type.
    template <class U> explicit VectorN(const VectorN<U, 4>& vector) :
        x_(static_cast<T>(vector.x_)),
        y_(static_cast<T>(vector.y_)),
        z_(static_cast<T>(vector.z_)),
        w_(static_cast<T>(vector.w_))
    {
    }
    /// Construct from a 3-dimensional vector and the W coordinate.
    VectorN(const VectorN<T, 3>& vector, T w) :
        x_(vector.x_),
        y_(vector.y_),
        z_(vector.z_),
        w_(w)
    {
    }
    /// Construct from coordinates.
    VectorN(T x, T y, T z, T w) : x_(x), y_(y), z_(z), w_(w) { }
    /// Construct from a float array.
    explicit VectorN(const T* data) : x_(data[0]), y_(data[1]), z_(data[2]), w_(data[3]) { }
    /// Assign from another vector.
    VectorN& operator =(const VectorN& rhs) = default;

    /// Return vector data.
    const T* Data() const { return &x_; }
    /// Return vector element by index.
    T& operator[](unsigned idx) { return (&x_)[idx]; }
    /// Return const vector element by index.
    T operator[](unsigned idx) const { return (&x_)[idx]; }

    /// Make vector orthogonal to the axis.
    VectorN Orthogonalize(const VectorN& axis) const { return axis.CrossProduct(*this).CrossProduct(axis).Normalized(); }
    /// Calculate cross product.
    VectorN CrossProduct(const VectorN& rhs) const
    {
        return VectorN(
            y_ * rhs.z_ - z_ * rhs.y_,
            z_ * rhs.x_ - x_ * rhs.z_,
            x_ * rhs.y_ - y_ * rhs.x_
        );
    }

    /// X coordinate.
    T x_ = T{};
    /// Y coordinate.
    T y_ = T{};
    /// Z coordinate.
    T z_ = T{};
    /// W coordinate.
    T w_ = T{};

    /// Zero vector.
    static const VectorN ZERO;
    /// (1,1,1) vector.
    static const VectorN ONE;

public:
    // Copy-paste from VectorN<T, N>
    void Normalize();
    float Length() const;
    T LengthSquared() const;
    T DotProduct(const VectorN& rhs) const;
    T AbsDotProduct(const VectorN& rhs) const;
    float ProjectOntoAxis(const VectorN& axis) const;
    VectorN Abs() const;
    VectorN Lerp(const VectorN& rhs, float t) const;
    bool Equals(const VectorN& rhs) const;
    float Angle(const VectorN& rhs) const;
    bool IsNaN() const;
    VectorN Normalized() const;
    String ToString() const;
    unsigned ToHash() const;

};

/// Test for equality with another vector without epsilon.
template <class T, unsigned N> bool operator ==(const VectorN<T, N>& lhs, const VectorN<T, N>& rhs)
{
    for (unsigned i = 0; i < N; ++i)
        if (lhs[i] != rhs[i])
            return false;
    return true;
}

/// Test for inequality with another vector without epsilon.
template <class T, unsigned N> bool operator !=(const VectorN<T, N>& lhs, const VectorN<T, N>& rhs)
{
    return !(lhs == rhs);
}

/// Add-assign a vector.
template <class T, unsigned N> VectorN<T, N>& operator +=(VectorN<T, N>& lhs, const VectorN<T, N>& rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] += rhs[i];
    return lhs;
}

/// Subtract-assign a vector.
template <class T, unsigned N> VectorN<T, N>& operator -=(VectorN<T, N>& lhs, const VectorN<T, N>& rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] -= rhs[i];
    return lhs;
}

/// Multiply-assign a scalar.
template <class T, unsigned N, class U> VectorN<T, N>& operator *=(VectorN<T, N>& lhs, U rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] *= static_cast<T>(rhs);
    return lhs;
}

/// Multiply-assign a vector.
template <class T, unsigned N> VectorN<T, N>& operator *=(VectorN<T, N>& lhs, const VectorN<T, N>& rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] *= rhs[i];
    return lhs;
}

/// Divide-assign a scalar.
template <class T, unsigned N, class U> VectorN<T, N>& operator /=(VectorN<T, N>& lhs, U rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] /= static_cast<T>(rhs);
    return lhs;
}

/// Divide-assign a vector.
template <class T, unsigned N> VectorN<T, N>& operator /=(VectorN<T, N>& lhs, const VectorN<T, N>& rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] /= rhs[i];
    return lhs;
}

/// Return negation.
template <class T, unsigned N> VectorN<T, N> operator -(VectorN<T, N> lhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] = -lhs[i];
    return lhs;
}

/// Add a vector.
template <class T, unsigned N> VectorN<T, N> operator +(VectorN<T, N> lhs, const VectorN<T, N>& rhs)
{
    lhs += rhs;
    return lhs;
}

/// Subtract a vector.
template <class T, unsigned N> VectorN<T, N> operator -(VectorN<T, N> lhs, const VectorN<T, N>& rhs)
{
    lhs -= rhs;
    return lhs;
}

/// Multiply with a scalar.
template <class T, unsigned N, class U> VectorN<T, N> operator *(VectorN<T, N> lhs, U rhs)
{
    lhs *= rhs;
    return lhs;
}

/// Multiply with a scalar.
template <class T, unsigned N, class U> VectorN<T, N> operator *(U lhs, VectorN<T, N> rhs)
{
    rhs *= lhs;
    return rhs;
}

/// Multiply with a vector.
template <class T, unsigned N> VectorN<T, N> operator *(VectorN<T, N> lhs, const VectorN<T, N>& rhs)
{
    lhs *= rhs;
    return lhs;
}

/// Divide by a scalar.
template <class T, unsigned N, class U> VectorN<T, N> operator /(VectorN<T, N> lhs, U rhs)
{
    lhs /= rhs;
    return lhs;
}

/// Divide by a vector.
template <class T, unsigned N> VectorN<T, N> operator /(VectorN<T, N> lhs, const VectorN<T, N>& rhs)
{
    lhs /= rhs;
    return lhs;
}

/// Per-component linear interpolation between two 2-vectors.
template <class T, unsigned N> VectorN<T, N> VectorLerp(const VectorN<T, N>& lhs, const VectorN<T, N>& rhs, const VectorN<T, N>& t)
{
    return lhs + (rhs - lhs) * t;
}

/// Per-component min of two 2-vectors.
template <class T, unsigned N> VectorN<T, N> VectorMin(VectorN<T, N> lhs, const VectorN<T, N>& rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] = Min(lhs[i], rhs[i]);
    return lhs;
}

/// Per-component max of two 2-vectors.
template <class T, unsigned N> VectorN<T, N> VectorMax(VectorN<T, N> lhs, const VectorN<T, N>& rhs)
{
    for (unsigned i = 0; i < N; ++i)
        lhs[i] = Max(lhs[i], rhs[i]);
    return lhs;
}

/// Per-component floor of 2-vector.
template <class T, unsigned N> VectorN<T, N> VectorFloor(VectorN<T, N> vec)
{
    for (unsigned i = 0; i < N; ++i)
        vec[i] = Floor(vec[i]);
    return vec;
}

/// Per-component round of 2-vector.
template <class T, unsigned N> VectorN<T, N> VectorRound(VectorN<T, N> vec)
{
    for (unsigned i = 0; i < N; ++i)
        vec[i] = Round(vec[i]);
    return vec;
}

/// Per-component ceil of 2-vector.
template <class T, unsigned N> VectorN<T, N> VectorCeil(VectorN<T, N> vec)
{
    for (unsigned i = 0; i < N; ++i)
        vec[i] = Ceil(vec[i]);
    return vec;
}

// Implement methods
/// Normalize to unit length.
template <class T, unsigned N> void VectorN<T, N>::Normalize()
{
    float lenSquared = LengthSquared();
    if (!Urho3D::Equals(lenSquared, 1.0f) && lenSquared > 0.0f)
    {
        float invLen = 1.0f / sqrtf(lenSquared);
        for (unsigned i = 0; i < N; ++i)
            (*this[i]) *= invLen;
    }
}

/// Return length.
template <class T, unsigned N> float VectorN<T, N>::Length() const { return sqrtf(LengthSquared()); }

/// Return squared length.
template <class T, unsigned N> T VectorN<T, N>::LengthSquared() const { return DotProduct(*this); }

/// Calculate dot product.
template <class T, unsigned N> T VectorN<T, N>::DotProduct(const VectorN<T, N>& rhs) const
{
    T dp = 0;
    for (unsigned i = 0; i < N; ++i)
        dp += (*this)[i] * rhs[i];
    return dp;
}

/// Calculate absolute dot product.
template <class T, unsigned N> T VectorN<T, N>::AbsDotProduct(const VectorN<T, N>& rhs) const
{
    T dp = 0;
    for (unsigned i = 0; i < N; ++i)
        dp += Abs((*this)[i] * rhs[i]);
    return dp;
}

/// Project vector onto axis.
template <class T, unsigned N> float VectorN<T, N>::ProjectOntoAxis(const VectorN<T, N>& axis) const
{
    return DotProduct(axis.Normalized());
}

/// Return absolute vector.
template <class T, unsigned N> VectorN<T, N> VectorN<T, N>::Abs() const
{
    VectorN<T, N> vec;
    for (unsigned i = 0; i < N; ++i)
        vec[i] = Abs((*this)[i]);
    return vec;
}

/// Linear interpolation with another vector.
template <class T, unsigned N> VectorN<T, N> VectorN<T, N>::Lerp(const VectorN<T, N>& rhs, float t) const
{
    return *this * (1.0f - t) + rhs * t;
}

/// Test for equality with another vector with epsilon.
template <class T, unsigned N> bool VectorN<T, N>::Equals(const VectorN<T, N>& rhs) const
{
    for (unsigned i = 0; i < N; ++i)
        if (!Urho3D::Equals((*this)[i], rhs[i]))
            return false;
    return true;
}

/// Returns the angle between this vector and another vector in degrees.
template <class T, unsigned N> float VectorN<T, N>::Angle(const VectorN<T, N>& rhs) const { return Urho3D::Acos(DotProduct(rhs) / (Length() * rhs.Length())); }

/// Return whether is NaN.
template <class T, unsigned N> bool VectorN<T, N>::IsNaN() const
{
    for (unsigned i = 0; i < N; ++i)
        if (IsNaN((*this)[i]))
            return true;
    return false;
}

/// Return normalized to unit length.
template <class T, unsigned N> VectorN<T, N> VectorN<T, N>::Normalized() const
{
    VectorN<T, N> vec = *this;
    vec.Normalize();
    return vec;
}

// Define constants
template <class T> const VectorN<T, 2> VectorN<T, 2>::ZERO(0, 0);
template <class T> const VectorN<T, 2> VectorN<T, 2>::LEFT(-1, 0);
template <class T> const VectorN<T, 2> VectorN<T, 2>::RIGHT(1, 0);
template <class T> const VectorN<T, 2> VectorN<T, 2>::UP(0, 1);
template <class T> const VectorN<T, 2> VectorN<T, 2>::DOWN(0, -1);
template <class T> const VectorN<T, 2> VectorN<T, 2>::ONE(1, 1);

template <class T> const VectorN<T, 3> VectorN<T, 3>::ZERO;
template <class T> const VectorN<T, 3> VectorN<T, 3>::LEFT(-1, 0, 0);
template <class T> const VectorN<T, 3> VectorN<T, 3>::RIGHT(1, 0, 0);
template <class T> const VectorN<T, 3> VectorN<T, 3>::UP(0, 1, 0);
template <class T> const VectorN<T, 3> VectorN<T, 3>::DOWN(0, -1, 0);
template <class T> const VectorN<T, 3> VectorN<T, 3>::FORWARD(0, 0, 1);
template <class T> const VectorN<T, 3> VectorN<T, 3>::BACK(0, 0, -1);
template <class T> const VectorN<T, 3> VectorN<T, 3>::ONE(1, 1, 1);

template <class T> const VectorN<T, 4> VectorN<T, 4>::ZERO;
template <class T> const VectorN<T, 4> VectorN<T, 4>::ONE(1, 1, 1, 1);

// Specializations
template <> String VectorN<int, 2>::ToString() const;
template <> String VectorN<int, 3>::ToString() const;
template <> String VectorN<int, 4>::ToString() const;
template <> String VectorN<float, 2>::ToString() const;
template <> String VectorN<float, 3>::ToString() const;
template <> String VectorN<float, 4>::ToString() const;

template <> unsigned VectorN<int, 2>::ToHash() const;
template <> unsigned VectorN<int, 3>::ToHash() const;
template <> unsigned VectorN<int, 4>::ToHash() const;
template <> unsigned VectorN<float, 2>::ToHash() const;
template <> unsigned VectorN<float, 3>::ToHash() const;
template <> unsigned VectorN<float, 4>::ToHash() const;

/// Two-dimensional vector with integer values.
using IntVector2 = VectorN<int, 2>;
/// Two-dimensional vector.
using Vector2 = VectorN<float, 2>;

/// Per-component floor of 2-vector. Returns IntVector2.
inline IntVector2 VectorFloorToInt(const Vector2& vec) { return IntVector2(FloorToInt(vec.x_), FloorToInt(vec.y_)); }

/// Per-component round of 2-vector. Returns IntVector2.
inline IntVector2 VectorRoundToInt(const Vector2& vec) { return IntVector2(RoundToInt(vec.x_), RoundToInt(vec.y_)); }

/// Per-component ceil of 2-vector. Returns IntVector2.
inline IntVector2 VectorCeilToInt(const Vector2& vec) { return IntVector2(CeilToInt(vec.x_), CeilToInt(vec.y_)); }

/// Return a random value from [0, 1) from 2-vector seed.
/// http://stackoverflow.com/questions/12964279/whats-the-origin-of-this-glsl-rand-one-liner
inline float StableRandom(const Vector2& seed) { return Fract(Sin(seed.DotProduct(Vector2(12.9898f, 78.233f)) * M_RADTODEG) * 43758.5453f); }

/// Return a random value from [0, 1) from scalar seed.
inline float StableRandom(float seed) { return StableRandom(Vector2(seed, seed)); }

}
