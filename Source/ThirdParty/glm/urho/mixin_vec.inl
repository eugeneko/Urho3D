/// Construct from data.
static vec FromData(const T* data)
{
    vec result;
    for (length_t i = 0; i < length(); ++i)
        result[i] = data[i];
    return result;
}

#if MIXIN_DIM == 3

/// Construct from 2D vector.
vec(const vec<2, T, Q>& rhs) : x(rhs.x), y(rhs.y), z(0) { }

#endif

/// Calculate dot product.
T DotProduct(const vec& rhs) const
{
    T dp = 0;
    for (length_t i = 0; i < length(); ++i)
        dp += (*this)[i] * rhs[i];
    return dp;
}

/// Return squared length.
T LengthSquared() const { return DotProduct(*this); }

/// Return length.
float Length() const { return sqrtf(LengthSquared()); }

/// Normalize to unit length.
void Normalize()
{
    const float lenSquared = static_cast<float>(LengthSquared());
    if (!Urho3D::Equals(lenSquared, 1.0f) && lenSquared > 0.0f)
        (*this) *= 1.0f / sqrtf(lenSquared);
}

/// Return normalized to unit length.
vec Normalized() const
{
    vec result = *this;
    result.Normalize();
    return result;
}

/// Calculate absolute dot product.
float AbsDotProduct(const vec& rhs) const
{
    T dp = 0;
    for (length_t i = 0; i < length(); ++i)
        dp += Urho3D::Abs((*this)[i] * rhs[i]);
    return dp;
}

/// Project vector onto axis.
float ProjectOntoAxis(const vec& axis) const { return DotProduct(axis.Normalized()); }

#if MIXIN_DIM == 3
/// Make vector orthogonal to the axis.
vec Orthogonalize(const vec& axis) const { return axis.CrossProduct(*this).CrossProduct(axis).Normalized(); }

/// Calculate cross product.
vec CrossProduct(const vec& rhs) const
{
    return vec(
        y * rhs.z - z * rhs.y,
        z * rhs.x - x * rhs.z,
        x * rhs.y - y * rhs.x
    );
}
#endif

/// Return absolute vector.
vec Abs() const
{
    vec result;
    for (length_t i = 0; i < length(); ++i)
        result[i] = Urho3D::Abs((*this)[i]);
    return result;
}

/// Linear interpolation with another vector.
vec Lerp(const vec& rhs, float t) const { return *this * (1.0f - t) + rhs * t; }

/// Test for equality with another vector with epsilon.
bool Equals(const vec& rhs) const
{
    for (length_t i = 0; i < length(); ++i)
        if (!Urho3D::Equals((*this)[i], rhs[i]))
            return false;
    return true;
}

/// Returns the angle between this vector and another vector in degrees.
float Angle(const vec& rhs) const { return Urho3D::Acos(DotProduct(rhs) / (Length() * rhs.Length())); }

/// Return whether is NaN.
bool IsNaN() const
{
    for (length_t i = 0; i < N; ++i)
        if (IsNaN(vec[i]))
            return true;
    return false;
}

/// Return float data.
const T* Data() const { return &x; }

/// Return as string.
Urho3D::String ToString() const;

/// Return hash value for HashSet & HashMap.
unsigned ToHash() const
{
    // \todo Implement me
    return 0;
//     unsigned hash = 37;
//     hash = 37 * hash + FloatToRawIntBits(x);
//     hash = 37 * hash + FloatToRawIntBits(y);
//     hash = 37 * hash + FloatToRawIntBits(z);
//
//     return hash;
}

/// Zero vector.
static const vec ZERO;
/// One vector.
static const vec ONE;

#if MIXIN_DIM == 2 || MIXIN_DIM == 3
/// (-1,0) vector.
static const vec LEFT;
/// (1,0) vector.
static const vec RIGHT;
/// (0,1) vector.
static const vec UP;
/// (0,-1) vector.
static const vec DOWN;
#endif

#if MIXIN_DIM == 3
/// (0,0,1) vector.
static const vec FORWARD;
/// (0,0,-1) vector.
static const vec BACK;
#endif
