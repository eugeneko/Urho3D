namespace glm
{

#if MIXIN_DIM == 2

template <class T, qualifier Q> const vec<2, T, Q> vec<2, T, Q>::ZERO(0, 0);
template <class T, qualifier Q> const vec<2, T, Q> vec<2, T, Q>::LEFT(-1, 0);
template <class T, qualifier Q> const vec<2, T, Q> vec<2, T, Q>::RIGHT(1, 0);
template <class T, qualifier Q> const vec<2, T, Q> vec<2, T, Q>::UP(0, 1);
template <class T, qualifier Q> const vec<2, T, Q> vec<2, T, Q>::DOWN(0, -1);
template <class T, qualifier Q> const vec<2, T, Q> vec<2, T, Q>::ONE(1, 1);

#elif MIXIN_DIM == 3

template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::ZERO(0, 0, 0);
template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::LEFT(-1, 0, 0);
template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::RIGHT(1, 0, 0);
template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::UP(0, 1, 0);
template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::DOWN(0, -1, 0);
template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::FORWARD(0, 0, 1);
template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::BACK(0, 0, -1);
template <class T, qualifier Q> const vec<3, T, Q> vec<3, T, Q>::ONE(1, 1, 1);

#elif MIXIN_DIM == 4

template <class T, qualifier Q> const vec<4, T, Q> vec<4, T, Q>::ZERO(0, 0, 0, 0);
template <class T, qualifier Q> const vec<4, T, Q> vec<4, T, Q>::ONE(1, 1, 1, 1);

#else
#error Invalid MIXIN_DIM value
#endif

}
