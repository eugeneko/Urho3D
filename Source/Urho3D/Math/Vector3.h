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

#include "../Math/Vector2.h"
#include "../Math/MathDefs.h"
#include <glm/vec3.hpp>

namespace Urho3D
{

/// Three-dimensional vector with integer values.
using IntVector3 = glm::ivec3;

/// Three-dimensional vector.
using Vector3 = glm::vec3;

/// Per-component linear interpolation between two 3-vectors.
inline Vector3 VectorLerp(const Vector3& lhs, const Vector3& rhs, const Vector3& t) { return lhs + (rhs - lhs) * t; }

/// Per-component min of two 3-vectors.
inline Vector3 VectorMin(const Vector3& lhs, const Vector3& rhs) { return Vector3(Min(lhs.x, rhs.x), Min(lhs.y, rhs.y), Min(lhs.z, rhs.z)); }

/// Per-component max of two 3-vectors.
inline Vector3 VectorMax(const Vector3& lhs, const Vector3& rhs) { return Vector3(Max(lhs.x, rhs.x), Max(lhs.y, rhs.y), Max(lhs.z, rhs.z)); }

/// Per-component floor of 3-vector.
inline Vector3 VectorFloor(const Vector3& vec) { return Vector3(Floor(vec.x), Floor(vec.y), Floor(vec.z)); }

/// Per-component round of 3-vector.
inline Vector3 VectorRound(const Vector3& vec) { return Vector3(Round(vec.x), Round(vec.y), Round(vec.z)); }

/// Per-component ceil of 3-vector.
inline Vector3 VectorCeil(const Vector3& vec) { return Vector3(Ceil(vec.x), Ceil(vec.y), Ceil(vec.z)); }

/// Per-component floor of 3-vector. Returns IntVector3.
inline IntVector3 VectorFloorToInt(const Vector3& vec) { return IntVector3(FloorToInt(vec.x), FloorToInt(vec.y), FloorToInt(vec.z)); }

/// Per-component round of 3-vector. Returns IntVector3.
inline IntVector3 VectorRoundToInt(const Vector3& vec) { return IntVector3(RoundToInt(vec.x), RoundToInt(vec.y), RoundToInt(vec.z)); }

/// Per-component ceil of 3-vector. Returns IntVector3.
inline IntVector3 VectorCeilToInt(const Vector3& vec) { return IntVector3(CeilToInt(vec.x), CeilToInt(vec.y), CeilToInt(vec.z)); }

/// Per-component min of two 3-vectors.
inline IntVector3 VectorMin(const IntVector3& lhs, const IntVector3& rhs) { return IntVector3(Min(lhs.x, rhs.x), Min(lhs.y, rhs.y), Min(lhs.z, rhs.z)); }

/// Per-component max of two 3-vectors.
inline IntVector3 VectorMax(const IntVector3& lhs, const IntVector3& rhs) { return IntVector3(Max(lhs.x, rhs.x), Max(lhs.y, rhs.y), Max(lhs.z, rhs.z)); }

/// Return a random value from [0, 1) from 3-vector seed.
inline float StableRandom(const Vector3& seed) { return StableRandom(Vector2(StableRandom(Vector2(seed.x, seed.y)), seed.z)); }

}
