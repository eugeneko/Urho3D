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

namespace Urho3D
{

/// Three-dimensional vector with integer values.
using IntVector3 = VectorN<int, 3>;

/// Three-dimensional vector.
using Vector3 = VectorN<float, 3>;

/// Per-component floor of 3-vector. Returns IntVector3.
inline IntVector3 VectorFloorToInt(const Vector3& vec) { return IntVector3(FloorToInt(vec.x_), FloorToInt(vec.y_), FloorToInt(vec.z_)); }

/// Per-component round of 3-vector. Returns IntVector3.
inline IntVector3 VectorRoundToInt(const Vector3& vec) { return IntVector3(RoundToInt(vec.x_), RoundToInt(vec.y_), RoundToInt(vec.z_)); }

/// Per-component ceil of 3-vector. Returns IntVector3.
inline IntVector3 VectorCeilToInt(const Vector3& vec) { return IntVector3(CeilToInt(vec.x_), CeilToInt(vec.y_), CeilToInt(vec.z_)); }

/// Return a random value from [0, 1) from 3-vector seed.
inline float StableRandom(const Vector3& seed) { return StableRandom(Vector2(StableRandom(Vector2(seed.x_, seed.y_)), seed.z_)); }

}
