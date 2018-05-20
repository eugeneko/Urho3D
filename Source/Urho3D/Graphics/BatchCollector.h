//
// Copyright (c) 2008-2018 the Urho3D project.
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

#include "../Core/Object.h"
#include "../Graphics/Batch.h"
#include "../Container/HashMap.h"
// #include "../Container/List.h"
// #include "../Graphics/Batch.h"
// #include "../Graphics/Light.h"
// #include "../Graphics/Zone.h"
// #include "../Math/Polyhedron.h"

namespace Urho3D
{

struct ScenePassInfo;
class WorkQueue;

// TODO(eugeneko) This is stub, replace after SceneGrid integration
struct SceneGridDrawableSoA
{
    unsigned size_{};
    Vector<Drawable*> drawable_;
    Vector<StringHash> drawableType_;
    Vector<bool> visible_;
    Vector<unsigned> numLights_;
    Vector<unsigned> firstLight_;

    bool IsValid() const
    {
        return size_ == drawable_.Size()
            && size_ == drawableType_.Size()
            && size_ == visible_.Size()
            && size_ == numLights_.Size()
            && size_ == firstLight_.Size()
            ;
    }
    void ClearVisible()
    {
        for (unsigned i = 0; i < size_; ++i)
            visible_[i] = false;
    }
    void ClearTemporary()
    {
        for (unsigned i = 0; i < size_; ++i)
        {
            numLights_[i] = 0;
            firstLight_[i] = 0;
        }
    }
    void Push(Drawable* drawable)
    {
        ++size_;

        drawable_.Push(drawable);
        drawableType_.Push(drawable->GetType());
        visible_.Push(false);
        numLights_.Push(0);
        firstLight_.Push(0);

        assert(IsValid());
    }
    void EraseSwap(unsigned index)
    {
        assert(0);
    }
    void Clear()
    {
        assert(0);
    }

};

struct LitGeometryDescIdx
{
    unsigned drawableIndex_{};
    float sortValue_{};
    unsigned short lightIndex_{};
    bool perVertex_{};
    bool negativeLight_{};
};

struct LitGeometryDescPacked
{
    LitGeometryDescPacked() = default;
    LitGeometryDescPacked(const LitGeometryDescIdx& source)
        : sortValue_(source.sortValue_)
        , lightIndex_(source.lightIndex_)
        , perVertex_(source.perVertex_)
        , negativeLight_(source.negativeLight_)
    {
    }
    float sortValue_{};
    unsigned short lightIndex_{};
    bool perVertex_{};
    bool negativeLight_{};
    bool operator < (const LitGeometryDescPacked& rhs) const
    {
        if (perVertex_ != rhs.perVertex_)
            return perVertex_ < rhs.perVertex_;
        return sortValue_ < rhs.sortValue_;
    }

};

struct ThreadedGroupedBatchQueue
{
    Vector<BatchQueue> threadQueues_;
};

class BatchCollector : public Object
{
    URHO3D_OBJECT(BatchCollector, Object);

public:
    BatchCollector(Context* context);
    void Initialize(bool threading, const PODVector<ScenePassInfo>& scenePasses);

    void Clear(unsigned frameNumber);
    /// Collect lights. Lights shall have the same indexes as used in LitGeometryDesc.
    void CollectLights(const PODVector<Light*>& lights);
    /// Collect visible geometries.
    void CollectVisibleGeometry(SceneGridDrawableSoA& sceneData);
    /// Collect lit geometries from given unsorted array of lit geometries.
    void CollectLitGeometries(const Vector<LitGeometryDescIdx>& litGeometries, SceneGridDrawableSoA& sceneData);

    const Vector<Drawable*>& GetVisibleGeometries() const { return visibleGeometries_; }
    const Vector<unsigned>& GetVisibleGeometriesNumLights() const { return numLightsPerVisibleGeometry_; }
    const Vector<LitGeometryDescPacked>& GetLitGeometries() const { return litGeometries_; }

private:
    WorkQueue* workQueue_{};
    bool threading_{};

    Vector<Light*> lights_;
    Vector<Drawable*> visibleGeometries_;
    Vector<unsigned> numLightsPerVisibleGeometry_;
    Vector<LitGeometryDescPacked> litGeometries_;

    Vector<UniquePtr<ThreadedGroupedBatchQueue>> scenePassQueues_;
};

}
