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

#include "../Core/WorkQueue.h"
#include "../Graphics/Drawable.h"
#include "../Graphics/StaticModel.h"
#include "../Graphics/Light.h"
#include "../Graphics/OcclusionBuffer.h"
#include "../Graphics/View.h"
#include "../Graphics/Camera.h"
#include "../Scene/Node.h"

namespace Urho3D
{

template <class T> static void ScheduleWork(WorkQueue* workQueue,
    unsigned threshold, unsigned numElements, unsigned bucketSize, const T& work)
{
    if (numElements > 0 && threshold != 0 && numElements <= threshold)
    {
        work(0, numElements, 0);
    }
    else
    {
        const unsigned numBuckets = Max(1u, numElements / bucketSize);
        const unsigned realBucketSize = (numElements - 1) / numBuckets + 1;
        for (unsigned i = 0; i < numBuckets; ++i)
        {
            const unsigned firstElementInBucket = i * realBucketSize;
            const unsigned remainingElements = numElements - firstElementInBucket;
            const unsigned numElementsInBucket = Min(realBucketSize, remainingElements);
            if (numElementsInBucket > 0)
            {
                workQueue->ScheduleWork([=](unsigned threadIndex)
                {
                    work(firstElementInBucket, firstElementInBucket + numElementsInBucket, threadIndex);
                });
            }
        }
    }
}

struct BasicDrawableProcessorSoA
{
    Vector<Drawable*> drawable_;

    Vector<bool> transformDirty_;
    Vector<Matrix3x4> transform_;
    Vector<BoundingBox> worldBoundingBox_;
    Vector<Sphere> worldBoundingSphere_;

    Vector<unsigned> drawableFlag_;
    Vector<unsigned> viewMask_;
    Vector<bool> occluder_;
    Vector<bool> occludee_;

    void Push(Drawable* drawable)
    {
        drawable_.Push(drawable);

        transformDirty_.Push(true);
        transform_.Push(Matrix3x4());
        worldBoundingBox_.Push(BoundingBox());
        worldBoundingSphere_.Push(Sphere());

        drawableFlag_.Push(drawable->GetDrawableFlags());
        viewMask_.Push(drawable->GetViewMask());
        occluder_.Push(drawable->IsOccluder());
        occludee_.Push(drawable->IsOccludee());
    }
    void EraseSwap(unsigned index)
    {
        drawable_.EraseSwap(index);
        transformDirty_.EraseSwap(index);
        transform_.EraseSwap(index);
        worldBoundingBox_.EraseSwap(index);
        worldBoundingSphere_.EraseSwap(index);
        drawableFlag_.EraseSwap(index);
        viewMask_.EraseSwap(index);
        occluder_.EraseSwap(index);
        occludee_.EraseSwap(index);
    }
    void Clear()
    {
        drawable_.Clear();
        transformDirty_.Clear();
        transform_.Clear();
        worldBoundingBox_.Clear();
        worldBoundingSphere_.Clear();
        drawableFlag_.Clear();
        viewMask_.Clear();
        occluder_.Clear();
        occludee_.Clear();
    }
    unsigned Size() const { return drawable_.Size(); }
};

class BasicDrawableProcessor
{
public:
    virtual ~BasicDrawableProcessor()
    {
        for (Drawable* drawable : data_.drawable_)
            drawable->SetDrawableIndex(DrawableIndex{});
        data_.Clear();
    }
    virtual unsigned GetSceneQueryThreshold() const { return 128; }
    const BasicDrawableProcessorSoA& GetData() const { return data_; }

    virtual void AddDrawable(Drawable* drawable)
    {
        assert(!drawable->GetDrawableIndex());

        DrawableIndex index;
        index.processor_ = this;
        index.index_ = data_.Size();

        data_.Push(drawable);

        drawable->SetDrawableIndex(index);
    }
    virtual void RemoveDrawable(Drawable* drawable)
    {
        const DrawableIndex drawableIndex = drawable->GetDrawableIndex();
        const unsigned index = drawableIndex.index_;

        assert(drawableIndex);
        assert(drawableIndex.processor_ == this);

        drawable->SetDrawableIndex(DrawableIndex{});

        data_.EraseSwap(index);
    }

    void MarkTransformDirty(unsigned index) { data_.transformDirty_[index] = true; }
    void SetViewMask(unsigned index, unsigned viewMask) { data_.viewMask_[index] = viewMask; }
    void SetOccluder(unsigned index, bool occluder) { data_.occluder_[index] = occluder; }
    void SetOccludee(unsigned index, bool occludee) { data_.occludee_[index] = occludee; }
    unsigned GetNumDrawables() const { return data_.Size(); }

    virtual void UpdateDirtyDrawables()
    {
        const unsigned numDrawables = GetNumDrawables();
        for (unsigned i = 0; i < numDrawables; ++i)
        {
            if (data_.transformDirty_[i])
            {
                Drawable* drawable = data_.drawable_[i];
                data_.transformDirty_[i] = false;
                data_.transform_[i] = drawable->GetNode()->GetWorldTransform();
                data_.worldBoundingBox_[i] = drawable->GetWorldBoundingBox();
                data_.worldBoundingSphere_[i] = Sphere(data_.worldBoundingBox_[i]);
            }
        }
    }

protected:
    bool CheckDrawableIsZone(unsigned index, unsigned viewMask) const
    {
        const unsigned flags = data_.drawableFlag_[index];
        return (data_.viewMask_[index] & viewMask) && (flags & DRAWABLE_ZONE);
    }
    bool CheckDrawableIsGeometry(unsigned index, unsigned viewMask) const
    {
        const unsigned flags = data_.drawableFlag_[index];
        return (data_.viewMask_[index] & viewMask) && (flags & DRAWABLE_GEOMETRY);
    }
    bool CheckDrawableIsLight(unsigned index, unsigned viewMask) const
    {
        const unsigned flags = data_.drawableFlag_[index];
        return (data_.viewMask_[index] & viewMask) && (flags & DRAWABLE_LIGHT);
    }
    bool CheckDrawableFlags(unsigned index, unsigned drawableFlags, unsigned viewMask) const
    {
        return (data_.drawableFlag_[index] & drawableFlags) && (data_.viewMask_[index] & viewMask);
    }
    bool CheckDrawableZoneOrOccluder(unsigned index, unsigned viewMask) const
    {
        const unsigned flags = data_.drawableFlag_[index];
        return (data_.viewMask_[index] & viewMask)
            && (flags == DRAWABLE_ZONE || (flags == DRAWABLE_GEOMETRY && data_.occluder_[index]));
    }
    bool CheckDrawableInFrustum(unsigned index, const Frustum& frustum) const
    {
        const Sphere& boundingSpere = data_.worldBoundingSphere_[index];
        const BoundingBox& boundingBox = data_.worldBoundingBox_[index];

        const Intersection sphereIntersection = frustum.IsInside(boundingSpere);
        return sphereIntersection != OUTSIDE
            && (sphereIntersection == INSIDE || frustum.IsInsideFast(boundingBox));
    }
public:
    BasicDrawableProcessorSoA data_;
};

struct QueryRangeOwner
{
    QueryRangeOwner() = default;
    QueryRangeOwner(BasicDrawableProcessor* processor, unsigned beginIndex, unsigned endIndex)
        : processor_(processor), beginIndex_(beginIndex), endIndex_(endIndex) {}
    BasicDrawableProcessor* processor_{};
    unsigned beginIndex_{};
    unsigned endIndex_{};
};

inline void MergeRangeOwnerVectors(Vector<QueryRangeOwner>& dest, Vector<QueryRangeOwner>& src, unsigned baseIndex)
{
    for (QueryRangeOwner rangeOwner : src)
    {
        rangeOwner.beginIndex_ += baseIndex;
        rangeOwner.endIndex_ += baseIndex;
    }
}

struct ZonesInFrustumQueryResult
{
    Vector<Zone*> zones_;
    void Clear()
    {
        zones_.Clear();
    }
    void Merge(Vector<ZonesInFrustumQueryResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (ZonesInFrustumQueryResult& result : array)
            {
                zones_.Push(result.zones_);
            }
        }
        else
        {
            Swap(zones_, array[0].zones_);
        }
    }
};

class ZoneProcessor : public BasicDrawableProcessor
{
public:
    virtual void QueryZonesInFrustum(unsigned viewMask, const Frustum& frustum,
        WorkQueue* workQueue, unsigned bucketSize, Vector<ZonesInFrustumQueryResult>& threadResults)
    {
        ScheduleWork(workQueue, GetSceneQueryThreshold(), GetNumDrawables(), bucketSize,
            [this, viewMask, frustum, &threadResults](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                if (CheckDrawableIsZone(index, viewMask) && CheckDrawableInFrustum(index, frustum))
                {
                    threadResults[threadIndex].zones_.Push(static_cast<Zone*>(data_.drawable_[index]));
                }
            }
        });
    }
};

struct LightsInFrustumQueryResult
{
    Vector<QueryRangeOwner> owners_;
    Vector<Light*> lights_;
    void Clear()
    {
        owners_.Clear();
        lights_.Clear();
    }
    void Merge(Vector<LightsInFrustumQueryResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (LightsInFrustumQueryResult& result : array)
            {
                MergeRangeOwnerVectors(owners_, result.owners_, lights_.Size());
                lights_.Push(result.lights_);
            }
        }
        else
        {
            Swap(owners_, array[0].owners_);
            Swap(lights_, array[0].lights_);
        }
    }
};

class LightProcessor : public BasicDrawableProcessor
{
public:
    virtual void QueryLightsInFrustum(unsigned viewMask, const Frustum& frustum,
        WorkQueue* workQueue, unsigned bucketSize, Vector<LightsInFrustumQueryResult>& threadResults)
    {
        ScheduleWork(workQueue, GetSceneQueryThreshold(), GetNumDrawables(), bucketSize,
            [this, viewMask, frustum, &threadResults](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            threadResults[threadIndex].owners_.Push(QueryRangeOwner{ this, beginIndex, endIndex });
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                if (CheckDrawableIsLight(index, viewMask) && CheckDrawableInFrustum(index, frustum))
                {
                    threadResults[threadIndex].lights_.Push(static_cast<Light*>(data_.drawable_[index]));
                }
            }
        });
    }
    virtual void QueryLightsOcclusion(OcclusionBuffer* buffer, LightsInFrustumQueryResult& lightsArray,
        WorkQueue* workQueue, unsigned beginIndex, unsigned endIndex, Vector<LightsInFrustumQueryResult>& threadResults) const
    {
        ScheduleWork(workQueue, 0, lightsArray.lights_.Size(), 0,
            [](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
        });
    }

};

struct GeometryOccludersInFrustumQueryResult
{
    Vector<Drawable*> occluders_;
    void Clear()
    {
        occluders_.Clear();
    }
    void Merge(Vector<GeometryOccludersInFrustumQueryResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (GeometryOccludersInFrustumQueryResult& result : array)
            {
                occluders_.Push(result.occluders_);
            }
        }
        else
        {
            Swap(occluders_, array[0].occluders_);
        }
    }
};

struct GeometriesInFrustumQueryResult
{
    Vector<QueryRangeOwner> owners_;
    Vector<Drawable*> geometries_;
    void Clear()
    {
        owners_.Clear();
        geometries_.Clear();
    }
    void Merge(Vector<GeometriesInFrustumQueryResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (GeometriesInFrustumQueryResult& result : array)
            {
                MergeRangeOwnerVectors(owners_, result.owners_, geometries_.Size());
                geometries_.Push(result.geometries_);
            }
        }
        else
        {
            Swap(owners_, array[0].owners_);
            Swap(geometries_, array[0].geometries_);
        }
    }
};

class GeometryProcessor : public BasicDrawableProcessor
{
public:
    virtual void QueryGeometryOccludersInFrustum(unsigned viewMask, const Frustum& frustum,
        WorkQueue* workQueue, unsigned bucketSize, Vector<GeometryOccludersInFrustumQueryResult>& threadResults)
    {
        ScheduleWork(workQueue, GetSceneQueryThreshold(), GetNumDrawables(), bucketSize,
            [this, viewMask, frustum, &threadResults](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                if (CheckDrawableIsGeometry(index, viewMask) && data_.occluder_[index] && CheckDrawableInFrustum(index, frustum))
                {
                    threadResults[threadIndex].occluders_.Push(data_.drawable_[index]);
                }
            }
        });
    }
    virtual void QueryGeometriesInFrustum(unsigned viewMask, const Frustum& frustum,
        WorkQueue* workQueue, unsigned bucketSize, Vector<GeometriesInFrustumQueryResult>& threadResults)
    {
        ScheduleWork(workQueue, GetSceneQueryThreshold(), GetNumDrawables(), bucketSize,
            [this, viewMask, frustum, &threadResults](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            threadResults[threadIndex].owners_.Push(QueryRangeOwner{ this, beginIndex, endIndex });
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                if (CheckDrawableIsGeometry(index, viewMask) && CheckDrawableInFrustum(index, frustum))
                {
                    threadResults[threadIndex].geometries_.Push(data_.drawable_[index]);
                }
            }
        });
    }
};

struct SceneFrustumQueryResult
{
    Vector<Drawable*> drawables_;
};

class StaticModelProcessor : public GeometryProcessor
{
public:
    virtual void QueryZonesInFrustum(unsigned viewMask, const Frustum& frustum,
        WorkQueue* workQueue, unsigned bucketSize, Vector<ZonesInFrustumQueryResult>& threadResults) const
    {
        // Nope
    }
    virtual void QueryLightsInFrustum(unsigned viewMask, const Frustum& frustum,
        WorkQueue* workQueue, unsigned bucketSize, Vector<LightsInFrustumQueryResult>& threadResults) const
    {
        // Nope
    }

};

class DrawableProcessor
{
public:
    DrawableProcessor()
    {
        AddProcessor<GeometryProcessor>(StringHash());
        AddProcessor<ZoneProcessor, Zone>();
        AddProcessor<LightProcessor, Light>();
        AddProcessor<StaticModelProcessor, StaticModel>();
    }

    template <class T> void AddProcessor(StringHash type)
    {
        assert(!specializedProcessors_.Contains(type));
        UniquePtr<BasicDrawableProcessor> processor(new T());
        if (type)
            specializedProcessors_[type] = processor.Get();
        else
            defaultProcessor_ = processor.Get();

        if (auto zoneProcessor = dynamic_cast<ZoneProcessor*>(processor.Get()))
            zoneProcessors_.Push(zoneProcessor);
        if (auto lightProcessor = dynamic_cast<LightProcessor*>(processor.Get()))
            lightProcessors_.Push(lightProcessor);
        if (auto geometryProcessor = dynamic_cast<GeometryProcessor*>(processor.Get()))
            geometryProcessors_.Push(geometryProcessor);
        processors_.Push(std::move(processor));
    }

    template <class T, class U> void AddProcessor() { AddProcessor<T>(U::GetTypeStatic()); }

    BasicDrawableProcessor& GetProcessor(const StringHash& type)
    {
        auto iter = specializedProcessors_.Find(type);
        return iter == specializedProcessors_.End() ? *defaultProcessor_ : *iter->second_;
    }

    void AddDrawable(Drawable* drawable)
    {
        const StringHash type = drawable->GetType();
        BasicDrawableProcessor& processor = GetProcessor(type);
        processor.AddDrawable(drawable);
    }

    void RemoveDrawable(Drawable* drawable)
    {
        const StringHash type = drawable->GetType();
        BasicDrawableProcessor& processor = GetProcessor(type);
        processor.RemoveDrawable(drawable);
    }

    /// Prepare for multi-threading. Calculate total amount of drawables and optimal chunk sizes.
    void PrepareForThreading(WorkQueue* workQueue)
    {
        workQueue_ = workQueue;
        numThreads_ = workQueue_->GetNumThreads() + 1;

        // Calculate numbers
        numDrawables_ = 0;
        maxNumDrawablesInProcessor_ = 0;
        for (const auto& processor : processors_)
        {
            const unsigned count = processor->GetNumDrawables();
            numDrawables_ += count;
            maxNumDrawablesInProcessor_ = Max(maxNumDrawablesInProcessor_, count);
        }

        // Compute optimal chunk size
        numDrawablesInBucket_ = maxNumDrawablesInProcessor_ / numThreads_ + 1;
    }
    void UpdateDirtyDrawables()
    {
        // This may cause Node transforms recalculation, so it mustn't be threaded
        for (const auto& processor : processors_)
            processor->UpdateDirtyDrawables();
    }

    template <class T> void ClearThreadCache(T& cache)
    {
        cache.Resize(numThreads_);
        for (auto& result : cache)
        {
            result.Clear();
        }
    }
    template <class T> void MergeThreadCache(T& result, Vector<T>& cache)
    {
        workQueue_->Complete(M_MAX_UNSIGNED);
        result.Merge(cache);
    }

    void QueryZonesInFrustum(unsigned viewMask, const Frustum& frustum)
    {
        ClearThreadCache(zonesInFrustumQueryCache_);
        for (const auto& processor : zoneProcessors_)
            processor->QueryZonesInFrustum(viewMask, frustum, workQueue_, numDrawablesInBucket_, zonesInFrustumQueryCache_);
        MergeThreadCache(zonesInFrustumQuery_, zonesInFrustumQueryCache_);
    }
    const ZonesInFrustumQueryResult& GetZonesInFrustum() const { return zonesInFrustumQuery_; }

    void QueryGeometryOccludersInFrustum(unsigned viewMask, const Frustum& frustum)
    {
        ClearThreadCache(geometryOccludersInFrustumQueryCache_);
        for (const auto& processor : geometryProcessors_)
            processor->QueryGeometryOccludersInFrustum(viewMask, frustum, workQueue_, numDrawablesInBucket_, geometryOccludersInFrustumQueryCache_);
        MergeThreadCache(geometryOccludersInFrustumQuery_, geometryOccludersInFrustumQueryCache_);
    }
    const GeometryOccludersInFrustumQueryResult& GetGeometryOccludersInFrustum() const { return geometryOccludersInFrustumQuery_; }

    void QueryLightsInFrustum(unsigned viewMask, const Frustum& frustum)
    {
        ClearThreadCache(lightsInFrustumQueryCache_);
        for (const auto& processor : lightProcessors_)
            processor->QueryLightsInFrustum(viewMask, frustum, workQueue_, numDrawablesInBucket_, lightsInFrustumQueryCache_);
        MergeThreadCache(lightsInFrustumQuery_, lightsInFrustumQueryCache_);
    }
    const LightsInFrustumQueryResult& GetLightsInFrustum() const { return lightsInFrustumQuery_; }

    void QueryGeometriesInFrustum(unsigned viewMask, const Frustum& frustum)
    {
        ClearThreadCache(geometriesInFrustumQueryCache_);
        for (const auto& processor : geometryProcessors_)
            processor->QueryGeometriesInFrustum(viewMask, frustum, workQueue_, numDrawablesInBucket_, geometriesInFrustumQueryCache_);
        MergeThreadCache(geometriesInFrustumQuery_, geometriesInFrustumQueryCache_);
    }
    const GeometriesInFrustumQueryResult& GetGeometriesInFrustum() const { return geometriesInFrustumQuery_; }

    void QueryLightsOcclusion(OcclusionBuffer* buffer)
    {
        ClearThreadCache(lightsInFrustumQueryCache_);

        for (const QueryRangeOwner& owner : lightsInFrustumQuery_.owners_)
        {
//             owner.processor_->QueryLightsOcclusion(buffer, lightsInFrustumQuery_,
//                 owner.beginIndex_, owner.endIndex_, lightsInFrustumQueryCache_);
        }

        MergeThreadCache(lightsInFrustumQuery_, lightsInFrustumQueryCache_);
    }

    void ProcessLights()
    {

    }

private:
    Vector<UniquePtr<BasicDrawableProcessor>> processors_;
    Vector<ZoneProcessor*> zoneProcessors_;
    Vector<LightProcessor*> lightProcessors_;
    Vector<GeometryProcessor*> geometryProcessors_;
    HashMap<StringHash, BasicDrawableProcessor*> specializedProcessors_;
    BasicDrawableProcessor* defaultProcessor_{};

    WorkQueue* workQueue_{};
    unsigned numThreads_{};
    unsigned numDrawables_{};
    unsigned maxNumDrawablesInProcessor_{};
    unsigned numDrawablesInBucket_{};

    Vector<ZonesInFrustumQueryResult> zonesInFrustumQueryCache_;
    ZonesInFrustumQueryResult zonesInFrustumQuery_;

    Vector<GeometryOccludersInFrustumQueryResult> geometryOccludersInFrustumQueryCache_;
    GeometryOccludersInFrustumQueryResult geometryOccludersInFrustumQuery_;

    Vector<LightsInFrustumQueryResult> lightsInFrustumQueryCache_;
    LightsInFrustumQueryResult lightsInFrustumQuery_;

    Vector<GeometriesInFrustumQueryResult> geometriesInFrustumQueryCache_;
    GeometriesInFrustumQueryResult geometriesInFrustumQuery_;
};

}
