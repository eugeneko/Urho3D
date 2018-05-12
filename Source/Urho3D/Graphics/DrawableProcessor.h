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
#include "../Graphics/Renderer.h"
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

struct SceneProcessorDrawableSoA
{
    unsigned size_ = 0;

    Vector<Drawable*> drawable_;
    Vector<StringHash> drawableType_;
    Vector<unsigned> drawableFlag_;

    Vector<bool> dirty_;
    Vector<Matrix3x4> transform_;
    Vector<BoundingBox> worldBoundingBox_;
    Vector<Sphere> worldBoundingSphere_;

    Vector<float> drawDistance_;
    Vector<unsigned> viewMask_;
    Vector<unsigned> lightMask_;
    Vector<unsigned> shadowMask_;
    Vector<unsigned> zoneMask_;
    Vector<bool> occluder_;
    Vector<bool> occludee_;

    Vector<Zone*> cachedZone_;
    Vector<bool> cachedZoneDirty_;
    Vector<unsigned> cachedZoneViewMask_;
    Vector<unsigned> cachedZoneLightMask_;
    Vector<unsigned> cachedZoneShadowMask_;

    bool IsValid() const
    {
        return size_ == drawable_.Size()
            && size_ == drawableType_.Size()
            && size_ == drawableFlag_.Size()

            && size_ == dirty_.Size()
            && size_ == transform_.Size()
            && size_ == worldBoundingBox_.Size()
            && size_ == worldBoundingSphere_.Size()

            && size_ == drawDistance_.Size()
            && size_ == viewMask_.Size()
            && size_ == lightMask_.Size()
            && size_ == shadowMask_.Size()
            && size_ == zoneMask_.Size()
            && size_ == occluder_.Size()
            && size_ == occludee_.Size()

            && size_ == cachedZone_.Size()
            && size_ == cachedZoneDirty_.Size()
            && size_ == cachedZoneViewMask_.Size()
            && size_ == cachedZoneLightMask_.Size()
            && size_ == cachedZoneShadowMask_.Size();
    }
    void Push(Drawable* drawable)
    {
        ++size_;

        drawable_.Push(drawable);
        drawableType_.Push(drawable->GetType());
        drawableFlag_.Push(drawable->GetDrawableFlags());

        dirty_.Push(true);
        transform_.Push(Matrix3x4());
        worldBoundingBox_.Push(BoundingBox());
        worldBoundingSphere_.Push(Sphere());

        drawDistance_.Push(drawable->GetDrawDistance()); // TODO(eugeneko) Add update
        viewMask_.Push(drawable->GetViewMask());
        lightMask_.Push(drawable->GetLightMask()); // TODO(eugeneko) Add update
        shadowMask_.Push(drawable->GetLightMask()); // TODO(eugeneko) Add update
        zoneMask_.Push(drawable->GetZoneMask()); // TODO(eugeneko) Add update
        occluder_.Push(drawable->IsOccluder());
        occludee_.Push(drawable->IsOccludee());

        cachedZone_.Push(nullptr);
        cachedZoneDirty_.Push(true);
        cachedZoneViewMask_.Push(0); // TODO(eugeneko) Add update
        cachedZoneLightMask_.Push(0); // TODO(eugeneko) Add update
        cachedZoneShadowMask_.Push(0); // TODO(eugeneko) Add update

        assert(IsValid());
    }
    void EraseSwap(unsigned index)
    {
        assert(index < size_);
        --size_;

        drawable_.EraseSwap(index);
        drawableType_.EraseSwap(index);
        drawableFlag_.EraseSwap(index);

        dirty_.EraseSwap(index);
        transform_.EraseSwap(index);
        worldBoundingBox_.EraseSwap(index);
        worldBoundingSphere_.EraseSwap(index);

        drawDistance_.EraseSwap(index);
        viewMask_.EraseSwap(index);
        lightMask_.EraseSwap(index);
        shadowMask_.EraseSwap(index);
        zoneMask_.EraseSwap(index);
        occluder_.EraseSwap(index);
        occludee_.EraseSwap(index);

        cachedZone_.EraseSwap(index);
        cachedZoneDirty_.EraseSwap(index);
        cachedZoneViewMask_.EraseSwap(index);
        cachedZoneLightMask_.EraseSwap(index);
        cachedZoneShadowMask_.EraseSwap(index);

        assert(IsValid());
    }
    void Clear()
    {
        size_ = 0;

        drawable_.Clear();
        drawableType_.Clear();
        drawableFlag_.Clear();

        dirty_.Clear();
        transform_.Clear();
        worldBoundingBox_.Clear();
        worldBoundingSphere_.Clear();

        drawDistance_.Clear();
        viewMask_.Clear();
        lightMask_.Clear();
        shadowMask_.Clear();
        zoneMask_.Clear();
        occluder_.Clear();
        occludee_.Clear();

        cachedZone_.Clear();
        cachedZoneDirty_.Clear();
        cachedZoneViewMask_.Clear();
        cachedZoneLightMask_.Clear();
        cachedZoneShadowMask_.Clear();

        assert(IsValid());
    }

    bool MatchViewMask(unsigned index, unsigned viewMask) const { return !!(viewMask_[index] & viewMask); }
    bool IsZone(unsigned index) const { return !!(drawableFlag_[index] & DRAWABLE_ZONE); }
    bool IsLight(unsigned index) const { return !!(drawableFlag_[index] & DRAWABLE_LIGHT); }
    bool IsGeometry(unsigned index) const { return !!(drawableFlag_[index] & DRAWABLE_GEOMETRY); }
    Intersection IntersectSphereFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInside(worldBoundingSphere_[index]); }
    Intersection IntersectBoxFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInsideFast(worldBoundingBox_[index]); }
    bool IsInFrustum(unsigned index, const Frustum& frustum) const
    {
        const Intersection sphereIntersection = IntersectSphereFrustum(index, frustum);
        return sphereIntersection != OUTSIDE
            && (sphereIntersection == INSIDE || IntersectBoxFrustum(index, frustum) != OUTSIDE);
    }
    bool IsCenterInFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInside(worldBoundingSphere_[index].center_) != OUTSIDE; }

};

/// Contains persistent scene state.
class SceneProcessor
{
public:
    virtual ~SceneProcessor()
    {
        for (Drawable* drawable : data_.drawable_)
            drawable->SetDrawableIndex(DrawableIndex{});
        data_.Clear();
    }

    void MarkDirty(unsigned index) { data_.dirty_[index] = true; }
    void SetViewMask(unsigned index, unsigned viewMask) { data_.viewMask_[index] = viewMask; }
    void SetOccluder(unsigned index, bool occluder) { data_.occluder_[index] = occluder; }
    void SetOccludee(unsigned index, bool occludee) { data_.occludee_[index] = occludee; }
    void SetCachedZone(unsigned index, Zone* zone) { data_.cachedZone_[index] = zone; }

    virtual void AddDrawable(Drawable* drawable)
    {
        assert(!drawable->GetDrawableIndex());

        DrawableIndex index;
        index.processor_ = this;
        index.index_ = data_.size_;

        data_.Push(drawable);

        drawable->SetDrawableIndex(index);
    }
    virtual void RemoveDrawable(Drawable* drawable)
    {
        const DrawableIndex drawableIndex = drawable->GetDrawableIndex();
        const unsigned index = drawableIndex.index_;

        assert(drawableIndex);
        assert(drawableIndex.processor_ == this);

        if (data_.IsZone(index))
        {
            // Remove cached zone
            for (unsigned i = 0; i < data_.size_; ++i)
            {
                if (data_.cachedZone_[i] == drawable)
                {
                    data_.cachedZone_[i] = nullptr;
                    data_.cachedZoneDirty_[i] = true;
                }
            }
        }
        drawable->SetDrawableIndex(DrawableIndex{});

        data_.EraseSwap(index);
    }
    virtual void MarkZonesDirty()
    {
        dirtyZoneCache_ = true;
    }

    /// Update dirty values. Call it before using persistent data.
    virtual void UpdateDirtyDrawables()
    {
        bool allZonesDirty = false;
        for (unsigned i = 0; i < data_.size_; ++i)
        {
            if (data_.dirty_[i])
            {
                Drawable* drawable = data_.drawable_[i];
                data_.dirty_[i] = false;
                data_.transform_[i] = drawable->GetNode()->GetWorldTransform();
                data_.worldBoundingBox_[i] = drawable->GetWorldBoundingBox();
                data_.worldBoundingSphere_[i] = Sphere(data_.worldBoundingBox_[i]);

                // Mark zone dirty
                data_.cachedZoneDirty_[i] = true;

                // Mark dirty zones
                if (data_.IsZone(i))
                    allZonesDirty = true;
            }
        }

        // Reset cached zones if any is dirty
        if (allZonesDirty)
        {
            for (unsigned i = 0; i < data_.size_; ++i)
                data_.cachedZoneDirty_[i] = true;
        }
    }

    SceneProcessorDrawableSoA& GetData() { return data_; }

private:
    SceneProcessorDrawableSoA data_;
    bool dirtyZoneCache_ = false;
};

struct SceneQuerySpiritsResult
{
    Vector<Light*> lights_;
    Vector<Drawable*> occluders_;
    Vector<Zone*> zones_;
    void Clear()
    {
        lights_.Clear();
        occluders_.Clear();
        zones_.Clear();
    }
    void Merge(Vector<SceneQuerySpiritsResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (SceneQuerySpiritsResult& result : array)
            {
                lights_.Push(result.lights_);
                occluders_.Push(result.occluders_);
                zones_.Push(result.zones_);
            }
        }
        else
        {
            Swap(lights_, array[0].lights_);
            Swap(occluders_, array[0].occluders_);
            Swap(zones_, array[0].zones_);
        }
    }
};

struct SceneQueryGeometriesResult
{
    Vector<Drawable*> geometries_;
    Vector<Zone*> zones_;
    Vector<unsigned> zoneLightMasks_;
    Vector<unsigned> zoneShadowMasks_;
    Vector<float> distances_;
    Vector<Vector2> minmaxZ_;

    float minZ_{};
    float maxZ_{};
    bool IsValid() const
    {
        return geometries_.Size() == zones_.Size()
            && geometries_.Size() == zoneLightMasks_.Size()
            && geometries_.Size() == zoneShadowMasks_.Size()
            && geometries_.Size() == distances_.Size()
            && geometries_.Size() == minmaxZ_.Size();
    }
    void Clear()
    {
        geometries_.Clear();
        zones_.Clear();
        zoneLightMasks_.Clear();
        zoneShadowMasks_.Clear();
        distances_.Clear();
        minmaxZ_.Clear();
        minZ_ = M_INFINITY;
        maxZ_ = 0.0f;

        assert(IsValid());
    }
    void Merge(Vector<SceneQueryGeometriesResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (SceneQueryGeometriesResult& result : array)
            {
                geometries_.Push(result.geometries_);
                zones_.Push(result.zones_);
                zoneLightMasks_.Push(result.zoneLightMasks_);
                zoneShadowMasks_.Push(result.zoneShadowMasks_);
                distances_.Push(result.distances_);
                minmaxZ_.Push(result.minmaxZ_);
                minZ_ = Min(minZ_, result.minZ_);
                maxZ_ = Max(maxZ_, result.maxZ_);
            }
        }
        else
        {
            Swap(geometries_, array[0].geometries_);
            Swap(zones_, array[0].zones_);
            Swap(zoneLightMasks_, array[0].zoneLightMasks_);
            Swap(zoneShadowMasks_, array[0].zoneShadowMasks_);
            Swap(distances_, array[0].distances_);
            Swap(minmaxZ_, array[0].minmaxZ_);
            minZ_ = array[0].minZ_;
            maxZ_ = array[0].maxZ_;
        }
        if (minZ_ == M_INFINITY)
            minZ_ = 0.0f;

        assert(IsValid());
    }
};

struct ViewCookedZonesData
{
    bool cameraZoneOverride_{};
    Vector<Zone*> zones_;
    Zone* cameraZone_{};
    Zone* farClipZone_{};
};

/// Contains view-specific scene data.
class ViewProcessor
{
public:
    virtual unsigned GetSceneQueryThreshold() { return 64; }
    /// Query zones, lights and occluders.
    virtual void QuerySceneSpirits(WorkQueue* workQueue, const SceneProcessorDrawableSoA& sceneData,
        unsigned viewMask, const Frustum& frustum)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        sceneSpiritsQueryThreadResults_.Resize(numThreads);
        for (SceneQuerySpiritsResult& result : sceneSpiritsQueryThreadResults_)
            result.Clear();

        workQueue->ScheduleWork(GetSceneQueryThreshold(), sceneData.size_, numThreads,
            [this, viewMask, frustum, &sceneData](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            SceneQuerySpiritsResult& result = sceneSpiritsQueryThreadResults_[threadIndex];
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                if (!sceneData.MatchViewMask(index, viewMask))
                    continue;

                Drawable* drawable = sceneData.drawable_[index];

                // TODO(eugeneko) Get rid of branching
                if (sceneData.IsZone(index))
                {
                    result.zones_.Push(static_cast<Zone*>(drawable));
                }
                else if (sceneData.IsLight(index))
                {
                    if (sceneData.IsInFrustum(index, frustum))
                        result.lights_.Push(static_cast<Light*>(drawable));
                }
                else if (sceneData.IsGeometry(index) && sceneData.occluder_[index])
                {
                    if (sceneData.IsInFrustum(index, frustum))
                        result.occluders_.Push(drawable);
                }
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
        sceneSpiritsQueryResult_.Merge(sceneSpiritsQueryThreadResults_);
    }
    /// Prepare zones for further usage.
    virtual void CookZones(Camera* cullCamera, Zone* defaultZone)
    {
        zones_.zones_.Clear();
        zones_.zones_.Push(sceneSpiritsQueryResult_.zones_);
        zones_.cameraZoneOverride_ = false;
        zones_.cameraZone_ = nullptr;
        zones_.farClipZone_ = nullptr;

        // Sort zones
        Sort(zones_.zones_.Begin(), zones_.zones_.End(),
            [](Zone* lhs, Zone* rhs) { return lhs->GetPriority() > rhs->GetPriority(); });

        // Find camera and far clip zones
        Node* cullCameraNode = cullCamera->GetNode();
        const Vector3 cameraPos = cullCameraNode->GetWorldPosition();
        const Vector3 farClipPos = cameraPos + cullCameraNode->GetWorldDirection() * Vector3(0.0f, 0.0f, cullCamera->GetFarClip());
        for (Zone* zone : zones_.zones_)
        {
            if (!zones_.cameraZone_ && zone->IsInside(cameraPos))
                zones_.cameraZone_ = zone;
            if (!zones_.farClipZone_ && zone->IsInside(farClipPos))
                zones_.farClipZone_ = zone;
            if (zones_.cameraZone_ && zones_.farClipZone_)
                break;
        }

        if (!zones_.cameraZone_)
            zones_.cameraZone_ = defaultZone;

        if (!zones_.farClipZone_)
            zones_.farClipZone_ = defaultZone;
    }
    /// Query geometries.
    virtual void QuerySceneGeometries(WorkQueue* workQueue, SceneProcessorDrawableSoA& sceneData,
        Camera* cullCamera, OcclusionBuffer* occlusionBuffer)
    {
        const unsigned viewMask = cullCamera->GetViewMask();
        const Frustum& frustum = cullCamera->GetFrustum();
        const Matrix3x4& viewMatrix = cullCamera->GetView();
        const Vector3 viewZ = Vector3(viewMatrix.m20_, viewMatrix.m21_, viewMatrix.m22_);
        const Vector3 absViewZ = viewZ.Abs();

        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        sceneGeometriesQueryThreadResults_.Resize(numThreads);
        for (SceneQueryGeometriesResult& result : sceneGeometriesQueryThreadResults_)
            result.Clear();

        workQueue->ScheduleWork(GetSceneQueryThreshold(), sceneData.size_, numThreads,
            [=, &sceneData](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            SceneQueryGeometriesResult& result = sceneGeometriesQueryThreadResults_[threadIndex];
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                // Discard drawables from other views
                if (!sceneData.MatchViewMask(index, viewMask))
                    continue;

                // Discard spirits and drawables outside the frustum
                if (!sceneData.IsGeometry(index) || !sceneData.IsInFrustum(index, frustum))
                    continue;

                // Calculate drawable distance
                const float drawDistance = sceneData.drawDistance_[index];
                const float distance = cullCamera->GetDistance(sceneData.worldBoundingSphere_[index].center_);

                // Discard if drawable is too far
                if (drawDistance > 0.0f && distance > drawDistance)
                    continue;

                // Discard using occlusion buffer
                const BoundingBox& boundingBox = sceneData.worldBoundingBox_[index];
                if (occlusionBuffer && sceneData.occludee_[index] && !occlusionBuffer->IsVisible(boundingBox))
                    continue;

                // Update drawable zone
                if (!zones_.cameraZoneOverride_ && !zones_.zones_.Empty())
                {
                    UpdateDirtyZone(sceneData, index, viewMask, frustum);
                }

                // Update min and max Z
                float minZ = M_LARGE_VALUE;
                float maxZ = M_LARGE_VALUE;

                const Vector3 center = boundingBox.Center();
                const Vector3 edge = boundingBox.Size() * 0.5f;

                // Do not add "infinite" objects like skybox to prevent shadow map focusing behaving erroneously
                if (edge.LengthSquared() < M_LARGE_VALUE * M_LARGE_VALUE)
                {
                    const float viewCenterZ = viewZ.DotProduct(center) + viewMatrix.m23_;
                    const float viewEdgeZ = absViewZ.DotProduct(edge);
                    minZ = viewCenterZ - viewEdgeZ;
                    maxZ = viewCenterZ + viewEdgeZ;
                    result.minZ_ = Min(result.minZ_, minZ);
                    result.maxZ_ = Max(result.maxZ_, maxZ);
                }

                // Push results
                result.geometries_.Push(sceneData.drawable_[index]);
                result.zones_.Push(sceneData.cachedZone_[index]);
                result.zoneLightMasks_.Push(sceneData.cachedZoneLightMask_[index]);
                result.zoneShadowMasks_.Push(sceneData.cachedZoneShadowMask_[index]);
                result.distances_.Push(distance);
                result.minmaxZ_.Push(Vector2(minZ, maxZ));
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
        sceneGeometriesQueryResult_.Merge(sceneGeometriesQueryThreadResults_);
    }
    /// Update geometries.
    virtual void UpdateGeometries(WorkQueue* workQueue, const FrameInfo& frame)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        workQueue->ScheduleWork(GetSceneQueryThreshold(), sceneGeometriesQueryResult_.geometries_.Size(), numThreads,
            [this, &frame](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                // TODO(eugeneko) REMOVE!!
                Drawable* drawable = sceneGeometriesQueryResult_.geometries_[index];
                drawable->UpdateBatches(frame);
                drawable->SetZone(sceneGeometriesQueryResult_.zones_[index]);
                drawable->MarkInView(frame);
                drawable->SetMinMaxZ(sceneGeometriesQueryResult_.minmaxZ_[index].x_, sceneGeometriesQueryResult_.minmaxZ_[index].y_);
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
    }

    const SceneQuerySpiritsResult& GetSpirits() const { return sceneSpiritsQueryResult_; }
    const SceneQueryGeometriesResult& GetGeometries() const { return sceneGeometriesQueryResult_; }

protected:
    void UpdateDirtyZone(SceneProcessorDrawableSoA& sceneData, unsigned index, unsigned viewMask, const Frustum& frustum) const
    {
        const Vector3 drawableCenter = sceneData.worldBoundingSphere_[index].center_;

        Zone* cachedZone = sceneData.cachedZone_[index];
        const bool cachedZoneDirty = sceneData.cachedZoneDirty_[index];
        const unsigned cachedZoneViewMask = sceneData.cachedZoneViewMask_[index];

        // TODO(eugeneko) Is branch is not optimized for multiple zones
        if (cachedZoneDirty || !cachedZone || !(cachedZoneViewMask & viewMask))
        {
            // Find new zone
            const bool temporary = sceneData.IsCenterInFrustum(index, frustum);
            const int highestZonePriority = zones_.zones_.Empty() ? M_MIN_INT : zones_.zones_[0]->GetPriority();

            Zone* newZone = nullptr;

            // First check if the current zone remains a conclusive result
            if (cachedZone && (cachedZone->GetViewMask() & viewMask)
                && cachedZone->GetPriority() >= highestZonePriority
                && (sceneData.zoneMask_[index] & cachedZone->GetZoneMask())
                && cachedZone->IsInside(drawableCenter))
            {
                newZone = cachedZone;
            }
            else
            {
                // Search for appropriate zone
                for (Zone* zone : zones_.zones_)
                {
                    if ((sceneData.zoneMask_[index] & zone->GetZoneMask()) && zone->IsInside(drawableCenter))
                    {
                        newZone = zone;
                        break;
                    }
                }
            }

            // Setup new zone
            sceneData.cachedZoneDirty_[index] = !temporary;
            sceneData.cachedZone_[index] = newZone;
            sceneData.cachedZoneLightMask_[index] = newZone->GetLightMask();
            sceneData.cachedZoneShadowMask_[index] = newZone->GetShadowMask();
        }
    }

private:
    Vector<SceneQuerySpiritsResult> sceneSpiritsQueryThreadResults_;
    SceneQuerySpiritsResult sceneSpiritsQueryResult_;

    ViewCookedZonesData zones_;

    Vector<SceneQueryGeometriesResult> sceneGeometriesQueryThreadResults_;
    SceneQueryGeometriesResult sceneGeometriesQueryResult_;


};

struct BasicDrawableProcessorSoA
{
    Vector<Drawable*> drawable_;

    Vector<bool> transformDirty_;
    Vector<Matrix3x4> transform_;
    Vector<BoundingBox> worldBoundingBox_;
    Vector<Sphere> worldBoundingSphere_;

    Vector<unsigned> drawableFlag_;
    Vector<unsigned> viewMask_;
    Vector<unsigned> lightMask_;
    Vector<unsigned> zoneMask_;
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
        lightMask_.Push(drawable->GetLightMask()); // TODO(eugeneko) Add update
        zoneMask_.Push(drawable->GetZoneMask()); // TODO(eugeneko) Add update
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
        lightMask_.EraseSwap(index);
        zoneMask_.EraseSwap(index);
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
        lightMask_.Clear();
        zoneMask_.Clear();
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
//         assert(!drawable->GetDrawableIndex());
//
//         DrawableIndex index;
//         index.processor_ = this;
//         index.index_ = data_.Size();
//
//         data_.Push(drawable);
//
//         drawable->SetDrawableIndex(index);
    }
    virtual void RemoveDrawable(Drawable* drawable)
    {
//         const DrawableIndex drawableIndex = drawable->GetDrawableIndex();
//         const unsigned index = drawableIndex.index_;
//
//         assert(drawableIndex);
//         assert(drawableIndex.processor_ == this);
//
//         drawable->SetDrawableIndex(DrawableIndex{});
//
//         data_.EraseSwap(index);
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
        dest.Push(rangeOwner);
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

class ZoneDrawableProcessor : public BasicDrawableProcessor
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

class LightDrawableProcessor : public BasicDrawableProcessor
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
    Vector<Sphere> worldBoundingSpheres_;
    Vector<unsigned> lightMasks_;
    Vector<unsigned> zoneMasks_;
    void Clear()
    {
        owners_.Clear();
        geometries_.Clear();
        worldBoundingSpheres_.Clear();
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
                worldBoundingSpheres_.Push(result.worldBoundingSpheres_);
            }
        }
        else
        {
            Swap(owners_, array[0].owners_);
            Swap(geometries_, array[0].geometries_);
            Swap(worldBoundingSpheres_, array[0].worldBoundingSpheres_);
        }
    }
};

class GeometryDrawableProcessor : public BasicDrawableProcessor
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
            GeometriesInFrustumQueryResult& result = threadResults[threadIndex];
            result.owners_.Push(QueryRangeOwner{ this, beginIndex, endIndex });
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                if (CheckDrawableIsGeometry(index, viewMask) && CheckDrawableInFrustum(index, frustum))
                {
                    result.geometries_.Push(data_.drawable_[index]);
                    result.worldBoundingSpheres_.Push(data_.worldBoundingSphere_[index]);
                    result.lightMasks_.Push(data_.lightMask_[index]);
                    result.zoneMasks_.Push(data_.zoneMask_[index]);
                }
            }
        });
    }
};

struct SceneFrustumQueryResult
{
    Vector<Drawable*> drawables_;
};

class StaticModelDrawableProcessor : public GeometryDrawableProcessor
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

class ZoneProcessor
{
public:
    ZoneProcessor(Zone* defaultZone) : defaultZone_(defaultZone), cameraZone_(defaultZone), farClipZone_(defaultZone) {}
    void ProcessQueryResult(ZonesInFrustumQueryResult& query, Camera* cullCamera)
    {
        Swap(zones_, query.zones_);
        highestZonePriority_ = M_MIN_INT;

        // Update zone priorities
        int bestPriority = M_MIN_INT;
        Node* cameraNode = cullCamera->GetNode();
        Vector3 cameraPos = cameraNode->GetWorldPosition();

        for (Zone* zone : zones_)
        {
            zones_.Push(zone);
            int priority = zone->GetPriority();
            if (priority > highestZonePriority_)
                highestZonePriority_ = priority;
            if (priority > bestPriority && zone->IsInside(cameraPos))
            {
                cameraZone_ = zone;
                bestPriority = priority;
            }
        }

        // Determine the zone at far clip distance. If not found, or camera zone has override mode, use camera zone
        cameraZoneOverride_ = cameraZone_->GetOverride();
        if (!cameraZoneOverride_)
        {
            Vector3 farClipPos = cameraPos + cameraNode->GetWorldDirection() * Vector3(0.0f, 0.0f, cullCamera->GetFarClip());
            bestPriority = M_MIN_INT;

            for (Zone* zone : zones_)
            {
                int priority = zone->GetPriority();
                if (priority > bestPriority && zone->IsInside(farClipPos))
                {
                    farClipZone_ = zone;
                    bestPriority = priority;
                }
            }
        }
        if (farClipZone_ == defaultZone_)
            farClipZone_ = cameraZone_;

    }
public:
    Zone* defaultZone_{};
    Vector<Zone*> zones_;
    int highestZonePriority_{};
    Zone* cameraZone_{};
    bool cameraZoneOverride_{};
    Zone* farClipZone_{};
};

struct LitGeometriesQueryResult
{
    /// Lit geometries.
    PODVector<Drawable*> litGeometries_;

    void Clear()
    {
        litGeometries_.Clear();
    }
};

class LightProcessor
{
public:
    void ProcessQueryResult(LightsInFrustumQueryResult& query,
        const FrameInfo& frameInfo, OcclusionBuffer* occlusionBuffer, Camera* cullCamera)
    {
        // Cull lights
        const unsigned numLights = query.lights_.Size();
        for (unsigned i = 0; i < numLights; ++i)
        {
            Light* light = query.lights_[i];
            if (occlusionBuffer && light->IsOccludee() && !occlusionBuffer->IsVisible(light->GetWorldBoundingBox()))
                continue;

            light->UpdateBatches(frameInfo);
            // If draw distance non-zero, update and check it
            float maxDistance = light->GetDrawDistance();
            if (maxDistance > 0.0f)
            {
                if (light->GetDistance() > maxDistance)
                    continue;
            }

            light->MarkInView(frameInfo);

            // Skip lights with zero brightness or black color
            if (!light->GetEffectiveColor().Equals(Color::BLACK))
                lights_.Push(light);
        }

        // Sort the lights to brightest/closest first, and per-vertex lights first so that per-vertex base pass can be evaluated first
        for (unsigned i = 0; i < lights_.Size(); ++i)
        {
            Light* light = lights_[i];
            light->SetIntensitySortValue(cullCamera->GetDistance(light->GetNode()->GetWorldPosition()));
            light->SetLightQueue(nullptr);
        }

        Sort(lights_.Begin(), lights_.End(), CompareLights);
    }
    void QueryLitGeometries(GeometriesInFrustumQueryResult& geometriesQuery, WorkQueue* workQueue)
    {
        litGeometries_.Resize(lights_.Size());
        ScheduleWork(workQueue, 2, lights_.Size(), 1,
            [this, &geometriesQuery](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            for (unsigned i = beginIndex; i < endIndex; ++i)
            {
                LitGeometriesQueryResult& result = litGeometries_[i];
                Light* light = lights_[i];
                const LightType type = light->GetLightType();
                const unsigned lightMask = light->GetLightMask();

                result.Clear();
//                 switch (type)
//                 {
//                 case LIGHT_DIRECTIONAL:
//                     const unsigned numGeometries = geometriesQuery.geometries_.Size();
//                     for (unsigned i = 0; i < numGeometries; ++i)
//                     {
//                         if (GetLightMask(geometries_[i]) & lightMask)
//                             query.litGeometries_.Push(geometries_[i]);
//                     }
//                     break;
//
//                 case LIGHT_SPOT:
//                     {
//                         FrustumOctreeQuery octreeQuery(tempDrawables, light->GetFrustum(), DRAWABLE_GEOMETRY,
//                             cullCamera_->GetViewMask());
//                         octree_->GetDrawables(octreeQuery);
//                         for (unsigned i = 0; i < tempDrawables.Size(); ++i)
//                         {
//                             if (tempDrawables[i]->IsInView(frame_) && (GetLightMask(tempDrawables[i]) & lightMask))
//                                 query.litGeometries_.Push(tempDrawables[i]);
//                         }
//                     }
//                     break;
//
//                 case LIGHT_POINT:
//                     {
//                         Sphere lightSphere(light->GetNode()->GetWorldPosition(), light->GetRange());
//                         SphereOctreeQuery octreeQuery(tempDrawables, lightSphere,
//                             DRAWABLE_GEOMETRY, cullCamera_->GetViewMask());
//                         octree_->GetDrawables(octreeQuery);
//                         for (unsigned i = 0; i < tempDrawables.Size(); ++i)
//                         {
//                             if (tempDrawables[i]->IsInView(frame_) && (GetLightMask(tempDrawables[i]) & lightMask))
//                                 query.litGeometries_.Push(tempDrawables[i]);
//                         }
//
//                         auto& geomsInFrustum = drawableProcessor_->GetGeometriesInFrustum();
//                         unsigned numGeometries = geomsInFrustum.geometries_.Size();
//                         for (unsigned i = 0; i < numGeometries; ++i)
//                         {
//                             if (lightSphere.IsInside(geomsInFrustum.boundingSpheres_[i]))
//                                 query.litGeometries_.Push(geomsInFrustum.geometries_[i]);
//                         }
//                     }
//                     break;
//                 }
            }
        });
    }

    const Vector<Light*>& GetLights() const { return lights_; }

private:
    Vector<Light*> lights_;
    Vector<LitGeometriesQueryResult> litGeometries_;
};

class DrawableProcessor
{
public:
    UniquePtr<SceneProcessor> sceneProcessor_;
    UniquePtr<ViewProcessor> viewProcessor_;

    DrawableProcessor()
    {
        sceneProcessor_ = MakeUnique<SceneProcessor>();
        viewProcessor_ = MakeUnique<ViewProcessor>();

        AddProcessor<GeometryDrawableProcessor>(StringHash());
        AddProcessor<ZoneDrawableProcessor, Zone>();
        AddProcessor<LightDrawableProcessor, Light>();
        AddProcessor<StaticModelDrawableProcessor, StaticModel>();
    }

    template <class T> void AddProcessor(StringHash type)
    {
        assert(!specializedProcessors_.Contains(type));
        UniquePtr<BasicDrawableProcessor> processor(new T());
        if (type)
            specializedProcessors_[type] = processor.Get();
        else
            defaultProcessor_ = processor.Get();

        if (auto zoneProcessor = dynamic_cast<ZoneDrawableProcessor*>(processor.Get()))
            zoneProcessors_.Push(zoneProcessor);
        if (auto lightProcessor = dynamic_cast<LightDrawableProcessor*>(processor.Get()))
            lightProcessors_.Push(lightProcessor);
        if (auto geometryProcessor = dynamic_cast<GeometryDrawableProcessor*>(processor.Get()))
            geometryProcessors_.Push(geometryProcessor);
        processors_.Push(std::move(processor));
    }

    template <class T, class U> void AddProcessor() { AddProcessor<T>(U::GetTypeStatic()); }

    void Update(WorkQueue* workQueue, View* view)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        Camera* cullCamera = view->GetCullCamera();
        Zone* defaultZone = view->GetRenderer()->GetDefaultZone();

        sceneProcessor_->UpdateDirtyDrawables();
        viewProcessor_->QuerySceneSpirits(workQueue, sceneProcessor_->GetData(),
            cullCamera->GetViewMask(), cullCamera->GetFrustum());
        viewProcessor_->CookZones(cullCamera, defaultZone);
        viewProcessor_->QuerySceneGeometries(workQueue, sceneProcessor_->GetData(), cullCamera, nullptr);
        viewProcessor_->UpdateGeometries(workQueue, view->GetFrameInfo());
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
    ZonesInFrustumQueryResult& GetZonesInFrustum() { return zonesInFrustumQuery_; }

    void QueryGeometryOccludersInFrustum(unsigned viewMask, const Frustum& frustum)
    {
        ClearThreadCache(geometryOccludersInFrustumQueryCache_);
        for (const auto& processor : geometryProcessors_)
            processor->QueryGeometryOccludersInFrustum(viewMask, frustum, workQueue_, numDrawablesInBucket_, geometryOccludersInFrustumQueryCache_);
        MergeThreadCache(geometryOccludersInFrustumQuery_, geometryOccludersInFrustumQueryCache_);
    }
    GeometryOccludersInFrustumQueryResult& GetGeometryOccludersInFrustum() { return geometryOccludersInFrustumQuery_; }

    void QueryLightsInFrustum(unsigned viewMask, const Frustum& frustum)
    {
        ClearThreadCache(lightsInFrustumQueryCache_);
        for (const auto& processor : lightProcessors_)
            processor->QueryLightsInFrustum(viewMask, frustum, workQueue_, numDrawablesInBucket_, lightsInFrustumQueryCache_);
        MergeThreadCache(lightsInFrustumQuery_, lightsInFrustumQueryCache_);
    }
    LightsInFrustumQueryResult& GetLightsInFrustum() { return lightsInFrustumQuery_; }

    void QueryGeometriesInFrustum(unsigned viewMask, const Frustum& frustum)
    {
        ClearThreadCache(geometriesInFrustumQueryCache_);
        for (const auto& processor : geometryProcessors_)
            processor->QueryGeometriesInFrustum(viewMask, frustum, workQueue_, numDrawablesInBucket_, geometriesInFrustumQueryCache_);
        MergeThreadCache(geometriesInFrustumQuery_, geometriesInFrustumQueryCache_);
    }
    GeometriesInFrustumQueryResult& GetGeometriesInFrustum() { return geometriesInFrustumQuery_; }

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
    Vector<ZoneDrawableProcessor*> zoneProcessors_;
    Vector<LightDrawableProcessor*> lightProcessors_;
    Vector<GeometryDrawableProcessor*> geometryProcessors_;
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
