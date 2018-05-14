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
    Intersection IntersectBoxFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInside(worldBoundingBox_[index]); }
    bool IsInFrustum(unsigned index, const Frustum& frustum) const
    {
        const Intersection sphereIntersection = IntersectSphereFrustum(index, frustum);
        return sphereIntersection != OUTSIDE
            && (sphereIntersection == INSIDE || IntersectBoxFrustum(index, frustum) != OUTSIDE);
    }
    bool IsCenterInFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInside(worldBoundingSphere_[index].center_) != OUTSIDE; }
    bool IsOccludedByBuffer(unsigned index, OcclusionBuffer* buffer) const
    {
        return buffer && occludee_[index] && !buffer->IsVisible(worldBoundingBox_[index]);

    }

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

struct SceneQueryZonesAndOccludersResult
{
    Vector<Drawable*> occluders_;
    Vector<Zone*> zones_;
    void Clear()
    {
        occluders_.Clear();
        zones_.Clear();
    }
    void Merge(Vector<SceneQueryZonesAndOccludersResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (SceneQueryZonesAndOccludersResult& result : array)
            {
                occluders_.Push(result.occluders_);
                zones_.Push(result.zones_);
            }
        }
        else
        {
            Swap(occluders_, array[0].occluders_);
            Swap(zones_, array[0].zones_);
        }
    }
};

struct SceneQueryGeometriesAndLightsResult
{
    Vector<Light*> lights_;

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
        lights_.Clear();

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
    void Merge(Vector<SceneQueryGeometriesAndLightsResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (SceneQueryGeometriesAndLightsResult& result : array)
            {
                lights_.Push(result.lights_);

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
            Swap(lights_, array[0].lights_);

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
    Zone* cameraZone_{};
    Zone* farClipZone_{};
};

/// Contains view-specific scene data.
class ViewProcessor
{
public:
    virtual unsigned GetSceneQueryThreshold() { return 64; }
    virtual unsigned GetGeometryUpdateThreshold() { return 64; }
    virtual unsigned GetLightUpdateThreshold() { return 64; }
    /// Query zones and occluders.
    virtual void QuerySceneZonesAndOccluders(WorkQueue* workQueue, const SceneProcessorDrawableSoA& sceneData,
        unsigned viewMask, const Frustum& frustum)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        sceneZonesAndOccludersQueryThreadResults_.Resize(numThreads);
        for (SceneQueryZonesAndOccludersResult& result : sceneZonesAndOccludersQueryThreadResults_)
            result.Clear();

        workQueue->ScheduleWork(GetSceneQueryThreshold(), sceneData.size_, numThreads,
            [this, viewMask, frustum, &sceneData](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            SceneQueryZonesAndOccludersResult& result = sceneZonesAndOccludersQueryThreadResults_[threadIndex];
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
                else if (sceneData.IsGeometry(index) && sceneData.occluder_[index])
                {
                    if (sceneData.IsInFrustum(index, frustum))
                        result.occluders_.Push(drawable);
                }
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
        sceneZonesAndOccludersQueryResult_.Merge(sceneZonesAndOccludersQueryThreadResults_);
    }
    /// Prepare zones for further usage.
    virtual void CookZones(Camera* cullCamera, Zone* defaultZone)
    {
        Vector<Zone*>& zones = sceneZonesAndOccludersQueryResult_.zones_;
        zones_.cameraZoneOverride_ = false;
        zones_.cameraZone_ = nullptr;
        zones_.farClipZone_ = nullptr;

        // Sort zones
        Sort(zones.Begin(), zones.End(),
            [](Zone* lhs, Zone* rhs) { return lhs->GetPriority() > rhs->GetPriority(); });

        // Find camera and far clip zones
        Node* cullCameraNode = cullCamera->GetNode();
        const Vector3 cameraPos = cullCameraNode->GetWorldPosition();
        const Vector3 farClipPos = cameraPos + cullCameraNode->GetWorldDirection() * Vector3(0.0f, 0.0f, cullCamera->GetFarClip());
        for (Zone* zone : zones)
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
    virtual void QuerySceneGeometriesAndLights(WorkQueue* workQueue, SceneProcessorDrawableSoA& sceneData,
        Camera* cullCamera, OcclusionBuffer* occlusionBuffer)
    {
        const unsigned viewMask = cullCamera->GetViewMask();
        const Frustum& frustum = cullCamera->GetFrustum();
        const Matrix3x4& viewMatrix = cullCamera->GetView();
        const Vector3 viewZ = Vector3(viewMatrix.m20_, viewMatrix.m21_, viewMatrix.m22_);
        const Vector3 absViewZ = viewZ.Abs();

        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        sceneGeometriesAndLightsQueryThreadResults_.Resize(numThreads);
        for (SceneQueryGeometriesAndLightsResult& result : sceneGeometriesAndLightsQueryThreadResults_)
            result.Clear();

        workQueue->ScheduleWork(GetSceneQueryThreshold(), sceneData.size_, numThreads,
            [=, &sceneData](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            // TODO(eugeneko) Process lights here
            SceneQueryGeometriesAndLightsResult& result = sceneGeometriesAndLightsQueryThreadResults_[threadIndex];
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                // Discard drawables from other views
                if (!sceneData.MatchViewMask(index, viewMask))
                    continue;

                // Discard drawables except lights and geometries
                const bool isGeometry = sceneData.IsGeometry(index);
                const bool isLigth = sceneData.IsLight(index);
                if (!isGeometry && !isLigth)
                    continue;

                // Discard invisible
                if (!sceneData.IsInFrustum(index, frustum))
                    continue;

                // Calculate drawable distance
                const float drawDistance = sceneData.drawDistance_[index];
                const float distance = cullCamera->GetDistance(sceneData.worldBoundingSphere_[index].center_);

                // Discard if drawable is too far
                if (drawDistance > 0.0f && distance > drawDistance)
                    continue;

                // Discard using occlusion buffer
                if (sceneData.IsOccludedByBuffer(index, occlusionBuffer))
                    continue;

                // Update drawable zone
                Drawable* drawable = sceneData.drawable_[index];
                if (isGeometry)
                {
                    if (!zones_.cameraZoneOverride_ && !GetZones().Empty())
                    {
                        UpdateDirtyZone(sceneData, index, viewMask, frustum);
                    }

                    // Update min and max Z
                    float minZ = M_LARGE_VALUE;
                    float maxZ = M_LARGE_VALUE;

                    const BoundingBox& boundingBox = sceneData.worldBoundingBox_[index];
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
                    result.geometries_.Push(drawable);
                    result.zones_.Push(sceneData.cachedZone_[index]);
                    result.zoneLightMasks_.Push(sceneData.cachedZoneLightMask_[index]);
                    result.zoneShadowMasks_.Push(sceneData.cachedZoneShadowMask_[index]);
                    result.distances_.Push(distance);
                    result.minmaxZ_.Push(Vector2(minZ, maxZ));
                }
                else //if (isLight)
                {
                    result.lights_.Push(static_cast<Light*>(drawable));
                }
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
        sceneGeometriesAndLightsQueryResult_.Merge(sceneGeometriesAndLightsQueryThreadResults_);
    }
    /// Update visible lights.
    virtual void UpdateLights(WorkQueue* workQueue, const FrameInfo& frame)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;

        workQueue->ScheduleWork(GetLightUpdateThreshold(), sceneGeometriesAndLightsQueryResult_.lights_.Size(), numThreads,
            [this, &frame](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                // TODO(eugeneko) REMOVE!!
                Drawable* drawable = sceneGeometriesAndLightsQueryResult_.lights_[index];
                drawable->UpdateBatches(frame);
                drawable->MarkInView(frame);
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
    }
    /// Prepare lights for further usage.
    virtual void CookLights(Camera* cullCamera)
    {
        Vector<Light*> lights = sceneGeometriesAndLightsQueryResult_.lights_;
        // Sort the lights to brightest/closest first, and per-vertex lights first so that per-vertex base pass can be evaluated first
        for (Light* light : lights)
        {
            light->SetIntensitySortValue(cullCamera->GetDistance(light->GetNode()->GetWorldPosition()));
            light->SetLightQueue(nullptr);
        }

        Sort(lights.Begin(), lights.End(), CompareLights);
    }
    /// Update visible geometries.
    virtual void UpdateGeometries(WorkQueue* workQueue, const FrameInfo& frame)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;

        workQueue->ScheduleWork(GetGeometryUpdateThreshold(), sceneGeometriesAndLightsQueryResult_.geometries_.Size(), numThreads,
            [this, &frame](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                // TODO(eugeneko) REMOVE!!
                Drawable* drawable = sceneGeometriesAndLightsQueryResult_.geometries_[index];
                drawable->UpdateBatches(frame);
                drawable->SetZone(sceneGeometriesAndLightsQueryResult_.zones_[index]);
                drawable->MarkInView(frame);
                drawable->SetMinMaxZ(sceneGeometriesAndLightsQueryResult_.minmaxZ_[index].x_, sceneGeometriesAndLightsQueryResult_.minmaxZ_[index].y_);
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
    }

    const SceneQueryZonesAndOccludersResult& GetZonesAndOccluders() const { return sceneZonesAndOccludersQueryResult_; }
    const SceneQueryGeometriesAndLightsResult& GetGeometriesAndLights() const { return sceneGeometriesAndLightsQueryResult_; }
    const ViewCookedZonesData& GetCookedZonesInfo() const { return zones_; }
    const Vector<Zone*>& GetZones() const { return sceneZonesAndOccludersQueryResult_.zones_; }

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
            const Vector<Zone*>& zones = GetZones();
            const bool temporary = !sceneData.IsCenterInFrustum(index, frustum);
            const int highestZonePriority = zones.Empty() ? M_MIN_INT : zones[0]->GetPriority();

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
                for (Zone* zone : zones)
                {
                    if ((sceneData.zoneMask_[index] & zone->GetZoneMask()) && zone->IsInside(drawableCenter))
                    {
                        newZone = zone;
                        break;
                    }
                }
            }

            // Setup new zone
            sceneData.cachedZoneDirty_[index] = temporary;
            sceneData.cachedZone_[index] = newZone;
            sceneData.cachedZoneViewMask_[index] = newZone->GetViewMask();
            sceneData.cachedZoneLightMask_[index] = newZone->GetLightMask();
            sceneData.cachedZoneShadowMask_[index] = newZone->GetShadowMask();
        }
    }

private:
    Vector<SceneQueryZonesAndOccludersResult> sceneZonesAndOccludersQueryThreadResults_;
    SceneQueryZonesAndOccludersResult sceneZonesAndOccludersQueryResult_;

    ViewCookedZonesData zones_;

    Vector<SceneQueryGeometriesAndLightsResult> sceneGeometriesAndLightsQueryThreadResults_;
    SceneQueryGeometriesAndLightsResult sceneGeometriesAndLightsQueryResult_;


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
    }

    void Update(WorkQueue* workQueue, View* view)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        Camera* cullCamera = view->GetCullCamera();
        Zone* defaultZone = view->GetRenderer()->GetDefaultZone();

        sceneProcessor_->UpdateDirtyDrawables();
        viewProcessor_->QuerySceneZonesAndOccluders(workQueue, sceneProcessor_->GetData(),
            cullCamera->GetViewMask(), cullCamera->GetFrustum());
        viewProcessor_->CookZones(cullCamera, defaultZone);
        viewProcessor_->QuerySceneGeometriesAndLights(workQueue, sceneProcessor_->GetData(), cullCamera, nullptr);
        viewProcessor_->UpdateLights(workQueue, view->GetFrameInfo());
        viewProcessor_->CookLights(cullCamera);
        viewProcessor_->UpdateGeometries(workQueue, view->GetFrameInfo());
    }

#if 0
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
#endif
};

}
