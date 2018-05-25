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
#include "../Graphics/DrawableProcessor.h"
#include "../Container/HashMap.h"
// #include "../Container/List.h"
// #include "../Graphics/Batch.h"
// #include "../Graphics/Light.h"
// #include "../Graphics/Zone.h"
// #include "../Math/Polyhedron.h"

namespace Urho3D
{

struct ScenePassInfo;
class Renderer;
class WorkQueue;

struct LitGeometryDescIdx
{
    unsigned drawableIndex_{};
    float sortValue_{};
    unsigned short lightIndex_{};
    bool perVertex_{};
    bool negativeLight_{};
};

struct VisibleLightsPerThreadData
{
    /// Lit geometry.
    Vector<LitGeometryDescIdx> litGeometry_;
    /// Clear.
    void Clear()
    {
        litGeometry_.Clear();
    }
    /// Append.
    void Append(const VisibleLightsPerThreadData& other)
    {
        litGeometry_.Push(other.litGeometry_);
    }
};

struct LightBatchQueueEx : LightBatchQueue
{
    /// Whether the light is per-vertex.
    bool isPerVertex_{};
    /// Whether the light is shadowed.
    bool isShadowed_{};
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

struct BatchQueueData
{
    HashMap<BatchGroupKey, BatchGroup> batchGroups_;
    Vector<Batch> batches_;
    void ShallowClear();
    void AddBatch(const Batch& batch, bool allowInstancing);
    void AppendBatchGroups(const BatchQueueData& other);
    void ExportBatches(BatchQueue& queue);
    void ExportBatchGroups(BatchQueue& queue);
    static void ShallowClear(Vector<BatchQueueData*>& queues);
};

struct LightBatchQueueData
{
    BatchQueueData* litBaseQueueData_{};
    BatchQueueData* lightQueueData_{};
};

struct BatchCollectorPerThreadData
{
    /// Batch queues for each scene pass. Indexed via passIndex. Null for non-existing pass.
    Vector<BatchQueueData*> scenePassQueueData_;

    /// Persistent mapping for light batch queues data.
    HashMap<Light*, LightBatchQueueData> lightBatchQueueDataMap_;
    /// Batch queues for each light "litbase" pass. Indexed via lightIndex. Null for vertex lights.
    Vector<BatchQueueData*> litBaseQueueData_;
    /// Batch queues for each light "light" pass. Indexed via lightIndex. Null for vertex lights.
    Vector<BatchQueueData*> lightQueueData_;

    /// Clear light arrays.
    void ClearLightArrays(unsigned numLights)
    {
        litBaseQueueData_.Clear();
        litBaseQueueData_.Resize(numLights);
        lightQueueData_.Clear();
        lightQueueData_.Resize(numLights);
    }
    /// Perform shallow clear w/o deallocation.
    void ShallowClear()
    {
        BatchQueueData::ShallowClear(scenePassQueueData_);
        BatchQueueData::ShallowClear(litBaseQueueData_);
        BatchQueueData::ShallowClear(lightQueueData_);
    }
};

class BatchCollector : public Object
{
    URHO3D_OBJECT(BatchCollector, Object);

public:
    BatchCollector(Context* context);
    void Initialize(bool threading, const PODVector<ScenePassInfo>& scenePasses);

    /// Clear the state before processing the frame.
    // TODO(eugeneko) Pass sceneGrid here
    void Clear(Camera* cullCamera, const FrameInfo& frame);
    /// Collect zones and occluders.
    void CollectZonesAndOccluders(SceneGrid* sceneGrid);
    /// Process zones.
    void ProcessZones();
    /// Collect geometries and lights.
    void CollectGeometriesAndLights(SceneGrid* sceneGrid, OcclusionBuffer* occlusionBuffer);
    /// Update lights.
    void UpdateAndSortLights();
    /// Update visible geometries and shadow casters.
    void UpdateVisibleGeometriesAndShadowCasters();
    /// Process lights.
    void ProcessLights(SceneGrid* sceneGrid);
    /// Sort lit geometries.
    void SortLitGeometries(SceneGrid* sceneGrid);

    /// Add scene pass batch.
    void AddScenePassBatch(unsigned threadIndex, unsigned passIndex, const Batch& batch, bool grouped);
    /// Add light batch for "litbase" pass.
    void AddLitBaseBatch(unsigned threadIndex, unsigned lightIndex, const Batch& batch, bool grouped);
    /// Add light batch for "light" pass.
    void AddLightBatch(unsigned threadIndex, unsigned lightIndex, const Batch& batch, bool grouped);

    /// Merge threaded results.
    void MergeThreadedResults();
    /// Fill batch queues.
    void FillBatchQueues();
    /// Compose threaded results.
    // TODO(eugeneko) Remove this parameter
    void FinalizeBatches(unsigned alphaPassIndex);

    const Vector<Drawable*>& GetVisibleGeometries() const { return visibleGeometries_; }
    const Vector<unsigned>& GetVisibleGeometriesNumLights() const { return numLightsPerVisibleGeometry_; }
    const Vector<LitGeometryDescPacked>& GetLitGeometries() const { return sortedLitGeometries_; }

    BatchQueue* GetScenePassQueue(unsigned passIndex)
    {
        return scenePassQueues_[passIndex].Get();
    }
    LightBatchQueueEx* GetLightBatchQueue(unsigned lightIndex)
    {
        return lightBatchQueues_[lightIndex];
    }

    /// Return whether the calculations are threaded.
    bool IsThreaded() const { return threading_; }
    /// Return whether there are visible zones.
    bool HasVisibleZones() const { return !zonesAndOccluders_[0].zones_.Empty(); }
    /// Return all visible zones.
    const Vector<Zone*>& GetVisibleZones() const { return zonesAndOccluders_[0].zones_; }
    /// Return all visible lights.
    const Vector<Light*>& GetVisibleLights() const { return geometriesAndLights_[0].lights_; }
    /// Return visible light.
    Light* GetVisibleLight(unsigned lightIndex) const { return geometriesAndLights_[0].lights_[lightIndex]; }
    /// Return zone for specified drawable.
    Zone* GetDrawableZone(Drawable* drawable) const;
    /// Return far clip zone.
    Zone* GetFarClipZone() const { return zonesData_.farClipZone_; }

    /// Get all light batch queues. Indexed via lightIndex. Not null.
    const Vector<LightBatchQueueEx*>& GetLightBatchQueues() const { return lightBatchQueues_; }
    /// Return whether there is any light batch queue.
    bool HasLightBatchQueues() const { return !lightBatchQueues_.Empty(); }

private:
    BatchQueueData* AllocateStaticQueueData()
    {
        staticQueueDataPool_.Emplace(MakeUnique<BatchQueueData>());
        return staticQueueDataPool_.Back().Get();
    }
    BatchQueueData* AllocateDynamicQueueData()
    {
        dynamicQueueDataPool_.Emplace(MakeUnique<BatchQueueData>());
        return dynamicQueueDataPool_.Back().Get();
    }
    LightBatchQueueEx* AllocateLightBatchQueue()
    {
        lightBatchQueuePool_.Emplace(MakeUnique<LightBatchQueueEx>());
        return lightBatchQueuePool_.Back().Get();
    }
    /// Get or create light batch queue.
    LightBatchQueueEx* GetOrCreateLightBatchQueue(Light* light);

    void FinalizeBatchQueue(BatchQueue& queue, bool allowShadows);
    void FinalizeBatch(Batch& batch, bool allowShadows, const BatchQueue& queue);
    void FinalizeBatchGroup(BatchGroup& batchGroup, bool allowShadows, const BatchQueue& queue);

private:
    /// Clear vector per-element.
    template <class T> static void ClearVector(T& vector)
    {
        for (auto& item : vector)
            item.Clear();
    }
    /// Append all array elements to the first one.
    template <class T> static void AppendVectorToFirst(T& vector)
    {
        for (unsigned i = 1; i < vector.Size(); ++i)
            vector[0].Append(vector[i]);
    }
    /// Update drawable zone
    static void UpdateDirtyZone(const Vector<Zone*>& zones, SceneGridCellDrawableSoA& cellData, unsigned index,
        unsigned viewMask, const Frustum& frustum);
    /// Check whether the light has shadow.
    static bool IsLightShadowed(Light* light);

private:
    /// @name Constants since initialization
    /// @{

    Renderer* renderer_{};
    WorkQueue* workQueue_{};
    bool threading_{};
    unsigned maxScenePassIndex_{};
    unsigned numThreads_{};
    unsigned maxSortedInstances_{};

    /// @}

private:
    /// Per-frame constants
    /// @{

    Camera* cullCamera_{};
    unsigned viewMask_{};
    FrameInfo frame_{};
    unsigned frustumQueryThreadingThreshold_{ 32 };
    unsigned lightUpdateThreshold_{ 16 };
    unsigned geometryUpdateThreshold_{ 16 };

    /// @}

    /// Scene query for zones and occluders.
    SceneGridQueryResult zonesAndOccludersQuery_;
    /// Temporary buffer for visible zones and occluders. Size is equal to number of threads.
    Vector<SceneQueryZonesAndOccludersResult> zonesAndOccluders_;
    /// Processed zone data.
    ViewCookedZonesData zonesData_;
    /// Scene query for geometries and lights.
    SceneGridQueryResult geometriesAndLightsQuery_;
    /// Temporary buffer for visible geometries and lights. Size is equal to number of threads. Defines lightIndex.
    Vector<SceneQueryGeometriesAndLightsResult> geometriesAndLights_;
    /// Visible geometries. Ordered by gridIndex.
    Vector<Drawable*> visibleGeometries_;
    /// Temporary buffer for visible lights data. Size is equal to number of threads.
    Vector<VisibleLightsPerThreadData> lightsData_;

    /// Persistent mapping for light batch queues.
    HashMap<Light*, LightBatchQueueEx*> lightBatchQueueMap_;

    /// Light batch queues. Indexed via lightIndex. Not null.
    Vector<LightBatchQueueEx*> lightBatchQueues_;

    Vector<unsigned> numLightsPerVisibleGeometry_;
    Vector<LitGeometryDescPacked> sortedLitGeometries_;

    Vector<BatchCollectorPerThreadData> perThreadData_;

    Vector<UniquePtr<BatchQueue>> scenePassQueues_;

    Vector<UniquePtr<BatchQueueData>> staticQueueDataPool_;
    Vector<UniquePtr<BatchQueueData>> dynamicQueueDataPool_;
    Vector<UniquePtr<LightBatchQueueEx>> lightBatchQueuePool_;
};

}
