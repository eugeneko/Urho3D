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
class Renderer;
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
    void Clear(unsigned frameNumber);
    /// Process lights. Lights shall be ordered according to lightIndex.
    void ProcessLights(const PODVector<Light*>& lights);
    /// Collect visible geometries.
    void CollectVisibleGeometry(SceneGridDrawableSoA& sceneData);
    /// Collect lit geometries from given unsorted array of lit geometries.
    void CollectLitGeometries(const Vector<LitGeometryDescIdx>& litGeometries, SceneGridDrawableSoA& sceneData);

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
    const Vector<LitGeometryDescPacked>& GetLitGeometries() const { return litGeometries_; }

    BatchQueue* GetScenePassQueue(unsigned passIndex)
    {
        return scenePassQueues_[passIndex].Get();
    }
    LightBatchQueue* GetLightBatchQueue(unsigned lightIndex)
    {
        return lightBatchQueues_[lightIndex];
    }

    /// Return whether the calculations are threaded.
    bool IsThreaded() const { return threading_; }
    /// Get all light batch queues. Indexed via lightIndex. Not null.
    const Vector<LightBatchQueue*>& GetLightBatchQueues() const { return lightBatchQueues_; }
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
    LightBatchQueue* AllocateLightBatchQueue()
    {
        lightBatchQueuePool_.Emplace(MakeUnique<LightBatchQueue>());
        return lightBatchQueuePool_.Back().Get();
    }
    void ProcessLight(unsigned lightIndex);
    void FinalizeBatchQueue(BatchQueue& queue, bool allowShadows);
    void FinalizeBatch(Batch& batch, bool allowShadows, const BatchQueue& queue);
    void FinalizeBatchGroup(BatchGroup& batchGroup, bool allowShadows, const BatchQueue& queue);

private:
    Renderer* renderer_{};
    WorkQueue* workQueue_{};
    bool threading_{};
    unsigned maxScenePassIndex_{};
    unsigned numThreads_{};
    unsigned maxSortedInstances_{};

    /// Persistent mapping for light batch queues.
    HashMap<Light*, LightBatchQueue*> lightBatchQueueMap_;

    /// Currently visible lights. Indexed via lightIndex. Not null.
    Vector<Light*> lights_;
    /// Light batch queues. Indexed via lightIndex. Not null.
    Vector<LightBatchQueue*> lightBatchQueues_;

    Vector<Drawable*> visibleGeometries_;
    Vector<unsigned> numLightsPerVisibleGeometry_;
    Vector<LitGeometryDescPacked> litGeometries_;

    Vector<BatchCollectorPerThreadData> perThreadData_;

    Vector<UniquePtr<BatchQueue>> scenePassQueues_;

    Vector<UniquePtr<BatchQueueData>> staticQueueDataPool_;
    Vector<UniquePtr<BatchQueueData>> dynamicQueueDataPool_;
    Vector<UniquePtr<LightBatchQueue>> lightBatchQueuePool_;
};

}
