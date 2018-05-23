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

#include "../Precompiled.h"

// #include "../Container/Sort.h"
#include "../Core/Context.h"
#include "../Core/WorkQueue.h"
#include "../Graphics/BatchCollector.h"
#include "../Graphics/Geometry.h"
#include "../Graphics/Renderer.h"

// TODO(eugeneko) Remove these dependencies
#include "../Graphics/View.h"
// #include "../Core/Profiler.h"
// #include "../Graphics/Camera.h"
// #include "../Graphics/DebugRenderer.h"
// #include "../Graphics/Graphics.h"
// #include "../Graphics/GraphicsEvents.h"
// #include "../Graphics/GraphicsImpl.h"
// #include "../Graphics/Material.h"
// #include "../Graphics/OcclusionBuffer.h"
// #include "../Graphics/Octree.h"
// #include "../Graphics/RenderPath.h"
// #include "../Graphics/ShaderVariation.h"
// #include "../Graphics/Skybox.h"
// #include "../Graphics/Technique.h"
// #include "../Graphics/Texture2D.h"
// #include "../Graphics/Texture2DArray.h"
// #include "../Graphics/Texture3D.h"
// #include "../Graphics/TextureCube.h"
// #include "../Graphics/VertexBuffer.h"
// #include "../Graphics/View.h"
// #include "../IO/FileSystem.h"
// #include "../IO/Log.h"
// #include "../Resource/ResourceCache.h"
// #include "../Scene/Scene.h"
// #include "../UI/UI.h"

#include "../DebugNew.h"

namespace Urho3D
{

void BatchQueueData::ShallowClear()
{
    batches_.Clear();
    for (auto& group : batchGroups_)
    {
        group.second_.Clear();
    }
}

void BatchQueueData::AddBatch(const Batch& batch, bool allowInstancing)
{
    if (allowInstancing)
    {
        BatchGroupKey key(batch);

        // Find new group or spawn empty
        auto iterGroup = batchGroups_.Find(key);
        if (iterGroup == batchGroups_.End())
            iterGroup = batchGroups_.Insert(MakePair(key, BatchGroup()));

        // Reset batch if group was empty
        BatchGroup& group = iterGroup->second_;
        if (group.IsEmpty())
            group.ResetBatch(batch);

        // Add instance
        group.AddTransforms(batch);
    }
    else
    {
        // If batch is static with multiple world transforms and cannot instance, we must push copies of the batch individually
        if (batch.geometryType_ == GEOM_STATIC && batch.numWorldTransforms_ > 1)
        {
            Batch batchCopy = batch;
            batchCopy.numWorldTransforms_ = 1;
            for (unsigned i = 0; i < batch.numWorldTransforms_; ++i)
            {
                batchCopy.worldTransform_ = &batch.worldTransform_[i];
                batches_.Push(batchCopy);
            }
        }
        else
            batches_.Push(batch);

    }
}

void BatchQueueData::AppendBatchGroups(const BatchQueueData& other)
{
    for (const auto& otherGroup : other.batchGroups_)
    {
        // Skip empty groups
        if (otherGroup.second_.IsEmpty())
            continue;

        const BatchGroupKey& key = otherGroup.first_;
        auto iterGroup = batchGroups_.Find(key);
        if (iterGroup == batchGroups_.End())
        {
            // Copy entire group
            iterGroup = batchGroups_.Insert(MakePair(key, BatchGroup()));
            iterGroup->second_ = otherGroup.second_;
        }
        else
        {
            // Copy instances only
            iterGroup->second_.instances_.Push(otherGroup.second_.instances_);
        }
    }
}

void BatchQueueData::ExportBatches(BatchQueue& queue)
{
    for (Batch& batch : batches_)
        queue.sortedBatches_.Push(&batch);
}

void BatchQueueData::ExportBatchGroups(BatchQueue& queue)
{
    for (auto& batchGroup : batchGroups_)
        queue.sortedBatchGroups_.Push(&batchGroup.second_);
}

void BatchQueueData::ShallowClear(Vector<BatchQueueData*>& queues)
{
    for (BatchQueueData* queue : queues)
    {
        if (queue)
            queue->ShallowClear();
    }
}

//////////////////////////////////////////////////////////////////////////
BatchCollector::BatchCollector(Context* context)
    : Object(context)
    , renderer_(context->GetSubsystem<Renderer>())
    , workQueue_(context->GetSubsystem<WorkQueue>())
{

}

void BatchCollector::Initialize(bool threading, const PODVector<ScenePassInfo>& scenePasses)
{
    // Calculate constants
    threading_ = threading && workQueue_->GetNumThreads() > 0;
    numThreads_ = threading_ ? workQueue_->GetNumThreads() + 1 : 1;
    perThreadData_.Resize(numThreads_);

    maxSortedInstances_ = static_cast<unsigned>(renderer_->GetMaxSortedInstances());

    // Allocate query storages
    zonesAndOccluders_.Resize(numThreads_);
    geometriesAndLights_.Resize(numThreads_);

    // Calculate max scene pass
    maxScenePassIndex_ = 0;
    for (const ScenePassInfo& info : scenePasses)
        maxScenePassIndex_ = Max(maxScenePassIndex_, info.passIndex_);

    // Create batch queues for scene passes
    scenePassQueues_.Clear();
    scenePassQueues_.Resize(maxScenePassIndex_ + 1);
    for (const ScenePassInfo& info : scenePasses)
        scenePassQueues_[info.passIndex_] = MakeUnique<BatchQueue>();

    // Allocate queue data for scene passes
    staticQueueDataPool_.Clear();
    for (unsigned threadIndex = 0; threadIndex < numThreads_; ++threadIndex)
    {
        BatchCollectorPerThreadData& threadData = perThreadData_[threadIndex];
        threadData.scenePassQueueData_.Clear();
        threadData.scenePassQueueData_.Resize(maxScenePassIndex_ + 1);
        for (const ScenePassInfo& info : scenePasses)
            threadData.scenePassQueueData_[info.passIndex_] = AllocateStaticQueueData();
    }
}

void BatchCollector::Clear(Camera* cullCamera, unsigned frameNumber)
{
    cullCamera_ = cullCamera;
    viewMask_ = cullCamera->GetViewMask();

    lights_.Clear();
    visibleGeometries_.Clear();
    litGeometries_.Clear();
    for (BatchCollectorPerThreadData& perThread : perThreadData_)
        perThread.ShallowClear();
}

void BatchCollector::CollectZonesAndOccluders(SceneGrid* sceneGrid)
{
    const Frustum frustum = cullCamera_->GetFrustum();

    sceneGrid->QueryCellsInFrustum(frustum, frustumQueryThreadingThreshold_, numThreads_, zonesAndOccludersQuery_);
    ClearVector(zonesAndOccluders_);

    zonesAndOccludersQuery_.ScheduleWork(workQueue_,
        [=](SceneGridCellRef& cellRef, unsigned threadIndex)
    {
        SceneGridCellDrawableSoA& cell = *cellRef.data_;
        SceneQueryZonesAndOccludersResult& result = zonesAndOccluders_[threadIndex];
        for (unsigned index = cellRef.beginIndex_; index < cellRef.endIndex_; ++index)
        {
            if (!cell.MatchViewMask(index, viewMask_))
                continue;

            Drawable* drawable = cell.drawable_[index];

            // TODO(eugeneko) Get rid of branching
            if (cell.IsZone(index))
            {
                result.zones_.Push(static_cast<Zone*>(drawable));
            }
            else if (cell.IsGeometry(index) && cell.occluder_[index])
            {
                if (cellRef.intersection_ == INSIDE || cell.IsInFrustum(index, frustum))
                    result.occluders_.Push(drawable);
            }
        }
    });

    workQueue_->Complete(M_MAX_UNSIGNED);
    AppendVectorToFirst(zonesAndOccluders_);
}

void BatchCollector::ProcessZones()
{
    Zone* defaultZone = renderer_->GetDefaultZone();

    ZoneVector& zones = zonesAndOccluders_[0].zones_;
    zonesData_.cameraZoneOverride_ = false;
    zonesData_.cameraZone_ = nullptr;
    zonesData_.farClipZone_ = nullptr;

    // Sort zones
    Sort(zones.Begin(), zones.End(),
        [](Zone* lhs, Zone* rhs) { return lhs->GetPriority() > rhs->GetPriority(); });

    // Find camera and far clip zones
    Node* cullCameraNode = cullCamera_->GetNode();
    const Vector3 cameraPos = cullCameraNode->GetWorldPosition();
    const Vector3 farClipPos = cameraPos + cullCameraNode->GetWorldDirection() * Vector3(0.0f, 0.0f, cullCamera_->GetFarClip());
    for (Zone* zone : zones)
    {
        if (!zonesData_.cameraZone_ && zone->IsInside(cameraPos))
            zonesData_.cameraZone_ = zone;
        if (!zonesData_.farClipZone_ && zone->IsInside(farClipPos))
            zonesData_.farClipZone_ = zone;
        if (zonesData_.cameraZone_ && zonesData_.farClipZone_)
            break;
    }

    if (!zonesData_.cameraZone_)
        zonesData_.cameraZone_ = defaultZone;

    if (!zonesData_.farClipZone_)
        zonesData_.farClipZone_ = defaultZone;
}

void BatchCollector::CollectGeometriesAndLights(SceneGrid* sceneGrid, OcclusionBuffer* occlusionBuffer)
{
    const Frustum frustum = cullCamera_->GetFrustum();

    sceneGrid->QueryCellsInFrustum(frustum, frustumQueryThreadingThreshold_, numThreads_, geometriesAndLightsQuery_);

    const unsigned viewMask = cullCamera_->GetViewMask();
    const Matrix3x4 viewMatrix = cullCamera_->GetView();
    const Vector3 viewZ = Vector3(viewMatrix.m20_, viewMatrix.m21_, viewMatrix.m22_);
    const Vector3 absViewZ = viewZ.Abs();
    const unsigned cameraZoneLightMask = zonesData_.cameraZone_->GetLightMask();
    const unsigned cameraZoneShadowMask = zonesData_.cameraZone_->GetShadowMask();

    ClearVector(geometriesAndLights_);

    zonesAndOccludersQuery_.ScheduleWork(workQueue_,
        [=](SceneGridCellRef& cellRef, unsigned threadIndex)
    {
        SceneGridCellDrawableSoA& cell = *cellRef.data_;
        SceneQueryGeometriesAndLightsResult& result = geometriesAndLights_[threadIndex];

        // TODO(eugeneko) Process lights here
        for (unsigned index = cellRef.beginIndex_; index < cellRef.endIndex_; ++index)
        {
            Vector<bool>& resultVisible = sceneGrid->GetGlobalData().visible_;

            // Discard drawables from other views
            if (!cell.MatchViewMask(index, viewMask))
                continue;

            // Discard drawables except lights and geometries
            const bool isGeometry = cell.IsGeometry(index);
            const bool isLigth = cell.IsLight(index);
            if (!isGeometry && !isLigth)
                continue;

            // Discard invisible
            if (cellRef.intersection_ == INTERSECTS && !cell.IsInFrustum(index, frustum))
                continue;

            // Calculate drawable distance
            const float drawDistance = cell.drawDistance_[index];
            const float distance = cullCamera_->GetDistance(cell.boundingSphere_[index].center_);

            // Discard if drawable is too far
            if (drawDistance > 0.0f && distance > drawDistance)
                continue;

            // Discard using occlusion buffer
            if (cell.IsOccludedByBuffer(index, occlusionBuffer))
                continue;

            // Update drawable zone
            Drawable* drawable = cell.drawable_[index];
            if (isGeometry)
            {
                if (!zonesData_.cameraZoneOverride_ && HasActiveZones())
                {
                    UpdateDirtyZone(GetActiveZones(), cell, index, viewMask, frustum);
                }

                // Update min and max Z
                float minZ = M_LARGE_VALUE;
                float maxZ = M_LARGE_VALUE;

                const BoundingBox& boundingBox = cell.boundingBox_[index];
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

                // Get actual zone
                Zone* cachedZone = cell.cachedZone_[index];
                Zone* actualZone = cachedZone;
                unsigned zoneLightMask = cell.cachedZoneLightMask_[index];
                unsigned zoneShadowMask = cell.cachedZoneShadowMask_[index];

                if (zonesData_.cameraZoneOverride_ || !cachedZone)
                {
                    actualZone = zonesData_.cameraZone_;
                    zoneLightMask = cameraZoneLightMask;
                    zoneShadowMask = cameraZoneShadowMask;
                }

                // Get masks
                const unsigned drawableLightMask = cell.lightMask_[index];
                const unsigned drawableShadowMask = cell.shadowMask_[index];

                // Mask object as visible

                const unsigned gridIndex = drawable->GetDrawableIndex().gridIndex_;
                resultVisible[gridIndex] = true;
                // TODO(eugeneko) Get rid of drawable access
            }
            else //if (isLight)
            {
                result.lights_.Push(static_cast<Light*>(drawable));
            }
        }
    });

    workQueue_->Complete(M_MAX_UNSIGNED);

    AppendVectorToFirst(geometriesAndLights_);
    if (geometriesAndLights_[0].minZ_ == M_INFINITY)
        geometriesAndLights_[0].minZ_ = 0.0f;
}

void BatchCollector::ProcessLights(const PODVector<Light*>& lights)
{
    // Assign lights
    lights_.Insert(lights_.Begin(), lights.Begin(), lights.End());
    const unsigned numLights = lights_.Size();

    // Resize related arrays
    lightBatchQueues_.Resize(numLights);
    for (BatchCollectorPerThreadData& perThread : perThreadData_)
        perThread.ClearLightArrays(numLights);

    // Process each light individually
    for (unsigned lightIndex = 0; lightIndex < numLights; ++lightIndex)
        ProcessLight(lightIndex);
}

void BatchCollector::CollectVisibleGeometry(SceneGridDrawableSoA& sceneData)
{
    assert(sceneData.IsValid());

    visibleGeometries_.Clear();
    numLightsPerVisibleGeometry_.Clear();

    for (unsigned i = 0; i < sceneData.size_; ++i)
    {
        if (sceneData.visible_[i])
            visibleGeometries_.Push(sceneData.drawable_[i]);
    }
}

void BatchCollector::CollectLitGeometries(const Vector<LitGeometryDescIdx>& litGeometries, SceneGridDrawableSoA& sceneData)
{
    litGeometries_.Clear();
    sceneData.ClearTemporary();

    // 1st pass: calculate number of lights for each geometry
    for (const LitGeometryDescIdx& lit : litGeometries)
        ++sceneData.numLights_[lit.drawableIndex_];

    // 2nd pass: calculate ranges in destination array
    unsigned firstLight = 0;
    for (unsigned i = 0; i < sceneData.size_; ++i)
    {
        sceneData.firstLight_[i] = firstLight;
        firstLight += sceneData.numLights_[i];
    }

    // 3rd pass: reset number of lights to reuse it for filling array
    for (unsigned i = 0; i < sceneData.size_; ++i)
        sceneData.numLights_[i] = 0;

    // 4th pass: fill destination array
    litGeometries_.Resize(litGeometries.Size());
    for (const LitGeometryDescIdx& lit : litGeometries)
    {
        const unsigned drawableIndex = lit.drawableIndex_;
        const unsigned litIndex = sceneData.firstLight_[drawableIndex] + sceneData.numLights_[drawableIndex];
        litGeometries_[litIndex] = LitGeometryDescPacked(lit);
        ++sceneData.numLights_[drawableIndex];
    }

    // 5th pass: partially sort lit geometries
    for (unsigned i = 0; i < sceneData.size_; ++i)
    {
        if (sceneData.visible_[i])
        {
            const unsigned numLights = sceneData.numLights_[i];
            numLightsPerVisibleGeometry_.Push(numLights);
            if (numLights > 1)
            {
                auto begin = litGeometries_.Begin() + sceneData.firstLight_[i];
                auto end = begin + numLights;
                Sort(begin, end);
            }
        }
    }
}

void BatchCollector::AddScenePassBatch(unsigned threadIndex, unsigned passIndex, const Batch& batch, bool grouped)
{
    BatchQueueData* queueData = perThreadData_[threadIndex].scenePassQueueData_[passIndex];
    assert(queueData);
    queueData->AddBatch(batch, grouped);
}

void BatchCollector::AddLitBaseBatch(unsigned threadIndex, unsigned lightIndex, const Batch& batch, bool grouped)
{
    BatchQueueData* queueData = perThreadData_[threadIndex].litBaseQueueData_[lightIndex];
    assert(queueData);
    queueData->AddBatch(batch, grouped);
}

void BatchCollector::AddLightBatch(unsigned threadIndex, unsigned lightIndex, const Batch& batch, bool grouped)
{
    BatchQueueData* queueData = perThreadData_[threadIndex].lightQueueData_[lightIndex];
    assert(queueData);
    queueData->AddBatch(batch, grouped);
}

void BatchCollector::MergeThreadedResults()
{
    if (!threading_)
        return;

    // Append all batch groups for scene passes
    for (unsigned passIndex = 0; passIndex < maxScenePassIndex_; ++passIndex)
    {
        BatchQueueData* destData = perThreadData_[0].scenePassQueueData_[passIndex];
        if (!destData)
            continue;

        for (unsigned threadIndex = 1; threadIndex < numThreads_; ++threadIndex)
        {
            BatchQueueData* sourceData = perThreadData_[threadIndex].scenePassQueueData_[passIndex];
            assert(sourceData);
            destData->AppendBatchGroups(*sourceData);
        }
    }

    // Append all batch groups for light passes.
    for (unsigned lightIndex = 0; lightIndex < lights_.Size(); ++lightIndex)
    {
        BatchQueueData* destLitBaseData = perThreadData_[0].litBaseQueueData_[lightIndex];
        BatchQueueData* destLightData = perThreadData_[0].lightQueueData_[lightIndex];
        assert(!!destLitBaseData == !!destLightData);
        if (!destLightData || !destLitBaseData)
            continue;

        for (unsigned threadIndex = 1; threadIndex < numThreads_; ++threadIndex)
        {
            BatchQueueData* sourceLitBaseData = perThreadData_[threadIndex].litBaseQueueData_[lightIndex];
            BatchQueueData* sourceLightData = perThreadData_[threadIndex].lightQueueData_[lightIndex];
            assert(sourceLitBaseData && sourceLightData);
            destLitBaseData->AppendBatchGroups(*sourceLitBaseData);
            destLightData->AppendBatchGroups(*sourceLightData);
        }
    }
}

void BatchCollector::FillBatchQueues()
{
    // Fill batch queues for scene passes
    for (unsigned passIndex = 0; passIndex < maxScenePassIndex_; ++passIndex)
    {
        BatchQueue* queue = scenePassQueues_[passIndex].Get();
        if (!queue)
            continue;

        // Clear batch queue
        queue->Clear(maxSortedInstances_);

        // Copy pointers to batch groups from first threaded storage
        BatchQueueData* mergedQueueData = perThreadData_[0].scenePassQueueData_[passIndex];
        assert(mergedQueueData);
        mergedQueueData->ExportBatchGroups(*queue);

        // Copy pointers to batches from each threaded storage
        for (unsigned threadIndex = 0; threadIndex < numThreads_; ++threadIndex)
        {
            BatchQueueData* queueData = perThreadData_[threadIndex].scenePassQueueData_[passIndex];
            assert(queueData);
            queueData->ExportBatches(*queue);
        }
    }

    // Fill batch queues for lights
    for (unsigned lightIndex = 0; lightIndex < lights_.Size(); ++lightIndex)
    {
        Light* light = lights_[lightIndex];
        if (light->GetPerVertex())
            continue;

        LightBatchQueue* lightBatchQueue = lightBatchQueues_[lightIndex];
        assert(lightBatchQueue);

        BatchQueue& lightQueue = lightBatchQueue->litBatches_;
        BatchQueue& litBaseQueue = lightBatchQueue->litBaseBatches_;

        // Clear batch queues
        lightQueue.Clear(maxSortedInstances_);
        litBaseQueue.Clear(maxSortedInstances_);

        // Copy pointers to batch groups from first threaded storage
        BatchQueueData* mergedLightQueueData = perThreadData_[0].lightQueueData_[lightIndex];
        BatchQueueData* mergedLitBaseQueueData = perThreadData_[0].litBaseQueueData_[lightIndex];
        assert(mergedLightQueueData && mergedLitBaseQueueData);
        mergedLightQueueData->ExportBatchGroups(lightQueue);
        mergedLitBaseQueueData->ExportBatchGroups(litBaseQueue);

        // Copy pointers to batches from each threaded storage
        for (unsigned threadIndex = 0; threadIndex < numThreads_; ++threadIndex)
        {
            BatchQueueData* lightQueueData = perThreadData_[threadIndex].lightQueueData_[lightIndex];
            BatchQueueData* litBaseQueueData = perThreadData_[threadIndex].litBaseQueueData_[lightIndex];
            assert(lightQueueData && litBaseQueueData);
            lightQueueData->ExportBatches(lightQueue);
            litBaseQueueData->ExportBatches(litBaseQueue);
        }
    }
}

void BatchCollector::FinalizeBatches(unsigned alphaPassIndex)
{
    const bool reuseShadowMaps = renderer_->GetReuseShadowMaps();

    // Update batch queues for scene passes
    for (unsigned passIndex = 0; passIndex < maxScenePassIndex_; ++passIndex)
    {
        BatchQueue* queue = scenePassQueues_[passIndex].Get();
        if (!queue)
            continue;

        const bool allowShadows = passIndex != alphaPassIndex || !reuseShadowMaps;
        FinalizeBatchQueue(*queue, allowShadows);
    }

    // Update batch queues for lights
    for (unsigned lightIndex = 0; lightIndex < lightBatchQueues_.Size(); ++lightIndex)
    {
        Light* light = lights_[lightIndex];
        if (light->GetPerVertex())
            continue;

        LightBatchQueue* lightBatchQueue = lightBatchQueues_[lightIndex];
        assert(lightBatchQueue);
        FinalizeBatchQueue(lightBatchQueue->litBatches_, true);
        FinalizeBatchQueue(lightBatchQueue->litBaseBatches_, true);
    }
}

void BatchCollector::ProcessLight(unsigned lightIndex)
{
    Light* light = lights_[lightIndex];

    // Find or create light batch queue
    auto iterLightBatchQueue = lightBatchQueueMap_.Find(light);
    if (iterLightBatchQueue == lightBatchQueueMap_.End())
    {
        LightBatchQueue* queue = AllocateLightBatchQueue();

        queue->light_ = light;
        queue->negative_ = light->IsNegative();
        queue->shadowMap_ = nullptr;

        iterLightBatchQueue = lightBatchQueueMap_.Insert(MakePair(light, queue));
    }

    // Assign queue
    lightBatchQueues_[lightIndex] = iterLightBatchQueue->second_;

    // Find or create queues data for each thread
    if (!light->GetPerVertex())
    {
        for (BatchCollectorPerThreadData& perThread : perThreadData_)
        {
            auto iterQueueData = perThread.lightBatchQueueDataMap_.Find(light);
            if (iterQueueData == perThread.lightBatchQueueDataMap_.End())
            {
                LightBatchQueueData desc;
                desc.lightQueueData_ = AllocateDynamicQueueData();
                desc.litBaseQueueData_ = AllocateDynamicQueueData();
                iterQueueData = perThread.lightBatchQueueDataMap_.Insert(MakePair(light, desc));
            }

            perThread.litBaseQueueData_[lightIndex] = iterQueueData->second_.litBaseQueueData_;
            perThread.lightQueueData_[lightIndex] = iterQueueData->second_.lightQueueData_;
        }
    }
}

void BatchCollector::FinalizeBatchQueue(BatchQueue& queue, bool allowShadows)
{
    for (Batch* batch : queue.sortedBatches_)
        FinalizeBatch(*batch, allowShadows, queue);
    for (BatchGroup* batchGroup : queue.sortedBatchGroups_)
        FinalizeBatchGroup(*batchGroup, allowShadows, queue);
}

void BatchCollector::FinalizeBatch(Batch& batch, bool allowShadows, const BatchQueue& queue)
{
    if (!batch.material_)
        batch.material_ = renderer_->GetDefaultMaterial();

    renderer_->SetBatchShaders(batch, nullptr, allowShadows, queue);
    batch.CalculateSortKey();
}

void BatchCollector::FinalizeBatchGroup(BatchGroup& batchGroup, bool allowShadows, const BatchQueue& queue)
{
    const int minInstances_ = renderer_->GetMinInstances();

    // Convert to instanced if possible
    if (batchGroup.geometryType_ == GEOM_STATIC && batchGroup.geometry_->GetIndexBuffer()
        && (int)batchGroup.instances_.Size() >= minInstances_)
    {
        batchGroup.geometryType_ = GEOM_INSTANCED;
    }

    // Finalize as batch
    FinalizeBatch(batchGroup, allowShadows, queue);
}

void BatchCollector::UpdateDirtyZone(const Vector<Zone*>& zones, SceneGridCellDrawableSoA& cellData, unsigned index,
    unsigned viewMask, const Frustum& frustum)
{
    const Vector3 drawableCenter = cellData.boundingSphere_[index].center_;

    Zone* cachedZone = cellData.cachedZone_[index];
    const bool cachedZoneDirty = cellData.cachedZoneDirty_[index];
    const unsigned cachedZoneViewMask = cellData.cachedZoneViewMask_[index];

    // TODO(eugeneko) Is branch is not optimized for multiple zones
    if (cachedZoneDirty || !cachedZone || !(cachedZoneViewMask & viewMask))
    {
        // Find new zone
        const bool temporary = !cellData.IsCenterInFrustum(index, frustum);
        const int highestZonePriority = zones.Empty() ? M_MIN_INT : zones[0]->GetPriority();

        Zone* newZone = nullptr;

        // First check if the current zone remains a conclusive result
        if (cachedZone && (cachedZone->GetViewMask() & viewMask)
            && cachedZone->GetPriority() >= highestZonePriority
            && (cellData.zoneMask_[index] & cachedZone->GetZoneMask())
            && cachedZone->IsInside(drawableCenter))
        {
            newZone = cachedZone;
        }
        else
        {
            // Search for appropriate zone
            for (Zone* zone : zones)
            {
                if ((cellData.zoneMask_[index] & zone->GetZoneMask()) && zone->IsInside(drawableCenter))
                {
                    newZone = zone;
                    break;
                }
            }
        }

        // Setup new zone
        cellData.cachedZoneDirty_[index] = temporary;
        cellData.cachedZone_[index] = newZone;
        cellData.cachedZoneViewMask_[index] = newZone->GetViewMask();
        cellData.cachedZoneLightMask_[index] = newZone->GetLightMask();
        cellData.cachedZoneShadowMask_[index] = newZone->GetShadowMask();
    }
}

}
