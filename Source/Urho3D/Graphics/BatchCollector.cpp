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
#include "../Graphics/Technique.h"
#include "../Graphics/Texture2D.h"
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
    {
        if (!batchGroup.second_.IsEmpty())
            queue.sortedBatchGroups_.Push(&batchGroup.second_);
    }
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
    lightsData_.Resize(numThreads_);

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

void BatchCollector::Clear(Camera* cullCamera, const FrameInfo& frame)
{
    cullCamera_ = cullCamera;
    viewMask_ = cullCamera->GetViewMask();
    frame_ = frame;

    visibleGeometries_.Clear();
    sortedLitGeometries_.Clear();
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

    geometriesAndLightsQuery_.ScheduleWork(workQueue_,
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
                if (!zonesData_.cameraZoneOverride_ && HasVisibleZones())
                {
                    UpdateDirtyZone(GetVisibleZones(), cell, index, viewMask, frustum);
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

                // Write result
                const unsigned gridIndex = drawable->GetDrawableIndex().gridIndex_;
                resultVisible[gridIndex] = true;
                cell.visible_[index] = true;
                cell.minmaxZ_[index] = Vector2(minZ, maxZ);
                // TODO(eugeneko) Get rid of drawable access
                drawable->SetZone(actualZone);
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
    geometriesAndLights_[0].maxZ_ = Max(geometriesAndLights_[0].minZ_ + M_LARGE_VALUE, geometriesAndLights_[0].maxZ_);

    // Collect visible geometries
    visibleGeometries_.Clear();
    numLightsPerVisibleGeometry_.Clear();
    SceneGridDrawableSoA& globalData = sceneGrid->GetGlobalData();
    for (unsigned gridIndex = 0; gridIndex < globalData.size_; ++gridIndex)
    {
        if (globalData.visible_[gridIndex])
            visibleGeometries_.Push(globalData.drawable_[gridIndex]);
    }
}

void BatchCollector::UpdateAndSortLights()
{
    // Update lights
    workQueue_->ScheduleWork(lightUpdateThreshold_, geometriesAndLights_[0].lights_.Size(), numThreads_,
        [=](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
    {
        for (unsigned index = beginIndex; index < endIndex; ++index)
        {
            // TODO(eugeneko) Refactor me
            Drawable* light = geometriesAndLights_[0].lights_[index];
            light->UpdateBatches(frame_);
            light->MarkInView(frame_);
        }
    });

    workQueue_->Complete(M_MAX_UNSIGNED);

    // Sort the lights to brightest/closest first, and per-vertex lights first so that per-vertex base pass can be evaluated first
    LightVector& lights = geometriesAndLights_[0].lights_;
    for (Light* light : lights)
    {
        light->SetIntensitySortValue(cullCamera_->GetDistance(light->GetNode()->GetWorldPosition()));
    }
    Sort(lights.Begin(), lights.End(), CompareLights);
}

void BatchCollector::UpdateVisibleGeometriesAndShadowCasters()
{
    // TODO(eugeneko) Update shadow casters
    workQueue_->ScheduleWork(geometryUpdateThreshold_, visibleGeometries_.Size(), numThreads_,
        [=](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
    {
        for (unsigned index = beginIndex; index < endIndex; ++index)
        {
            Drawable* drawable = visibleGeometries_[index];
            drawable->UpdateBatches(frame_);
            drawable->MarkInView(frame_);
        }
    });

    workQueue_->Complete(M_MAX_UNSIGNED);
}

void BatchCollector::ProcessLights(SceneGrid* sceneGrid)
{
    // Assign lights
    const unsigned numLights = GetVisibleLights().Size();

    // Resize related arrays
    lightBatchQueues_.Resize(numLights);
    for (BatchCollectorPerThreadData& perThread : perThreadData_)
        perThread.ClearLightArrays(numLights);

    ClearVector(lightsData_, numLights);

    // Step 1: Collect lit geometries and lit bounding boxes
    for (unsigned lightIndex = 0; lightIndex < numLights; ++lightIndex)
    {
        // Setup light
        Light* light = GetVisibleLight(lightIndex);
        const LightType lightType = light->GetLightType();
        const unsigned lightMask = light->GetLightMask();

        LightBatchQueueEx* lightBatchQueue = GetOrCreateLightBatchQueue(light);
        lightBatchQueues_[lightIndex] = lightBatchQueue;

        lightBatchQueue->Clear();
        lightBatchQueue->negative_ = light->IsNegative();
        lightBatchQueue->shadowMap_ = nullptr;
        lightBatchQueue->isShadowed_ = IsLightShadowed(light);
        lightBatchQueue->isPerVertex_ = light->GetPerVertex();

        // Update shadow splits
        if (lightBatchQueue->isShadowed_)
        {
            SetupShadowSplits(renderer_, lightBatchQueue, cullCamera_);
        }

        // Process light
        if (lightType == LIGHT_DIRECTIONAL)
        {
            const bool isFocused = light->GetShadowFocus().focus_;

            // For directional light:
            // 1. Iterate over visible geometry;
            // 2. Collect lit geometry;
            // 3. Calculate lit drawables volume if focused and shadowed.
            // Reuse frustum query
            geometriesAndLightsQuery_.ScheduleWork(workQueue_,
                [=](SceneGridCellRef& cellRef, unsigned threadIndex)
            {
                SceneGridCellDrawableSoA& cell = *cellRef.data_;
                VisibleLightsPerThreadData& result = lightsData_[threadIndex];
                const unsigned numSplits = lightBatchQueue->numSplits_;

                for (unsigned index = cellRef.beginIndex_; index < cellRef.endIndex_; ++index)
                {
                    // Skip non-geometries or invisible geometries
                    if (!(cell.drawableFlag_[index] & DRAWABLE_GEOMETRY) || !cell.visible_[index])
                        continue;

                    // Filter by light mask
                    // TODO(eugeneko) Optimize it?
                    const unsigned zoneLightMask = GetActualZone(cell.cachedZone_[index])->GetLightMask();
                    if (!(cell.lightMask_[index] & zoneLightMask & lightMask))
                        continue;

                    // TODO(eugeneko) Get grid index w/o pointer picking
                    Drawable* drawable = cell.drawable_[index];
                    const unsigned gridIndex = drawable->GetDrawableIndex().gridIndex_;

                    // Make and push lit geometry
                    LitGeometryDescIdx litGeometry;
                    litGeometry.drawableIndex_ = gridIndex;
                    litGeometry.lightIndex_ = lightIndex;
                    litGeometry.negativeLight_ = lightBatchQueue->negative_;
                    litGeometry.perVertex_ = lightBatchQueue->isPerVertex_;
                    // TODO(eugeneko) Use true sort value
                    litGeometry.sortValue_ = light->GetSortValue();
                    result.litGeometry_.Push(litGeometry);

                    if (isFocused && lightBatchQueue->isShadowed_)
                    {
                        const float geometryMinZ = cell.minmaxZ_[index].x_;
                        const float geometryMaxZ = cell.minmaxZ_[index].y_;
                        for (unsigned splitIndex = 0; splitIndex < numSplits; ++splitIndex)
                        {
                            const float splitNearZ = lightBatchQueue->shadowNearSplits_[splitIndex];
                            const float splitFarZ = lightBatchQueue->shadowFarSplits_[splitIndex];
                            if (geometryMinZ <= splitFarZ && geometryMaxZ >= splitNearZ)
                                result.lightData_[lightIndex].litGeometriesBox_[splitIndex].Merge(cell.boundingBox_[index]);
                        }
                    }
                }
            });
        }
        else if (lightType == LIGHT_POINT)
        {
            const Sphere lightSphere(light->GetNode()->GetWorldPosition(), light->GetRange());
            // For point light:
            // 1. Query cells in sphere;
            // 2. Iterate over objects in cells;
            // 3. Collect lit geometry;
            workQueue_->ScheduleWork([=](unsigned threadIndex)
            {
                sceneGrid->ProcessCellsInSphere(lightSphere,
                    [=](SceneGridCellDrawableSoA& cell, bool isInside)
                {
                    VisibleLightsPerThreadData& result = lightsData_[threadIndex];
                    for (unsigned index = 0; index < cell.size_; ++index)
                    {
                        // Skip drawables outside the sphere
                        if (!isInside && !cell.IsInSphere(index, lightSphere))
                            continue;

                        // Filter by light mask
                        // TODO(eugeneko) Optimize it?
                        const unsigned zoneLightMask = GetActualZone(cell.cachedZone_[index])->GetLightMask();
                        if (!(cell.lightMask_[index] & zoneLightMask & lightMask))
                            continue;

                        // Add visible to lit geometry
                        if (cell.visible_[index])
                        {
                            // TODO(eugeneko) Get grid index w/o pointer picking
                            Drawable* drawable = cell.drawable_[index];
                            const unsigned gridIndex = drawable->GetDrawableIndex().gridIndex_;

                            LitGeometryDescIdx litGeometry;
                            litGeometry.drawableIndex_ = gridIndex;
                            litGeometry.lightIndex_ = lightIndex;
                            litGeometry.negativeLight_ = lightBatchQueue->negative_;
                            litGeometry.perVertex_ = lightBatchQueue->isPerVertex_;
                            // TODO(eugeneko) Use true sort value
                            litGeometry.sortValue_ = light->GetSortValue();
                            result.litGeometry_.Push(litGeometry);
                        }
                    }
                });
            });
        }
        else if (lightType == LIGHT_SPOT)
        {
            const Frustum lightFrustum = light->GetFrustum();
            const Sphere lightSphere(light->GetNode()->GetWorldPosition(), light->GetRange());
            // For spot lights:
            // 1. Query cells in frustum;
            // 2. Iterate over objects in cells;
            // 3. Collect lit geometry;
            workQueue_->ScheduleWork([=](unsigned threadIndex)
            {
                sceneGrid->ProcessCellsInFrustum(lightFrustum,
                    [=](SceneGridCellDrawableSoA& cell, bool isInside)
                {
                    VisibleLightsPerThreadData& result = lightsData_[threadIndex];
                    for (unsigned index = 0; index < cell.size_; ++index)
                    {
                        if (!isInside && !cell.IsInSphere(index, lightSphere))
                            continue;

                        // Skip drawables outside the sphere
                        if (!isInside && !cell.IsInFrustum(index, lightFrustum))
                            continue;

                        // Filter by light mask
                        // TODO(eugeneko) Optimize it?
                        const unsigned zoneLightMask = GetActualZone(cell.cachedZone_[index])->GetLightMask();
                        if (!(cell.lightMask_[index] & zoneLightMask & lightMask))
                            continue;

                        // Add visible to lit geometry
                        if (cell.visible_[index])
                        {
                            // TODO(eugeneko) Get grid index w/o pointer picking
                            Drawable* drawable = cell.drawable_[index];
                            const unsigned gridIndex = drawable->GetDrawableIndex().gridIndex_;

                            LitGeometryDescIdx litGeometry;
                            litGeometry.drawableIndex_ = gridIndex;
                            litGeometry.lightIndex_ = lightIndex;
                            litGeometry.negativeLight_ = lightBatchQueue->negative_;
                            litGeometry.perVertex_ = lightBatchQueue->isPerVertex_;
                            // TODO(eugeneko) Use true sort value
                            litGeometry.sortValue_ = light->GetSortValue();
                            result.litGeometry_.Push(litGeometry);
                        }
                    }
                });
            });
        }

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

    // Step 1+. Complete threaded job and merge results.
    workQueue_->Complete(M_MAX_UNSIGNED);
    AppendVectorToFirst(lightsData_);

    // Step 2. Update shadow cameras.
    const float sceneMinZ = geometriesAndLights_[0].minZ_;
    const float sceneMaxZ = geometriesAndLights_[0].maxZ_;
    for (unsigned lightIndex = 0; lightIndex < numLights; ++lightIndex)
    {
        LightBatchQueueEx* lightBatchQueue = lightBatchQueues_[lightIndex];
        if (lightBatchQueue->numSplits_ == 0)
            continue;

        const LightPerThreadData& lightData = lightsData_[0].lightData_[lightIndex];
        SetupShadowCameras(lightBatchQueue, lightData, cullCamera_, sceneMinZ, sceneMaxZ);
    }

    // Step 3. Collect shadow casters
    for (unsigned lightIndex = 0; lightIndex < numLights; ++lightIndex)
    {
        LightBatchQueueEx* lightBatchQueue = lightBatchQueues_[lightIndex];
        if (lightBatchQueue->numSplits_ == 0)
            continue;

        Light* light = lightBatchQueue->light_;
        for (unsigned splitIndex = 0; splitIndex < lightBatchQueue->numSplits_; ++splitIndex)
        {
            Camera* shadowCamera = lightBatchQueue->shadowCameras_[splitIndex];

            // TODO(eugeneko) Extract this code
            // @{
            const unsigned lightMask = light->GetLightMask();
            const Frustum& shadowCameraFrustum = shadowCamera->GetFrustum();
            const Matrix3x4& lightView = shadowCamera->GetView();
            const Matrix4& lightProj = shadowCamera->GetProjection();
            const LightType type = light->GetLightType();

            lightBatchQueue->shadowCasterBox_[splitIndex].Clear();

            // Transform scene frustum into shadow camera's view space for shadow caster visibility check. For point & spot lights,
            // we can use the whole scene frustum. For directional lights, use the intersection of the scene frustum and the split
            // frustum, so that shadow casters do not get rendered into unnecessary splits
            Frustum lightViewFrustum;
            if (type != LIGHT_DIRECTIONAL)
                lightViewFrustum = cullCamera_->GetSplitFrustum(sceneMinZ, sceneMaxZ).Transformed(lightView);
            else
                lightViewFrustum = cullCamera_->GetSplitFrustum(Max(sceneMinZ, lightBatchQueue->shadowNearSplits_[splitIndex]),
                    Min(sceneMaxZ, lightBatchQueue->shadowFarSplits_[splitIndex])).Transformed(lightView);

            BoundingBox lightViewFrustumBox(lightViewFrustum);

            // Check for degenerate split frustum: in that case there is no need to get shadow casters
            if (lightViewFrustum.vertices_[0] == lightViewFrustum.vertices_[4])
                continue;
            // @}

            if (type == LIGHT_DIRECTIONAL)
            {
                workQueue_->ScheduleWork([=](unsigned threadIndex)
                {
                    Vector<Drawable*>& shadowCasters = lightBatchQueue->shadowCasters_[splitIndex];
                    BoundingBox& shadowCasterBox = lightBatchQueue->shadowCasterBox_[splitIndex];
                    assert(shadowCasters.Empty());

                    sceneGrid->ProcessCellsInFrustum(shadowCameraFrustum,
                        [=, &shadowCasters, &shadowCasterBox](SceneGridCellDrawableSoA& cell, bool isInside)
                    {
                        for (unsigned index = 0; index < cell.size_; ++index)
                        {
                            // Skip if view mask doesn't match or if not shadow caster
                            if (!cell.castShadows_[index] || !(cell.viewMask_[index] & viewMask_))
                                continue;

                            // Skip drawables outside the frustum
                            if (!isInside && !cell.IsInFrustum(index, shadowCameraFrustum))
                                continue;

                            // Skip if shadow mask doesn't match
                            const unsigned zoneShadowMask = GetActualZone(cell.cachedZone_[index])->GetShadowMask();
                            if (!(cell.shadowMask_[index] & zoneShadowMask & lightMask))
                                continue;

                            // Calculate draw distance
                            const float drawDistance = cell.drawDistance_[index];
                            const float shadowDistance = GetEffectiveShadowDistance(cell.shadowDistance_[index], drawDistance);
                            const float distance = cullCamera_->GetDistance(cell.boundingSphere_[index].center_);

                            // Discard if drawable is too far
                            if (shadowDistance > 0.0f && distance > shadowDistance)
                                continue;

                            // Project shadow caster bounding box to light view space for visibility check
                            const BoundingBox lightViewBox = cell.boundingBox_[index].Transformed(lightView);

                            if (IsShadowCasterVisible(cell.visible_[index], lightViewBox, shadowCamera, lightView, lightViewFrustum, lightViewFrustumBox))
                            {
                                // Merge to shadow caster bounding box (only needed for focused spot lights) and add to the list
                                if (type == LIGHT_SPOT && light->GetShadowFocus().focus_)
                                {
                                    const BoundingBox lightProjBox{ lightViewBox.Projected(lightProj) };
                                    shadowCasterBox.Merge(lightProjBox);
                                }
                                shadowCasters.Push(cell.drawable_[index]);
                            }
                        }
                    });
                });
            }
            else
            {
                assert(0);
            }
        }
    }
    workQueue_->Complete(M_MAX_UNSIGNED);
}

void BatchCollector::SortLitGeometries(SceneGrid* sceneGrid)
{
    SceneGridDrawableSoA& globalData = sceneGrid->GetGlobalData();
    Vector<LitGeometryDescIdx>& unsortedLitGeometry = lightsData_[0].litGeometry_;

    // 1st pass: calculate number of lights for each geometry
    for (const LitGeometryDescIdx& lit : unsortedLitGeometry)
        ++globalData.numLights_[lit.drawableIndex_];

    // 2nd pass: calculate ranges in destination array
    unsigned firstLight = 0;
    for (unsigned i = 0; i < globalData.size_; ++i)
    {
        globalData.firstLight_[i] = firstLight;
        firstLight += globalData.numLights_[i];
    }

    // 3rd pass: reset number of lights to reuse it for filling array
    for (unsigned i = 0; i < globalData.size_; ++i)
        globalData.numLights_[i] = 0;

    // 4th pass: fill destination array
    sortedLitGeometries_.Resize(unsortedLitGeometry.Size());
    for (const LitGeometryDescIdx& lit : unsortedLitGeometry)
    {
        const unsigned drawableIndex = lit.drawableIndex_;
        const unsigned litIndex = globalData.firstLight_[drawableIndex] + globalData.numLights_[drawableIndex];
        sortedLitGeometries_[litIndex] = LitGeometryDescPacked(lit);
        ++globalData.numLights_[drawableIndex];
    }

    // 5th pass: partially sort lit geometries
    for (unsigned i = 0; i < globalData.size_; ++i)
    {
        if (globalData.visible_[i])
        {
            const unsigned numLights = globalData.numLights_[i];
            numLightsPerVisibleGeometry_.Push(numLights);
            if (numLights > 1)
            {
                auto begin = sortedLitGeometries_.Begin() + globalData.firstLight_[i];
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

void BatchCollector::CollectShadowBatches(const IntVector2& viewSize, int materialQuality)
{
    const unsigned numLights = GetVisibleLights().Size();

    // Get shadow maps and finalize shadow cameras
    for (unsigned lightIndex = 0; lightIndex < numLights; ++lightIndex)
    {
        LightBatchQueueEx* lightBatchQueue = lightBatchQueues_[lightIndex];
        Light* light = lightBatchQueue->light_;
        if (lightBatchQueue->numSplits_ > 0)
        {
            lightBatchQueue->shadowMap_ = renderer_->GetShadowMap(
                light, cullCamera_, (unsigned)viewSize.x_, (unsigned)viewSize.y_);

            if (!lightBatchQueue->shadowMap_)
                lightBatchQueue->numSplits_ = 0;
        }

        for (unsigned splitIndex = 0; splitIndex < lightBatchQueue->numSplits_; ++splitIndex)
        {
            const IntRect shadowViewport = GetShadowMapViewport(light, splitIndex, lightBatchQueue->shadowMap_);
            lightBatchQueue->shadowViewports_[splitIndex] = shadowViewport;

            Camera* shadowCamera = lightBatchQueue->shadowCameras_[splitIndex];
            const BoundingBox& shaderCasterBox = lightBatchQueue->shadowCasterBox_[splitIndex];
            FinalizeShadowCamera(shadowCamera, light, shadowViewport, shaderCasterBox);
        }
    }

    for (unsigned lightIndex = 0; lightIndex < numLights; ++lightIndex)
    {
        LightBatchQueueEx* lightBatchQueue = lightBatchQueues_[lightIndex];
        const unsigned numSplits = lightBatchQueue->numSplits_;
        if (numSplits == 0)
            continue;

        // TODO(eugeneko) Remove this resize
        lightBatchQueue->shadowSplits_.Resize(numSplits);
        lightBatchQueue->shadowSplitsData_.Resize(numSplits);
        for (unsigned splitIndex = 0; splitIndex < numSplits; ++splitIndex)
        {
            workQueue_->ScheduleWork([=](unsigned threadIndex)
            {
                const Vector<Drawable*>& shadowCasters = lightBatchQueue->shadowCasters_[splitIndex];
                BatchQueueData& shadowQueueData = lightBatchQueue->shadowSplitsData_[splitIndex];
                ShadowBatchQueue& shadowQueue = lightBatchQueue->shadowSplits_[splitIndex];

                shadowQueue.shadowViewport_ = lightBatchQueue->shadowViewports_[splitIndex];
                shadowQueue.shadowCamera_ = lightBatchQueue->shadowCameras_[splitIndex];
                shadowQueue.nearSplit_ = lightBatchQueue->shadowNearSplits_[splitIndex];
                shadowQueue.farSplit_ = lightBatchQueue->shadowFarSplits_[splitIndex];

                shadowQueueData.ShallowClear();
                for (Drawable* drawable : shadowCasters)
                {
                    const Vector<SourceBatch>& batches = drawable->GetBatches();

                    for (unsigned i = 0; i < batches.Size(); ++i)
                    {
                        const SourceBatch& srcBatch = batches[i];

                        Technique* tech = GetTechnique(materialQuality, drawable->GetLodDistance(), srcBatch.material_, renderer_->GetDefaultMaterial());
                        if (!srcBatch.geometry_ || !srcBatch.numWorldTransforms_ || !tech)
                            continue;

                        Pass* pass = tech->GetSupportedPass(Technique::shadowPassIndex);
                        // Skip if material has no shadow pass
                        if (!pass)
                            continue;

                        Batch destBatch(srcBatch);
                        destBatch.pass_ = pass;
                        destBatch.zone_ = nullptr;

                        shadowQueueData.AddBatch(destBatch, true);
                    }
                }

                // Fill batches
                shadowQueue.shadowBatches_.Clear(maxSortedInstances_);
                shadowQueueData.ExportBatches(shadowQueue.shadowBatches_);
                shadowQueueData.ExportBatchGroups(shadowQueue.shadowBatches_);
            });
        }
    }

    workQueue_->Complete(M_MAX_UNSIGNED);
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
    for (unsigned lightIndex = 0; lightIndex < GetVisibleLights().Size(); ++lightIndex)
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
    const unsigned numLights = GetVisibleLights().Size();
    for (unsigned lightIndex = 0; lightIndex < numLights; ++lightIndex)
    {
        Light* light = GetVisibleLight(lightIndex);
        if (light->GetPerVertex())
            continue;

        LightBatchQueueEx* lightBatchQueue = lightBatchQueues_[lightIndex];
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
        Light* light = GetVisibleLight(lightIndex);
        if (light->GetPerVertex())
            continue;

        LightBatchQueueEx* lightBatchQueue = lightBatchQueues_[lightIndex];
        assert(lightBatchQueue);
        FinalizeBatchQueue(lightBatchQueue->litBatches_, true);
        FinalizeBatchQueue(lightBatchQueue->litBaseBatches_, true);

        for (unsigned splitIndex = 0; splitIndex < lightBatchQueue->numSplits_; ++splitIndex)
        {
            FinalizeBatchQueue(lightBatchQueue->shadowSplits_[splitIndex].shadowBatches_, true);
        }
    }
}

Zone* BatchCollector::GetDrawableZone(Drawable* drawable) const
{
    if (zonesData_.cameraZoneOverride_)
        return zonesData_.cameraZone_;
    Zone* drawableZone = drawable->GetZone();
    return drawableZone ? drawableZone : zonesData_.cameraZone_;
}

Zone* BatchCollector::GetActualZone(Zone* drawableZone) const
{
    if (zonesData_.cameraZoneOverride_)
        return zonesData_.cameraZone_;
    return drawableZone ? drawableZone : zonesData_.cameraZone_;
}

LightBatchQueueEx* BatchCollector::GetOrCreateLightBatchQueue(Light* light)
{
    auto iterLightBatchQueue = lightBatchQueueMap_.Find(light);
    if (iterLightBatchQueue == lightBatchQueueMap_.End())
    {
        LightBatchQueueEx* queue = AllocateLightBatchQueue();
        queue->light_ = light;
        iterLightBatchQueue = lightBatchQueueMap_.Insert(MakePair(light, queue));
    }
    return iterLightBatchQueue->second_;
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

bool BatchCollector::IsLightShadowed(Light* light)
{
    // Check if light should be shadowed
    bool isShadowed = light->GetCastShadows() && !light->GetPerVertex() && light->GetShadowIntensity() < 1.0f;
    // If shadow distance non-zero, check it
    if (isShadowed && light->GetShadowDistance() > 0.0f && light->GetDistance() > light->GetShadowDistance())
        isShadowed = false;
    // OpenGL ES can not support point light shadows
#ifdef GL_ES_VERSION_2_0
    if (isShadowed && light->GetLightType() == LIGHT_POINT)
        isShadowed = false;
#endif
    return isShadowed;
}

void BatchCollector::SetupShadowSplits(Renderer* renderer, LightBatchQueueEx* queue, Camera* cullCamera)
{
    Light* light = queue->light_;

    unsigned splits = 0;

    switch (light->GetLightType())
    {
    case LIGHT_DIRECTIONAL:
        {
            const CascadeParameters& cascade = light->GetShadowCascade();

            float nearSplit = cullCamera->GetNearClip();
            float farSplit = 0;
            int numSplits = light->GetNumShadowSplits();

            while (splits < numSplits)
            {
                // If split is completely beyond camera far clip, we are done
                if (nearSplit > cullCamera->GetFarClip())
                    break;

                farSplit = Min(cullCamera->GetFarClip(), cascade.splits_[splits]);
                if (farSplit <= nearSplit)
                    break;

                // Setup the shadow camera for the split
                queue->shadowCameras_[splits] = renderer->GetShadowCamera();
                queue->shadowNearSplits_[splits] = nearSplit;
                queue->shadowFarSplits_[splits] = farSplit;

                nearSplit = farSplit;
                ++splits;
            }
        }
        break;

    case LIGHT_SPOT:
        {
            queue->shadowCameras_[0] = renderer->GetShadowCamera();
            splits = 1;
        }
        break;

    case LIGHT_POINT:
        {
            static const Vector3* directions[] =
            {
                &Vector3::RIGHT,
                &Vector3::LEFT,
                &Vector3::UP,
                &Vector3::DOWN,
                &Vector3::FORWARD,
                &Vector3::BACK
            };

            for (unsigned i = 0; i < MAX_CUBEMAP_FACES; ++i)
            {
                queue->shadowCameras_[i] = renderer->GetShadowCamera();
            }

            splits = MAX_CUBEMAP_FACES;
        }
        break;
    }

    queue->numSplits_ = splits;
}

void BatchCollector::SetupShadowCameras(LightBatchQueueEx* queue, const LightPerThreadData& lightData,
    Camera* cullCamera, float minZ, float maxZ)
{
    Light* light = queue->light_;
    const unsigned numSplits = queue->numSplits_;

    switch (light->GetLightType())
    {
    case LIGHT_DIRECTIONAL:
        {
            for (unsigned splitIndex = 0; splitIndex < numSplits; ++splitIndex)
            {
                Camera* shadowCamera = queue->shadowCameras_[splitIndex];
                const float nearSplit = queue->shadowNearSplits_[splitIndex];
                const float farSplit = queue->shadowFarSplits_[splitIndex];
                const BoundingBox& litGeometryBox = lightData.litGeometriesBox_[splitIndex];
                SetupDirLightShadowCamera(shadowCamera, light, nearSplit, farSplit, litGeometryBox,
                    cullCamera, minZ, maxZ);
            }
        }
        break;

    case LIGHT_SPOT:
        {
            assert(numSplits == 1);

            Camera* shadowCamera = queue->shadowCameras_[0];
            Node* cameraNode = shadowCamera->GetNode();
            Node* lightNode = light->GetNode();

            cameraNode->SetTransform(lightNode->GetWorldPosition(), lightNode->GetWorldRotation());
            shadowCamera->SetNearClip(light->GetShadowNearFarRatio() * light->GetRange());
            shadowCamera->SetFarClip(light->GetRange());
            shadowCamera->SetFov(light->GetFov());
            shadowCamera->SetAspectRatio(light->GetAspectRatio());
        }
        break;

    case LIGHT_POINT:
        {
            assert(numSplits == MAX_CUBEMAP_FACES);

            static const Vector3* directions[] =
            {
                &Vector3::RIGHT,
                &Vector3::LEFT,
                &Vector3::UP,
                &Vector3::DOWN,
                &Vector3::FORWARD,
                &Vector3::BACK
            };

            for (unsigned splitIndex = 0; splitIndex < MAX_CUBEMAP_FACES; ++splitIndex)
            {
                Camera* shadowCamera = queue->shadowCameras_[splitIndex];
                Node* cameraNode = shadowCamera->GetNode();

                // When making a shadowed point light, align the splits along X, Y and Z axes regardless of light rotation
                cameraNode->SetPosition(light->GetNode()->GetWorldPosition());
                cameraNode->SetDirection(*directions[splitIndex]);
                shadowCamera->SetNearClip(light->GetShadowNearFarRatio() * light->GetRange());
                shadowCamera->SetFarClip(light->GetRange());
                shadowCamera->SetFov(90.0f);
                shadowCamera->SetAspectRatio(1.0f);
            }
        }
        break;
    }
}

void BatchCollector::SetupDirLightShadowCamera(Camera* shadowCamera, Light* light,
    float nearSplit, float farSplit, const BoundingBox& litGeometriesBox,
    Camera* cullCamera, float minZ, float maxZ)
{
    Node* shadowCameraNode = shadowCamera->GetNode();
    Node* lightNode = light->GetNode();
    float extrusionDistance = Min(cullCamera->GetFarClip(), light->GetShadowMaxExtrusion());
    const FocusParameters& parameters = light->GetShadowFocus();

    // Calculate initial position & rotation
    Vector3 pos = cullCamera->GetNode()->GetWorldPosition() - extrusionDistance * lightNode->GetWorldDirection();
    shadowCameraNode->SetTransform(pos, lightNode->GetWorldRotation());

    // Calculate main camera shadowed frustum in light's view space
    farSplit = Min(farSplit, cullCamera->GetFarClip());
    // Use the scene Z bounds to limit frustum size if applicable
    if (parameters.focus_)
    {
        nearSplit = Max(minZ, nearSplit);
        farSplit = Min(maxZ, farSplit);
    }

    Frustum splitFrustum = cullCamera->GetSplitFrustum(nearSplit, farSplit);
    Polyhedron frustumVolume;
    frustumVolume.Define(splitFrustum);
    // If focusing enabled, clip the frustum volume by the combined bounding box of the lit geometries within the frustum
    if (parameters.focus_ && litGeometriesBox.Defined())
    {
        frustumVolume.Clip(litGeometriesBox);
        // If volume became empty, restore it to avoid zero size
        if (frustumVolume.Empty())
            frustumVolume.Define(splitFrustum);
    }

    // Transform frustum volume to light space
    const Matrix3x4& lightView = shadowCamera->GetView();
    frustumVolume.Transform(lightView);

    // Fit the frustum volume inside a bounding box. If uniform size, use a sphere instead
    BoundingBox shadowBox;
    if (!parameters.nonUniform_)
        shadowBox.Define(Sphere(frustumVolume));
    else
        shadowBox.Define(frustumVolume);

    shadowCamera->SetOrthographic(true);
    shadowCamera->SetAspectRatio(1.0f);
    shadowCamera->SetNearClip(0.0f);
    shadowCamera->SetFarClip(shadowBox.max_.z_);

    // Center shadow camera on the bounding box. Can not snap to texels yet as the shadow map viewport is unknown
    QuantizeDirLightShadowCamera(shadowCamera, light, IntRect(0, 0, 0, 0), shadowBox);
}

void BatchCollector::FinalizeShadowCamera(Camera* shadowCamera, Light* light, const IntRect& shadowViewport,
    const BoundingBox& shadowCasterBox)
{
    const FocusParameters& parameters = light->GetShadowFocus();
    auto shadowMapWidth = (float)(shadowViewport.Width());
    LightType type = light->GetLightType();

    if (type == LIGHT_DIRECTIONAL)
    {
        BoundingBox shadowBox;
        shadowBox.max_.y_ = shadowCamera->GetOrthoSize() * 0.5f;
        shadowBox.max_.x_ = shadowCamera->GetAspectRatio() * shadowBox.max_.y_;
        shadowBox.min_.y_ = -shadowBox.max_.y_;
        shadowBox.min_.x_ = -shadowBox.max_.x_;

        // Requantize and snap to shadow map texels
        QuantizeDirLightShadowCamera(shadowCamera, light, shadowViewport, shadowBox);
    }

    if (type == LIGHT_SPOT && parameters.focus_)
    {
        float viewSizeX = Max(Abs(shadowCasterBox.min_.x_), Abs(shadowCasterBox.max_.x_));
        float viewSizeY = Max(Abs(shadowCasterBox.min_.y_), Abs(shadowCasterBox.max_.y_));
        float viewSize = Max(viewSizeX, viewSizeY);
        // Scale the quantization parameters, because view size is in projection space (-1.0 - 1.0)
        float invOrthoSize = 1.0f / shadowCamera->GetOrthoSize();
        float quantize = parameters.quantize_ * invOrthoSize;
        float minView = parameters.minView_ * invOrthoSize;

        viewSize = Max(ceilf(viewSize / quantize) * quantize, minView);
        if (viewSize < 1.0f)
            shadowCamera->SetZoom(1.0f / viewSize);
    }

    // Perform a finalization step for all lights: ensure zoom out of 2 pixels to eliminate border filtering issues
    // For point lights use 4 pixels, as they must not cross sides of the virtual cube map (maximum 3x3 PCF)
    if (shadowCamera->GetZoom() >= 1.0f)
    {
        if (light->GetLightType() != LIGHT_POINT)
            shadowCamera->SetZoom(shadowCamera->GetZoom() * ((shadowMapWidth - 2.0f) / shadowMapWidth));
        else
        {
#ifdef URHO3D_OPENGL
            shadowCamera->SetZoom(shadowCamera->GetZoom() * ((shadowMapWidth - 3.0f) / shadowMapWidth));
#else
            shadowCamera->SetZoom(shadowCamera->GetZoom() * ((shadowMapWidth - 4.0f) / shadowMapWidth));
#endif
        }
    }
}

void BatchCollector::QuantizeDirLightShadowCamera(Camera* shadowCamera, Light* light, const IntRect& shadowViewport,
    const BoundingBox& viewBox)
{
    Node* shadowCameraNode = shadowCamera->GetNode();
    const FocusParameters& parameters = light->GetShadowFocus();
    auto shadowMapWidth = (float)(shadowViewport.Width());

    float minX = viewBox.min_.x_;
    float minY = viewBox.min_.y_;
    float maxX = viewBox.max_.x_;
    float maxY = viewBox.max_.y_;

    Vector2 center((minX + maxX) * 0.5f, (minY + maxY) * 0.5f);
    Vector2 viewSize(maxX - minX, maxY - minY);

    // Quantize size to reduce swimming
    // Note: if size is uniform and there is no focusing, quantization is unnecessary
    if (parameters.nonUniform_)
    {
        viewSize.x_ = ceilf(sqrtf(viewSize.x_ / parameters.quantize_));
        viewSize.y_ = ceilf(sqrtf(viewSize.y_ / parameters.quantize_));
        viewSize.x_ = Max(viewSize.x_ * viewSize.x_ * parameters.quantize_, parameters.minView_);
        viewSize.y_ = Max(viewSize.y_ * viewSize.y_ * parameters.quantize_, parameters.minView_);
    }
    else if (parameters.focus_)
    {
        viewSize.x_ = Max(viewSize.x_, viewSize.y_);
        viewSize.x_ = ceilf(sqrtf(viewSize.x_ / parameters.quantize_));
        viewSize.x_ = Max(viewSize.x_ * viewSize.x_ * parameters.quantize_, parameters.minView_);
        viewSize.y_ = viewSize.x_;
    }

    shadowCamera->SetOrthoSize(viewSize);

    // Center shadow camera to the view space bounding box
    Quaternion rot(shadowCameraNode->GetWorldRotation());
    Vector3 adjust(center.x_, center.y_, 0.0f);
    shadowCameraNode->Translate(rot * adjust, TS_WORLD);

    // If the shadow map viewport is known, snap to whole texels
    if (shadowMapWidth > 0.0f)
    {
        Vector3 viewPos(rot.Inverse() * shadowCameraNode->GetWorldPosition());
        // Take into account that shadow map border will not be used
        float invActualSize = 1.0f / (shadowMapWidth - 2.0f);
        Vector2 texelSize(viewSize.x_ * invActualSize, viewSize.y_ * invActualSize);
        Vector3 snap(-fmodf(viewPos.x_, texelSize.x_), -fmodf(viewPos.y_, texelSize.y_), 0.0f);
        shadowCameraNode->Translate(rot * snap, TS_WORLD);
    }
}

bool BatchCollector::IsShadowCasterVisible(bool drawableVisible, BoundingBox lightViewBox, Camera* shadowCamera, const Matrix3x4& lightView, const Frustum& lightViewFrustum, const BoundingBox& lightViewFrustumBox)
{
    if (shadowCamera->IsOrthographic())
    {
        // Extrude the light space bounding box up to the far edge of the frustum's light space bounding box
        lightViewBox.max_.z_ = Max(lightViewBox.max_.z_, lightViewFrustumBox.max_.z_);
        return lightViewFrustum.IsInsideFast(lightViewBox) != OUTSIDE;
    }
    else
    {
        // If light is not directional, can do a simple check: if object is visible, its shadow is too
        if (drawableVisible)
            return true;

        // For perspective lights, extrusion direction depends on the position of the shadow caster
        Vector3 center = lightViewBox.Center();
        Ray extrusionRay(center, center);

        float extrusionDistance = shadowCamera->GetFarClip();
        float originalDistance = Clamp(center.Length(), M_EPSILON, extrusionDistance);

        // Because of the perspective, the bounding box must also grow when it is extruded to the distance
        float sizeFactor = extrusionDistance / originalDistance;

        // Calculate the endpoint box and merge it to the original. Because it's axis-aligned, it will be larger
        // than necessary, so the test will be conservative
        Vector3 newCenter = extrusionDistance * extrusionRay.direction_;
        Vector3 newHalfSize = lightViewBox.Size() * sizeFactor * 0.5f;
        BoundingBox extrudedBox(newCenter - newHalfSize, newCenter + newHalfSize);
        lightViewBox.Merge(extrudedBox);

        return lightViewFrustum.IsInsideFast(lightViewBox) != OUTSIDE;
    }
}

IntRect BatchCollector::GetShadowMapViewport(Light* light, int splitIndex, Texture2D* shadowMap)
{
    const int width = shadowMap->GetWidth();
    const int height = shadowMap->GetHeight();

    switch (light->GetLightType())
    {
    case LIGHT_DIRECTIONAL:
        {
            int numSplits = light->GetNumShadowSplits();
            if (numSplits == 1)
                return {0, 0, width, height};
            else if (numSplits == 2)
                return {splitIndex * width / 2, 0, (splitIndex + 1) * width / 2, height};
            else
                return {(splitIndex & 1) * width / 2, (splitIndex / 2) * height / 2,
                    ((splitIndex & 1) + 1) * width / 2, (splitIndex / 2 + 1) * height / 2};
        }

    case LIGHT_SPOT:
        return {0, 0, width, height};

    case LIGHT_POINT:
        return {(splitIndex & 1) * width / 2, (splitIndex / 2) * height / 3,
            ((splitIndex & 1) + 1) * width / 2, (splitIndex / 2 + 1) * height / 3};
    }

    return {};
}

Technique* BatchCollector::GetTechnique(int materialQuality, float lodDistance, Material* material, Material* defaultMaterial)
{
    if (!material)
        return defaultMaterial->GetTechniques()[0].technique_;

    const Vector<TechniqueEntry>& techniques = material->GetTechniques();
    // If only one technique, no choice
    if (techniques.Size() == 1)
        return techniques[0].technique_;
    else
    {
        // Check for suitable technique. Techniques should be ordered like this:
        // Most distant & highest quality
        // Most distant & lowest quality
        // Second most distant & highest quality
        // ...
        for (unsigned i = 0; i < techniques.Size(); ++i)
        {
            const TechniqueEntry& entry = techniques[i];
            Technique* tech = entry.technique_;

            if (!tech || (!tech->IsSupported()) || materialQuality < entry.qualityLevel_)
                continue;
            if (lodDistance >= entry.lodDistance_)
                return tech;
        }

        // If no suitable technique found, fallback to the last
        return techniques.Size() ? techniques.Back().technique_ : nullptr;
    }
}

}
