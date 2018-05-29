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

#include "../Core/Profiler.h"
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

using DrawableVector = Vector<Drawable*>;
using LightVector = Vector<Light*>;
using ZoneVector = Vector<Zone*>;

struct SceneGridCellDrawableSoA
{
    unsigned size_ = 0;

    DrawableVector drawable_;
    Vector<StringHash> drawableType_;
    Vector<unsigned> drawableFlag_;

    Vector<bool> dirty_;
    Vector<Matrix3x4> transform_;
    Vector<BoundingBox> boundingBox_;
    Vector<Sphere> boundingSphere_;

    Vector<float> drawDistance_;
    Vector<float> shadowDistance_;
    Vector<unsigned> viewMask_;
    Vector<unsigned> lightMask_;
    Vector<unsigned> shadowMask_;
    Vector<unsigned> zoneMask_;
    Vector<bool> occluder_;
    Vector<bool> occludee_;
    Vector<bool> castShadows_;

    ZoneVector cachedZone_;
    Vector<bool> cachedZoneDirty_;
    // TODO(eugeneko) Worth removing
    Vector<unsigned> cachedZoneViewMask_;
    Vector<unsigned> cachedZoneLightMask_;
    Vector<unsigned> cachedZoneShadowMask_;

    Vector<bool> visible_;
    Vector<Vector2> minmaxZ_;

    bool IsValid() const
    {
        return size_ == drawable_.Size()
            && size_ == drawableType_.Size()
            && size_ == drawableFlag_.Size()

            && size_ == dirty_.Size()
            && size_ == transform_.Size()
            && size_ == boundingBox_.Size()
            && size_ == boundingSphere_.Size()

            && size_ == drawDistance_.Size()
            && size_ == shadowDistance_.Size()
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
            && size_ == cachedZoneShadowMask_.Size()

            && size_ == visible_.Size()
            && size_ == minmaxZ_.Size()
            ;
    }
    void Push(Drawable* drawable)
    {
        drawable_.Push(drawable);
        drawableType_.Push(drawable->GetType());
        drawableFlag_.Push(drawable->GetDrawableFlags());

        dirty_.Push(true);
        transform_.Push(Matrix3x4());
        boundingBox_.Push(BoundingBox());
        boundingSphere_.Push(Sphere());

        drawDistance_.Push(0.0f);
        shadowDistance_.Push(0.0f);
        viewMask_.Push(0);
        lightMask_.Push(0);
        shadowMask_.Push(0);
        zoneMask_.Push(0);
        occluder_.Push(false);
        occludee_.Push(false);
        castShadows_.Push(false);
        UpdateDrawable(size_, drawable);

        cachedZone_.Push(nullptr);
        cachedZoneDirty_.Push(true);
        cachedZoneViewMask_.Push(0);
        cachedZoneLightMask_.Push(0);
        cachedZoneShadowMask_.Push(0);

        visible_.Push(false);
        minmaxZ_.Push(Vector2::ZERO);

        ++size_;
        assert(IsValid());
    }
    void EraseSwap(unsigned index)
    {
        assert(index < size_);
        --size_;

        if (size_ >= 1)
        {
            // Re-assign last drawable index
            drawable_.Back()->SetDrawableIndexInCell(index);
        }

        drawable_.EraseSwap(index);
        drawableType_.EraseSwap(index);
        drawableFlag_.EraseSwap(index);

        dirty_.EraseSwap(index);
        transform_.EraseSwap(index);
        boundingBox_.EraseSwap(index);
        boundingSphere_.EraseSwap(index);

        drawDistance_.EraseSwap(index);
        shadowDistance_.EraseSwap(index);
        viewMask_.EraseSwap(index);
        lightMask_.EraseSwap(index);
        shadowMask_.EraseSwap(index);
        zoneMask_.EraseSwap(index);
        occluder_.EraseSwap(index);
        occludee_.EraseSwap(index);
        castShadows_.EraseSwap(index);

        cachedZone_.EraseSwap(index);
        cachedZoneDirty_.EraseSwap(index);
        cachedZoneViewMask_.EraseSwap(index);
        cachedZoneLightMask_.EraseSwap(index);
        cachedZoneShadowMask_.EraseSwap(index);

        visible_.EraseSwap(index);
        minmaxZ_.EraseSwap(index);

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
        boundingBox_.Clear();
        boundingSphere_.Clear();

        drawDistance_.Clear();
        shadowDistance_.Clear();
        viewMask_.Clear();
        lightMask_.Clear();
        shadowMask_.Clear();
        zoneMask_.Clear();
        occluder_.Clear();
        occludee_.Clear();
        castShadows_.Clear();

        cachedZone_.Clear();
        cachedZoneDirty_.Clear();
        cachedZoneViewMask_.Clear();
        cachedZoneLightMask_.Clear();
        cachedZoneShadowMask_.Clear();

        visible_.Clear();
        minmaxZ_.Clear();

        assert(IsValid());
    }

    /// Clear temporary data.
    void ClearTemporary()
    {
        for (unsigned i = 0; i < size_; ++i)
        {
            visible_[i] = false;
        }
    }
    /// Update drawable data.
    void UpdateDrawable(unsigned index, Drawable* drawable)
    {
        drawDistance_[index] = drawable->GetDrawDistance();
        shadowDistance_[index] = drawable->GetShadowDistance();
        viewMask_[index] = drawable->GetViewMask();
        lightMask_[index] = drawable->GetLightMask();
        shadowMask_[index] = drawable->GetShadowMask();
        zoneMask_[index] = drawable->GetZoneMask();
        occluder_[index] = drawable->IsOccluder();
        occludee_[index] = drawable->IsOccludee();
        castShadows_[index] = drawable->GetCastShadows();
    }
    /// Reset specified zone for each drawable.
    void ResetZone(Zone* zone)
    {
        for (unsigned i = 0; i < size_; ++i)
        {
            if (cachedZone_[i] == zone)
            {
                cachedZone_[i] = nullptr;
                cachedZoneDirty_[i] = true;
            }
        }
    }
    /// Mark zones dirty.
    void MarkZonesDirty()
    {
        for (unsigned i = 0; i < size_; ++i)
            cachedZoneDirty_[i] = true;
    }

    bool MatchViewMask(unsigned index, unsigned viewMask) const { return !!(viewMask_[index] & viewMask); }
    bool IsZone(unsigned index) const { return !!(drawableFlag_[index] & DRAWABLE_ZONE); }
    bool IsLight(unsigned index) const { return !!(drawableFlag_[index] & DRAWABLE_LIGHT); }
    bool IsGeometry(unsigned index) const { return !!(drawableFlag_[index] & DRAWABLE_GEOMETRY); }
    Intersection IntersectSphereFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInside(boundingSphere_[index]); }
    Intersection IntersectBoxFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInsideFast(boundingBox_[index]); }
    Intersection IntersectSphereSphere(unsigned index, const Sphere& sphere) const { return sphere.IsInside(boundingSphere_[index]); }
    Intersection IntersectBoxSphere(unsigned index, const Sphere& sphere) const { return sphere.IsInsideFast(boundingBox_[index]); }
    bool IsInFrustum(unsigned index, const Frustum& frustum) const
    {
        const Intersection sphereIntersection = IntersectSphereFrustum(index, frustum);
        return sphereIntersection != OUTSIDE
            && (sphereIntersection == INSIDE || IntersectBoxFrustum(index, frustum) != OUTSIDE);
    }
    bool IsInSphere(unsigned index, const Sphere& sphere) const
    {
        const Intersection sphereIntersection = IntersectSphereSphere(index, sphere);
        return sphereIntersection != OUTSIDE
            && (sphereIntersection == INSIDE || IntersectBoxSphere(index, sphere) != OUTSIDE);
    }
    bool IsCenterInFrustum(unsigned index, const Frustum& frustum) const { return frustum.IsInside(boundingSphere_[index].center_) != OUTSIDE; }
    bool IsOccludedByBuffer(unsigned index, OcclusionBuffer* buffer) const
    {
        return buffer && occludee_[index] && !buffer->IsVisible(boundingBox_[index]);

    }

};

struct SceneGridCell
{
    IntVector3 index_;
    BoundingBox innerBoundingBox_;
    BoundingBox outerBoundingBox_;
    SceneGridCellDrawableSoA data_;
    bool IsInside(const BoundingBox& boundingBox) const
    {
        return innerBoundingBox_.IsInside(boundingBox.Center()) == INSIDE
            && outerBoundingBox_.IsInside(boundingBox) == INSIDE;
    }
};

struct SceneGridDrawableSoA
{
    unsigned size_{};
    DrawableVector drawable_;
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
    // TODO(eugeneko) remove me
    void ClearVisible()
    {
        ClearTemporary();
    }
    void ClearTemporary()
    {
        for (unsigned i = 0; i < size_; ++i)
        {
            visible_[i] = false;
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
        assert(index < size_);
        --size_;

        if (size_ >= 1)
        {
            // Re-assign last drawable index
            drawable_.Back()->SetDrawableIndexInGrid(index);
        }

        drawable_.EraseSwap(index);
        drawableType_.EraseSwap(index);
        visible_.EraseSwap(index);
        numLights_.EraseSwap(index);
        firstLight_.EraseSwap(index);

        assert(IsValid());
    }
    void Clear()
    {
        size_ = 0;

        drawable_.Clear();
        drawableType_.Clear();
        visible_.Clear();
        numLights_.Clear();
        firstLight_.Clear();

        assert(IsValid());
    }

};

struct SceneGridCellRef
{
    Intersection intersection_{};
    unsigned beginIndex_{};
    unsigned endIndex_{};
    SceneGridCellDrawableSoA* data_{};
};

struct SceneGridQueryResult
{
    void Clear()
    {
        cellMasks_.Clear();
        cells_.Clear();
        threadRanges_.Clear();
    }
    template <class T> void ScheduleWork(WorkQueue* workQueue, const T& work) const
    {
        const unsigned numThreads = threadRanges_.Size();
        for (unsigned i = 0; i < numThreads; ++i)
        {
            const unsigned beginCell = threadRanges_[i].first_;
            const unsigned endCell = threadRanges_[i].second_;
            workQueue->ScheduleWork(
                [=](unsigned threadIndex)
            {
                for (unsigned cellIndex = beginCell; cellIndex < endCell; ++cellIndex)
                {
                    work(cells_[cellIndex], threadIndex);
                }
            });
        }
    }
    Vector<Intersection> cellMasks_;
    Vector<SceneGridCellRef> cells_;
    Vector<Pair<unsigned, unsigned>> threadRanges_;
};

class SceneGrid
{
public:
    /// Construct.
    SceneGrid() = default;
    /// Non-copyable.
    SceneGrid(const SceneGrid& other) = delete;
    /// Destruct.
    ~SceneGrid()
    {
        RemoveAllDrawables();
    }
    /// Reset scene grid.
    void Reset(const BoundingBox& boundingBox, const IntVector3& numCellsXYZ, float pad)
    {
        boundingBox_ = boundingBox;
        numCellsXYZ_ = numCellsXYZ;
        numCells_ = static_cast<unsigned>(numCellsXYZ_.x_ * numCellsXYZ_.y_ * numCellsXYZ_.z_) + 1;
        floatNumCellsXYZ_ = Vector3(static_cast<float>(numCellsXYZ_.x_), static_cast<float>(numCellsXYZ_.y_), static_cast<float>(numCellsXYZ_.z_));
        cellSize_ = boundingBox_.Size() / floatNumCellsXYZ_;

        drawablesData_.Clear();
        cells_.Resize(numCells_);

        // Reset default cell
        SceneGridCell* defaultCell = GetDefaultCell();
        defaultCell->index_ = IntVector3::ONE * -1;
        defaultCell->innerBoundingBox_ = BoundingBox(-M_LARGE_VALUE, M_LARGE_VALUE);
        defaultCell->outerBoundingBox_ = BoundingBox(-M_LARGE_VALUE, M_LARGE_VALUE);
        defaultCell->data_.Clear();

        // Reset cells
        unsigned index = 0;
        for (int x = 0; x < numCellsXYZ.x_; ++x)
        {
            for (int y = 0; y < numCellsXYZ.y_; ++y)
            {
                for (int z = 0; z < numCellsXYZ.z_; ++z)
                {
                    Vector3 floatIndex(static_cast<float>(x), static_cast<float>(y), static_cast<float>(z));
                    SceneGridCell& cell = cells_[index];
                    cell.index_ = IntVector3(x, y, z);
                    cell.innerBoundingBox_.min_ = boundingBox_.min_ + floatIndex * cellSize_;
                    cell.innerBoundingBox_.max_ = cell.innerBoundingBox_.min_ + Vector3::ONE * cellSize_;
                    cell.outerBoundingBox_.min_ = cell.innerBoundingBox_.min_ - Vector3::ONE * pad;
                    cell.outerBoundingBox_.max_ = cell.innerBoundingBox_.max_ + Vector3::ONE * pad;
                    cell.data_.Clear();
                    ++index;
                }
            }
        }
    }
    /// Update dirty drawables.
    void UpdateDirtyDrawables()
    {
        for (SceneGridCell& cell : cells_)
            UpdateDirtyDrawablesInCell(cell);
    }
    /// Reset all temporary marks.
    void ClearTemporary()
    {
        drawablesData_.ClearTemporary();
        for (SceneGridCell& cell : cells_)
            cell.data_.ClearTemporary();
    }
    /// Remove all drawables.
    void RemoveAllDrawables()
    {
        for (SceneGridCell& cell : cells_)
            RemoveAllDrawablesInCell(cell);
    }

    /// Query cells in frustum.
    void QueryCellsInFrustum(const Frustum& frustum, unsigned threadingThreshold, unsigned numThreads,
        SceneGridQueryResult& result)
    {
        result.Clear();

        result.cellMasks_.Resize(numCells_);

        // Test intersection and count all drawables inside
        unsigned numDrawables = 0;
        for (unsigned cellIndex = 0; cellIndex < numCells_; ++cellIndex)
        {
            result.cellMasks_[cellIndex] = OUTSIDE;
            if (cells_[cellIndex].data_.size_ > 0)
            {
                const Intersection intersection = frustum.IsInside(cells_[cellIndex].outerBoundingBox_);
                result.cellMasks_[cellIndex] = intersection;
                if (intersection != OUTSIDE)
                    numDrawables += cells_[cellIndex].data_.size_;
            }
        }

        // Apply threshold
        if (numDrawables < threadingThreshold)
            numThreads = 1;
        else if (numDrawables < numThreads)
            numThreads = numDrawables;

        // Append cells
        const unsigned maxDrawablesPerThread = (numDrawables - 1) / numThreads + 1;

        unsigned rangeBegin = 0;
        unsigned cellIndex = 0;
        unsigned beginIndex = 0;
        for (unsigned threadIndex = 0; threadIndex < numThreads; ++threadIndex)
        {
            unsigned remainingDrawables = maxDrawablesPerThread;
            while (remainingDrawables > 0 && cellIndex < numCells_)
            {
                SceneGridCell& cell = cells_[cellIndex];
                const unsigned cellSize = cell.data_.size_;

                // Skip empty cell
                if (result.cellMasks_[cellIndex] == OUTSIDE || beginIndex >= cellSize)
                {
                    ++cellIndex;
                    beginIndex = 0;
                    continue;
                }

                // Append cell
                SceneGridCellRef ref;
                ref.data_ = &cell.data_;
                ref.beginIndex_ = beginIndex;
                ref.endIndex_ = beginIndex + Min(cell.data_.size_ - beginIndex, remainingDrawables);
                ref.intersection_ = result.cellMasks_[cellIndex];
                result.cells_.Push(ref);

                remainingDrawables -= ref.endIndex_ - ref.beginIndex_;
                beginIndex = ref.endIndex_;
            }

            result.threadRanges_.Push({ rangeBegin, result.cells_.Size() });
            rangeBegin = result.cells_.Size();
        }
    }
    /// Process cells in sphere. Expected signature is `function(cell, isInside)`.
    template <class T>
    void ProcessCellsInSphere(const Sphere& sphere, T function)
    {
        for (unsigned cellIndex = 0; cellIndex < numCells_; ++cellIndex)
        {
            SceneGridCell& cell = cells_[cellIndex];
            if (cell.data_.size_ > 0)
            {
                const Intersection intersection = sphere.IsInside(cell.outerBoundingBox_);
                if (intersection != OUTSIDE)
                    function(cell.data_, intersection == INSIDE);
            }
        }
    }
    /// Process cells in frustum. Expected signature is `function(cell, isInside)`.
    template <class T>
    void ProcessCellsInFrustum(const Frustum& frustum, T function)
    {
        for (unsigned cellIndex = 0; cellIndex < numCells_; ++cellIndex)
        {
            SceneGridCell& cell = cells_[cellIndex];
            if (cell.data_.size_ > 0)
            {
                const Intersection intersection = frustum.IsInside(cell.outerBoundingBox_);
                if (intersection != OUTSIDE)
                    function(cell.data_, intersection == INSIDE);
            }
        }
    }

    /// Add drawable.
    void AddDrawable(Drawable* drawable)
    {
        DrawableIndex drawableIndex = drawable->GetDrawableIndex();
        assert(!drawableIndex || drawableIndex.grid_ == this);

        // Push global data and assign global index
        drawablesData_.Push(drawable);
        drawableIndex.grid_ = this;
        drawableIndex.gridIndex_ = drawablesData_.size_ - 1;

        // Insert into new cell and assign local index
        SceneGridCell* cell = FindCell(drawable->GetWorldBoundingBox());
        cell->data_.Push(drawable);
        drawableIndex.cell_ = cell;
        drawableIndex.cellIndex_ = cell->data_.size_ - 1;

        // Update index
        drawable->SetDrawableIndex(drawableIndex);
    }
    /// Remove drawable.
    void RemoveDrawable(Drawable* drawable)
    {
        const DrawableIndex drawableIndex = drawable->GetDrawableIndex();
        assert(drawableIndex.grid_ == this);

        // Remove zone
        if (drawable->GetDrawableFlags() & DRAWABLE_ZONE)
        {
            Zone* zone = static_cast<Zone*>(drawable);
            for (SceneGridCell& cell : cells_)
                cell.data_.ResetZone(zone);
        }

        // Remove drawable from cell
        SceneGridCellDrawableSoA& cellData = drawableIndex.cell_->data_;
        const unsigned index = drawableIndex.cellIndex_;

        cellData.EraseSwap(index);
        drawablesData_.EraseSwap(drawableIndex.gridIndex_);
    }
    /// Mark drawable dirty.
    void MarkDrawableDirty(Drawable* drawable)
    {
        const DrawableIndex drawableIndex = drawable->GetDrawableIndex();
        assert(drawableIndex.grid_ == this);

        SceneGridCell& cell = *drawableIndex.cell_;
        cell.data_.dirty_[drawableIndex.cellIndex_] = true;
    }
    /// Update drawable parameters.
    void UpdateDrawableParameters(Drawable* drawable)
    {
        const DrawableIndex drawableIndex = drawable->GetDrawableIndex();
        assert(drawableIndex.grid_ == this);

        SceneGridCellDrawableSoA& cellData = drawableIndex.cell_->data_;
        const unsigned index = drawableIndex.cellIndex_;

        cellData.UpdateDrawable(index, drawable);
    }
    /// Mark zone dirty.
    void MarkZoneDirty(Zone* zone)
    {
        for (SceneGridCell& cell : cells_)
            cell.data_.MarkZonesDirty();
    }

    /// Get global data.
    SceneGridDrawableSoA& GetGlobalData() { return drawablesData_; }
    /// Get default cell.
    SceneGridCell* GetDefaultCell() { return &cells_[numCells_ - 1]; }
    /// Find cell for specified bounding box.
    SceneGridCell* FindCell(const BoundingBox& drawableBoundingBox)
    {
        const Vector3 drawableCenter = drawableBoundingBox.Center();
        const IntVector3 cellIndex = VectorFloorToInt((drawableCenter - boundingBox_.min_) / boundingBox_.Size() * floatNumCellsXYZ_);
        if (cellIndex.x_ >= 0 && cellIndex.y_ >= 0 && cellIndex.z_ >= 0
            && cellIndex.x_ < numCellsXYZ_.x_ && cellIndex.y_ < numCellsXYZ_.y_ && cellIndex.z_ < numCellsXYZ_.z_)
        {
            const unsigned index = static_cast<unsigned>(cellIndex.x_ * numCellsXYZ_.y_ * numCellsXYZ_.z_
                + cellIndex.y_ * numCellsXYZ_.z_ + cellIndex.z_);
            if (cells_[index].IsInside(drawableBoundingBox))
            {
                return &cells_[index];
            }
        }
        return GetDefaultCell();
    }
    // TODO(eugeneko) Kill me
    SceneGridCellDrawableSoA& GetData() { return cells_[0].data_; }

protected:
    void UpdateDirtyDrawablesInCell(SceneGridCell& cell)
    {
        SceneGridCellDrawableSoA& data = cell.data_;
        for (unsigned i = 0; i < data.size_; ++i)
        {
            if (!data.dirty_[i])
                continue;

            Drawable* drawable = data.drawable_[i];
            data.dirty_[i] = false;
            data.transform_[i] = drawable->GetNode()->GetWorldTransform();
            data.boundingBox_[i] = drawable->GetWorldBoundingBox();
            data.boundingSphere_[i] = Sphere(data.boundingBox_[i]);

            data.cachedZoneDirty_[i] = true;

            // If this cell is default or drawable doesn't fit, update cell
            if (&cell == GetDefaultCell() || !cell.IsInside(data.boundingBox_[i]))
            {
                SceneGridCell* newCell = FindCell(data.boundingBox_[i]);
                if (newCell != &cell)
                {
                    DrawableIndex drawableIndex = drawable->GetDrawableIndex();

                    // Move to other cell
                    newCell->data_.Push(drawable);
                    cell.data_.EraseSwap(drawableIndex.cellIndex_);

                    drawableIndex.cell_ = newCell;
                    drawableIndex.cellIndex_ = newCell->data_.size_ - 1;
                    drawable->SetDrawableIndex(drawableIndex);
                }
            }

        }
    }
    static void RemoveAllDrawablesInCell(SceneGridCell& cell)
    {
        for (Drawable* drawable : cell.data_.drawable_)
            drawable->SetDrawableIndex(DrawableIndex{});
        cell.data_.Clear();
    }

    BoundingBox boundingBox_;
    IntVector3 numCellsXYZ_;
    unsigned numCells_{};
    Vector3 floatNumCellsXYZ_;
    Vector3 cellSize_;

    Vector<SceneGridCell> cells_;

    SceneGridDrawableSoA drawablesData_;
};

struct SceneQueryZonesAndOccludersResult
{
    DrawableVector occluders_;
    ZoneVector zones_;
    void Clear()
    {
        occluders_.Clear();
        zones_.Clear();
    }
    void Append(const SceneQueryZonesAndOccludersResult& other)
    {
        occluders_.Push(other.occluders_);
        zones_.Push(other.zones_);
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
    LightVector lights_;
    float minZ_{};
    float maxZ_{};
    void Clear()
    {
        lights_.Clear();
        minZ_ = M_INFINITY;
        maxZ_ = 0.0f;
    }
    void Append(const SceneQueryGeometriesAndLightsResult& other)
    {
        lights_.Push(other.lights_);

        minZ_ = Min(minZ_, other.minZ_);
        maxZ_ = Max(maxZ_, other.maxZ_);
    }
};

struct OldSceneQueryGeometriesAndLightsResult
{
    LightVector lights_;

    DrawableVector geometries_;
    Vector<unsigned> globalIndices_;
    ZoneVector zones_;
    Vector<unsigned> lightMasks_;
    Vector<unsigned> shadowMasks_;
    Vector<float> distances_;
    Vector<Vector2> minmaxZ_;
    Vector<BoundingBox> boundingBoxes_; // TODO(eugeneko) This is needed only for focusing

    float minZ_{};
    float maxZ_{};
    bool IsValid() const
    {
        const unsigned size = geometries_.Size();
        return size == geometries_.Size()
            && size == zones_.Size()
            && size == globalIndices_.Size()
            && size == lightMasks_.Size()
            && size == shadowMasks_.Size()
            && size == distances_.Size()
            && size == minmaxZ_.Size()
            && size == boundingBoxes_.Size()
            ;
    }
    void Clear()
    {
        lights_.Clear();

        geometries_.Clear();
        globalIndices_.Clear();
        zones_.Clear();
        lightMasks_.Clear();
        shadowMasks_.Clear();
        distances_.Clear();
        minmaxZ_.Clear();
        boundingBoxes_.Clear();
        minZ_ = M_INFINITY;
        maxZ_ = 0.0f;

        assert(IsValid());
    }
    void Merge(Vector<OldSceneQueryGeometriesAndLightsResult>& array)
    {
        Clear();
        if (array.Size() > 1)
        {
            for (OldSceneQueryGeometriesAndLightsResult& result : array)
            {
                lights_.Push(result.lights_);

                geometries_.Push(result.geometries_);
                globalIndices_.Push(result.globalIndices_);
                zones_.Push(result.zones_);
                lightMasks_.Push(result.lightMasks_);
                shadowMasks_.Push(result.shadowMasks_);
                distances_.Push(result.distances_);
                minmaxZ_.Push(result.minmaxZ_);
                boundingBoxes_.Push(result.boundingBoxes_);
                minZ_ = Min(minZ_, result.minZ_);
                maxZ_ = Max(maxZ_, result.maxZ_);
            }
        }
        else
        {
            Swap(lights_, array[0].lights_);

            Swap(geometries_, array[0].geometries_);
            Swap(globalIndices_, array[0].globalIndices_);
            Swap(zones_, array[0].zones_);
            Swap(lightMasks_, array[0].lightMasks_);
            Swap(shadowMasks_, array[0].shadowMasks_);
            Swap(distances_, array[0].distances_);
            Swap(minmaxZ_, array[0].minmaxZ_);
            Swap(boundingBoxes_, array[0].boundingBoxes_);
            minZ_ = array[0].minZ_;
            maxZ_ = array[0].maxZ_;
        }
        if (minZ_ == M_INFINITY)
            minZ_ = 0.0f;

        assert(IsValid());
    }
};

struct ZoneContext
{
    bool cameraZoneOverride_{};
    Zone* cameraZone_{};
    Zone* farClipZone_{};
    Zone* GetActualZone(Zone* drawableZone) const
    {
        if (cameraZoneOverride_)
            return cameraZone_;
        return drawableZone ? drawableZone : cameraZone_;
    }
};

/// Per-view drawable data.
struct ViewProcessorDrawableSoA
{
    unsigned Size() { return visible_.Size(); }
    void Resize(unsigned size)
    {
        visible_.Resize(size);
        distance_.Resize(size);
    }
    Vector<bool> visible_;
    Vector<float> distance_;
};

/// Contains view-specific scene data.
class ViewProcessor
{
public:
    virtual unsigned GetSceneQueryThreshold() { return 64; }
    virtual unsigned GetLightUpdateThreshold() { return 64; }
    /// Query zones and occluders.
    virtual void QuerySceneZonesAndOccluders(WorkQueue* workQueue, SceneGrid* sceneGrid,
        unsigned viewMask, const Frustum& frustum)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        sceneGrid->QueryCellsInFrustum(frustum, GetSceneQueryThreshold(), numThreads, zonesAndOccludersQuery_);

        zonesAndOccludersThreadResults_.Resize(numThreads);
        for (SceneQueryZonesAndOccludersResult& result : zonesAndOccludersThreadResults_)
            result.Clear();

        zonesAndOccludersQuery_.ScheduleWork(workQueue,
            [=](const SceneGridCellRef& cellRef, unsigned threadIndex)
        {
            SceneGridCellDrawableSoA& data = *cellRef.data_;
            SceneQueryZonesAndOccludersResult& result = zonesAndOccludersThreadResults_[threadIndex];
            for (unsigned index = cellRef.beginIndex_; index < cellRef.endIndex_; ++index)
            {
                if (!data.MatchViewMask(index, viewMask))
                    continue;

                Drawable* drawable = data.drawable_[index];

                // TODO(eugeneko) Get rid of branching
                if (data.IsZone(index))
                {
                    result.zones_.Push(static_cast<Zone*>(drawable));
                }
                else if (data.IsGeometry(index) && data.occluder_[index])
                {
                    if (cellRef.intersection_ == INSIDE || data.IsInFrustum(index, frustum))
                        result.occluders_.Push(drawable);
                }
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
        zonesAndOccluders_.Merge(zonesAndOccludersThreadResults_);
    }
    /// Prepare zones for further usage.
    virtual void CookZones(Camera* cullCamera, Zone* defaultZone)
    {
        ZoneVector& zones = zonesAndOccluders_.zones_;
        zonesData_.cameraZoneOverride_ = false;
        zonesData_.cameraZone_ = nullptr;
        zonesData_.farClipZone_ = nullptr;

        // Sort zones
        Sort(zones.Begin(), zones.End(),
            [](Zone* lhs, Zone* rhs) { return lhs->GetPriority() > rhs->GetPriority(); });

        // Find camera and far clip zones
        Node* cullCameraNode = cullCamera->GetNode();
        const Vector3 cameraPos = cullCameraNode->GetWorldPosition();
        const Vector3 farClipPos = cameraPos + cullCameraNode->GetWorldDirection() * Vector3(0.0f, 0.0f, cullCamera->GetFarClip());
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
    /// Query geometries.
    virtual void QuerySceneGeometriesAndLights(WorkQueue* workQueue, SceneGrid* sceneGrid,
        Camera* cullCamera, OcclusionBuffer* occlusionBuffer)
    {
        const unsigned viewMask = cullCamera->GetViewMask();
        const Frustum& frustum = cullCamera->GetFrustum();
        const Matrix3x4& viewMatrix = cullCamera->GetView();
        const Vector3 viewZ = Vector3(viewMatrix.m20_, viewMatrix.m21_, viewMatrix.m22_);
        const Vector3 absViewZ = viewZ.Abs();
        const unsigned cameraZoneLightMask = zonesData_.cameraZone_->GetLightMask();
        const unsigned cameraZoneShadowMask = zonesData_.cameraZone_->GetShadowMask();

        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        sceneGrid->QueryCellsInFrustum(frustum, GetSceneQueryThreshold(), numThreads, geometriesAndLightsQuery_);

        geometriesAndLightsThreadResults_.Resize(numThreads);
        for (OldSceneQueryGeometriesAndLightsResult& result : geometriesAndLightsThreadResults_)
            result.Clear();

        zonesAndOccludersQuery_.ScheduleWork(workQueue,
            [=](const SceneGridCellRef& cellRef, unsigned threadIndex)
        {
            SceneGridCellDrawableSoA& data = *cellRef.data_;
            OldSceneQueryGeometriesAndLightsResult& result = geometriesAndLightsThreadResults_[threadIndex];

            // TODO(eugeneko) Process lights here
            for (unsigned index = cellRef.beginIndex_; index < cellRef.endIndex_; ++index)
            {
                // Discard drawables from other views
                if (!data.MatchViewMask(index, viewMask))
                    continue;

                // Discard drawables except lights and geometries
                const bool isGeometry = data.IsGeometry(index);
                const bool isLigth = data.IsLight(index);
                if (!isGeometry && !isLigth)
                    continue;

                // Discard invisible
                if (cellRef.intersection_ == INTERSECTS && !data.IsInFrustum(index, frustum))
                    continue;

                // Calculate drawable distance
                const float drawDistance = data.drawDistance_[index];
                const float distance = cullCamera->GetDistance(data.boundingSphere_[index].center_);

                // Discard if drawable is too far
                if (drawDistance > 0.0f && distance > drawDistance)
                    continue;

                // Discard using occlusion buffer
                if (data.IsOccludedByBuffer(index, occlusionBuffer))
                    continue;

                // Update drawable zone
                Drawable* drawable = data.drawable_[index];
                if (isGeometry)
                {
                    if (!zonesData_.cameraZoneOverride_ && !GetZones().Empty())
                    {
                        UpdateDirtyZone(data, index, viewMask, frustum);
                    }

                    // Update min and max Z
                    float minZ = M_LARGE_VALUE;
                    float maxZ = M_LARGE_VALUE;

                    const BoundingBox& boundingBox = data.boundingBox_[index];
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
                    Zone* cachedZone = data.cachedZone_[index];
                    Zone* actualZone = cachedZone;
                    unsigned zoneLightMask = data.cachedZoneLightMask_[index];
                    unsigned zoneShadowMask = data.cachedZoneShadowMask_[index];

                    if (zonesData_.cameraZoneOverride_ || !cachedZone)
                    {
                        actualZone = zonesData_.cameraZone_;
                        zoneLightMask = cameraZoneLightMask;
                        zoneShadowMask = cameraZoneShadowMask;
                    }

                    // Get masks
                    const unsigned drawableLightMask = data.lightMask_[index];
                    const unsigned drawableShadowMask = data.shadowMask_[index];

                    // Push results
                    result.geometries_.Push(drawable);
                    result.globalIndices_.Push(index);
                    result.zones_.Push(actualZone);
                    result.lightMasks_.Push(drawableLightMask & zoneLightMask);
                    result.shadowMasks_.Push(drawableShadowMask & zoneShadowMask);
                    result.distances_.Push(distance);
                    result.minmaxZ_.Push(Vector2(minZ, maxZ));
                    result.boundingBoxes_.Push(boundingBox);
                }
                else //if (isLight)
                {
                    result.lights_.Push(static_cast<Light*>(drawable));
                }
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
        geometriesAndLights_.Merge(geometriesAndLightsThreadResults_);
    }
    /// Update visible lights.
    virtual void UpdateLights(WorkQueue* workQueue, const FrameInfo& frame)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;

        workQueue->ScheduleWork(GetLightUpdateThreshold(), geometriesAndLights_.lights_.Size(), numThreads,
            [this, &frame](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
        {
            for (unsigned index = beginIndex; index < endIndex; ++index)
            {
                // TODO(eugeneko) REMOVE!!
                Drawable* drawable = geometriesAndLights_.lights_[index];
                drawable->UpdateBatches(frame);
                drawable->MarkInView(frame);
            }
        });

        workQueue->Complete(M_MAX_UNSIGNED);
    }
    /// Prepare lights for further usage.
    virtual void CookLights(Camera* cullCamera)
    {
        LightVector& lights = geometriesAndLights_.lights_;
        const unsigned numLights = lights.Size();

        // Sort the lights to brightest/closest first, and per-vertex lights first so that per-vertex base pass can be evaluated first
        for (Light* light : lights)
        {
            light->SetIntensitySortValue(cullCamera->GetDistance(light->GetNode()->GetWorldPosition()));
//             light->SetLightQueue(nullptr);
        }

        Sort(lights.Begin(), lights.End(), CompareLights);
    }
    /// Reset visible objects data.
    virtual void ResetVisibleObjects(const SceneGridCellDrawableSoA& sceneDrawables)
    {
        // Reset array
        const unsigned numDrawables = sceneDrawables.size_;
        drawablesData_.Resize(numDrawables);
        for (unsigned i = 0; i < numDrawables; ++i)
            drawablesData_.visible_[i] = false;

        // Update distances for visible geometries
        const unsigned numVisibleGeometries = geometriesAndLights_.geometries_.Size();
        for (unsigned i = 0; i < numVisibleGeometries; ++i)
        {
            const unsigned drawableIndex = geometriesAndLights_.globalIndices_[i];
            drawablesData_.visible_[drawableIndex] = true;
            drawablesData_.distance_[drawableIndex] = geometriesAndLights_.distances_[i];
        }
    }

    const SceneQueryZonesAndOccludersResult& GetZonesAndOccluders() const { return zonesAndOccluders_; }
    const OldSceneQueryGeometriesAndLightsResult& GetGeometriesAndLights() const { return geometriesAndLights_; }
    const ZoneContext& GetCookedZonesData() const { return zonesData_; }
    const ZoneVector& GetZones() const { return zonesAndOccluders_.zones_; }
    const DrawableVector& GetOccluders() const { return zonesAndOccluders_.occluders_; }
    const DrawableVector& GetGeometries() const { return geometriesAndLights_.geometries_; }
    const LightVector& GetLights() const { return geometriesAndLights_.lights_; }

    const ViewProcessorDrawableSoA& GetDrawablesData() { return drawablesData_; }

protected:
    void UpdateDirtyZone(SceneGridCellDrawableSoA& cellData, unsigned index, unsigned viewMask, const Frustum& frustum) const
    {
        const Vector3 drawableCenter = cellData.boundingSphere_[index].center_;

        Zone* cachedZone = cellData.cachedZone_[index];
        const bool cachedZoneDirty = cellData.cachedZoneDirty_[index];
        const unsigned cachedZoneViewMask = cellData.cachedZoneViewMask_[index];

        // TODO(eugeneko) Is branch is not optimized for multiple zones
        if (cachedZoneDirty || !cachedZone || !(cachedZoneViewMask & viewMask))
        {
            // Find new zone
            const ZoneVector& zones = GetZones();
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

private:
    SceneGridQueryResult zonesAndOccludersQuery_;
    Vector<SceneQueryZonesAndOccludersResult> zonesAndOccludersThreadResults_;
    SceneQueryZonesAndOccludersResult zonesAndOccluders_;

    ZoneContext zonesData_;

    SceneGridQueryResult geometriesAndLightsQuery_;
    Vector<OldSceneQueryGeometriesAndLightsResult> geometriesAndLightsThreadResults_;
    OldSceneQueryGeometriesAndLightsResult geometriesAndLights_;

    ViewProcessorDrawableSoA drawablesData_;
};

struct DrawableUpdateManagerQueueSoA
{
    bool IsValid() const
    {
        return Size() == drawables_.Size()
            && Size() == distances_.Size()
            && Size() == zones_.Size()
            ;
    }
    unsigned Size() const { return drawables_.Size(); }
    void Clear()
    {
        drawables_.Clear();
        distances_.Clear();
        zones_.Clear();

        assert(IsValid());
    }
    void Push(Drawable* drawable, float distance, Zone* zone)
    {
        drawables_.Push(drawable);
        distances_.Push(distance);
        zones_.Push(zone);

        assert(IsValid());
    }
    DrawableVector drawables_;
    Vector<float> distances_;
    Vector<Zone*> zones_;
};

class DrawableBatchProcessor
{
public:
    virtual void UpdateDrawables(const FrameInfo& frame,
        const DrawableUpdateManagerQueueSoA& queue, unsigned beginIndex, unsigned endIndex)
    {
        for (unsigned index = beginIndex; index < endIndex; ++index)
        {
            Drawable* drawable = queue.drawables_[index];
            drawable->UpdateBatches(frame);
            drawable->SetZone(queue.zones_[index]);
            drawable->MarkInView(frame);
        }
    }
};

class DrawableUpdateManager
{
public:
    DrawableUpdateManager()
    {
        AddProcesor<DrawableBatchProcessor>(StringHash());
    }
    template <class T> void AddProcesor(StringHash drawableType)
    {
        keys_.Push(drawableType);
        queues_.Emplace();
        processors_.Push(MakeUnique<DrawableBatchProcessor>());
    }
    virtual unsigned GetUpdateThreshold() { return 16; }
    virtual void UpdateDrawables(WorkQueue* workQueue, const FrameInfo& frame,
        const SceneGridCellDrawableSoA& sceneData, const ViewProcessorDrawableSoA& viewData)
    {
        assert(sceneData.size_ == viewData.visible_.Size());

        FillQueues(sceneData, viewData);
        UpdateDrawablesInQueues(workQueue, frame);
    }
    unsigned FindQueue(const StringHash& drawableType) const
    {
        const unsigned index = keys_.IndexOf(drawableType);
        return index < keys_.Size() ? index : 0;
    }

protected:
    // TODO(eugeneko) Make it parallel too?
    virtual void FillQueues(const SceneGridCellDrawableSoA& sceneData, const ViewProcessorDrawableSoA& viewData)
    {
        // Clear queues
        for (DrawableUpdateManagerQueueSoA& queue : queues_)
            queue.Clear();

        // Fill queues with drawables
        for (unsigned i = 0; i < sceneData.size_; ++i)
        {
            if (viewData.visible_[i])
            {
                Drawable* drawable = sceneData.drawable_[i];
                const StringHash drawableType = sceneData.drawableType_[i];
                const unsigned queueIndex = FindQueue(drawableType);
                DrawableUpdateManagerQueueSoA& queue = queues_[queueIndex];
                queue.Push(drawable, viewData.distance_[i], sceneData.cachedZone_[i]);
            }
        }
    }
    virtual void UpdateDrawablesInQueues(WorkQueue* workQueue, const FrameInfo& frame)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        for (unsigned i = 0; i < queues_.Size(); ++i)
        {
            DrawableUpdateManagerQueueSoA& queue = queues_[i];
            assert(queue.IsValid());
            DrawableBatchProcessor* processor = processors_[i].Get();
            workQueue->ScheduleWork(GetUpdateThreshold(), queue.Size(), numThreads,
                [=, &queue, &frame](unsigned beginIndex, unsigned endIndex, unsigned threadIndex)
            {
                processor->UpdateDrawables(frame, queue, beginIndex, endIndex);
            });
        }
        workQueue->Complete(M_MAX_UNSIGNED);
    }

private:
    Vector<StringHash> keys_;
    Vector<DrawableUpdateManagerQueueSoA> queues_;
    Vector<UniquePtr<DrawableBatchProcessor>> processors_;
};

struct LightProcessorLightsSoA
{
    unsigned size_ = 0;

    Vector<Light*> lights_;
    Vector<bool> shadowed_;

    bool IsValid() const
    {
        return size_ == lights_.Size()
            && size_ == shadowed_.Size()
            ;
    }
    void Push(Light* light, bool drawShadows)
    {
        ++size_;

        lights_.Push(light);
        shadowed_.Push(drawShadows && IsLightShadowed(light));

        assert(IsValid());
    }
    void EraseSwap(unsigned index)
    {
        assert(index < size_);
        --size_;

        lights_.EraseSwap(index);
        shadowed_.EraseSwap(index);

        assert(IsValid());
    }
    void Clear()
    {
        size_ = 0;

        lights_.Clear();
        shadowed_.Clear();

        assert(IsValid());
    }

    static bool IsLightShadowed(Light* light)
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
};

class LightProcessor
{
public:
    virtual unsigned GetThreshold() { return 1; }
    virtual void SetupLights(const LightVector& visibleLights, bool drawShadows)
    {
        lights_.Clear();
        for (Light* light : visibleLights)
            lights_.Push(light, drawShadows);
    }
    virtual void QueryLitAndShadowGeometries(WorkQueue* workQueue, Renderer* renderer, Camera* cullCamera,
        const OldSceneQueryGeometriesAndLightsResult& visibleGeometries, const SceneGridCellDrawableSoA& sceneData)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;

        lightProcessingResults_.Resize(lights_.size_);
        for (unsigned i = 0; i < lights_.size_; ++i)
        {
            LightProcessingResult& result = lightProcessingResults_[i];
            Light* light = lights_.lights_[i];
            const bool isShadowed = lights_.shadowed_[i];
            const LightType lightType = light->GetLightType();

            result.Clear(light);
            switch (lightType)
            {
            case LIGHT_DIRECTIONAL:
                workQueue->ScheduleWork([=, &result, &visibleGeometries](unsigned /*threadIndex*/)
                {
                    // Gather lit geometries
                    const unsigned lightMask = light->GetLightMask();
                    const unsigned numGeometries = visibleGeometries.geometries_.Size();
                    for (unsigned i = 0; i < numGeometries; ++i)
                    {
                        if (visibleGeometries.lightMasks_[i] & lightMask)
                            result.litGeometries_.Push(visibleGeometries.geometries_[i]);
                    }

                    // Return if not shadowed
                    if (result.litGeometries_.Empty() || !isShadowed)
                        return;

                    // Setup shadow cameras
                    SetupShadowCameras(renderer, cullCamera, visibleGeometries, result);
                });
                break;
            default:
                assert(0);
                break;
            }
        }

        workQueue->Complete(M_MAX_UNSIGNED);
    }

    const Vector<LightProcessingResult>& GetProcessedLights() const { return lightProcessingResults_; }

public:
    /// Set up initial shadow camera view(s).
    static void SetupShadowCameras(Renderer* renderer, Camera* cullCamera,
        const OldSceneQueryGeometriesAndLightsResult& visibleGeometries, LightProcessingResult& query);
    /// Set up a directional light shadow camera.
    static void SetupDirLightShadowCamera(Camera* cullCamera, const OldSceneQueryGeometriesAndLightsResult& visibleGeometries,
        Camera* shadowCamera, Light* light, float nearSplit, float farSplit);
    /// Quantize a directional light shadow camera view to eliminate swimming.
    static void QuantizeDirLightShadowCamera(Camera* shadowCamera, Light* light, const IntRect& shadowViewport, const BoundingBox& viewBox);

private:
    Vector<bool> visibleObjects_;
    LightProcessorLightsSoA lights_;

    Vector<LightProcessingResult> lightProcessingResults_;
};

class DrawableProcessor : public Object
{
    URHO3D_OBJECT(DrawableProcessor, Object);

public:
    SceneGrid* sceneGrid_{};
    UniquePtr<ViewProcessor> viewProcessor_;
    UniquePtr<LightProcessor> lightProcessor_;

    UniquePtr<DrawableUpdateManager> drawableUpdateManager_;

    DrawableProcessor(Context* context) : Object(context)
    {
        viewProcessor_ = MakeUnique<ViewProcessor>();
        lightProcessor_ = MakeUnique<LightProcessor>();
        drawableUpdateManager_ = MakeUnique<DrawableUpdateManager>();
    }

    void Update(WorkQueue* workQueue, View* view)
    {
        const unsigned numThreads = workQueue->GetNumThreads() + 1;
        Camera* cullCamera = view->GetCullCamera();
        Renderer* renderer = view->GetRenderer();
        Zone* defaultZone = renderer->GetDefaultZone();

        {
            URHO3D_PROFILE(UpdateDirtyDrawables);
            sceneGrid_->UpdateDirtyDrawables();
        }

        {
            URHO3D_PROFILE(GetDrawables);
            viewProcessor_->QuerySceneZonesAndOccluders(workQueue, sceneGrid_,
                cullCamera->GetViewMask(), cullCamera->GetFrustum());
            viewProcessor_->CookZones(cullCamera, defaultZone);
            viewProcessor_->QuerySceneGeometriesAndLights(workQueue, sceneGrid_, cullCamera, nullptr);
            viewProcessor_->UpdateLights(workQueue, view->GetFrameInfo());
            viewProcessor_->CookLights(cullCamera);
//             viewProcessor_->ResetVisibleObjects(sceneGrid_->GetData());
        }

        {
            URHO3D_PROFILE(ProcessLights);
            lightProcessor_->SetupLights(viewProcessor_->GetLights(), view->drawShadows_);
            lightProcessor_->QueryLitAndShadowGeometries(workQueue, renderer, cullCamera,
                viewProcessor_->GetGeometriesAndLights(), sceneGrid_->GetData());
        }

        {
            URHO3D_PROFILE(UpdateBatches);
//             drawableUpdateManager_->UpdateDrawables(
//                 workQueue, view->GetFrameInfo(), sceneGrid_->GetData(), viewProcessor_->GetDrawablesData());
        }
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
