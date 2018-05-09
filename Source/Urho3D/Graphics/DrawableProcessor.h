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

#include "../Graphics/Drawable.h"
#include "../Graphics/StaticModel.h"
#include "../Graphics/Light.h"
#include "../Scene/Node.h"

namespace Urho3D
{

struct DefaultDrawableProcessorSoA
{
    Vector<Drawable*> drawable_;

    Vector<bool> transformDirty_;
    Vector<Matrix3x4> transform_;
    Vector<BoundingBox> worldBoundingBox_;
    Vector<Sphere> worldBoundingSphere_;

    Vector<unsigned> drawableFlag_;
    Vector<unsigned> viewMask_;
    Vector<bool> occluder_;

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
    }
    unsigned Size() const { return drawable_.Size(); }
};

class DefaultDrawableProcessor
{
public:
    virtual ~DefaultDrawableProcessor()
    {
        Clear();
    }
    const DefaultDrawableProcessorSoA& GetData() const { return data_; }

    virtual void AddDrawable(Drawable* drawable)
    {
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

    void GetAllDrawables(Vector<Drawable*>& drawables)
    {
        drawables.Push(data_.drawable_);
    }

    void MarkTransformDirty(unsigned index) { data_.transformDirty_[index] = true; }
    void SetViewMask(unsigned index, unsigned viewMask) { data_.viewMask_[index] = viewMask; }
    void SetOccluder(unsigned index, bool occluder) { data_.occluder_[index] = occluder; }
    unsigned GetNumDrawables() const { return data_.Size(); }

    virtual void Clear()
    {
        for (Drawable* drawable : data_.drawable_)
            drawable->SetDrawableIndex(DrawableIndex{});
        data_.Clear();
    }
    virtual void GetDrawables(PODVector<Drawable*>& result, unsigned drawableFlags, unsigned viewMask, const Frustum& frustum)
    {
        const unsigned numDrawables = GetNumDrawables();
        for (unsigned i = 0; i < numDrawables; ++i)
        {
            if (CheckDrawableFlags(i, drawableFlags, viewMask)
                && CheckDrawableInFrustum(i, frustum))
            {
                result.Push(data_.drawable_[i]);
            }
        }
    }
    virtual void GetZonesAndOccluders(PODVector<Drawable*>& result, unsigned viewMask, const Frustum& frustum)
    {
        const unsigned numDrawables = GetNumDrawables();
        for (unsigned i = 0; i < numDrawables; ++i)
        {
            if (CheckDrawableZoneOrOccluder(i, viewMask)
                && CheckDrawableInFrustum(i, frustum))
            {
                result.Push(data_.drawable_[i]);
            }
        }
    }
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

    bool CheckDrawableFlags(unsigned index, unsigned drawableFlags, unsigned viewMask) const
    {
        return (data_.drawableFlag_[index] & drawableFlags) && (data_.viewMask_[index] & viewMask);
    }
    bool CheckDrawableZoneOrOccluder(unsigned index, unsigned viewMask)
    {
        const unsigned flags = data_.drawableFlag_[index];
        return (data_.viewMask_[index] & viewMask)
            && (flags == DRAWABLE_ZONE || (flags == DRAWABLE_GEOMETRY && data_.occluder_[index]));
    }
    bool CheckDrawableInFrustum(unsigned index, const Frustum& frustum)
    {
        const Sphere& boundingSpere = data_.worldBoundingSphere_[index];
        const BoundingBox& boundingBox = data_.worldBoundingBox_[index];

        const Intersection sphereIntersection = frustum.IsInside(boundingSpere);
        return sphereIntersection != OUTSIDE
            && (sphereIntersection == INSIDE || frustum.IsInsideFast(boundingBox));
    }
public:
    DefaultDrawableProcessorSoA data_;
};

class StaticModelProcessor : public DefaultDrawableProcessor
{

};

class LightModelProcessor : public DefaultDrawableProcessor
{

};

class DrawableProcessor
{
public:
    bool AddDrawable(Drawable* drawable)
    {
        if (drawable->IsInstanceOf<StaticModel>())
        {
            staticModelProcessor_.AddDrawable(drawable);
            return true;
        }
        else
        {
            defaultProcessor_.AddDrawable(drawable);
            return true;
        }
        return false;
    }

    bool RemoveDrawable(Drawable* drawable)
    {
        if (drawable->IsInstanceOf<StaticModel>())
        {
            staticModelProcessor_.RemoveDrawable(drawable);
            return true;
        }
        else
        {
            defaultProcessor_.RemoveDrawable(drawable);
            return true;
        }
    }

    void GetDrawables(PODVector<Drawable*>& result, unsigned drawableFlags, unsigned viewMask, const Frustum& frustum)
    {
        staticModelProcessor_.GetDrawables(result, drawableFlags, viewMask, frustum);
        defaultProcessor_.GetDrawables(result, drawableFlags, viewMask, frustum);
    }
    void GetZonesAndOccluders(PODVector<Drawable*>& result, unsigned viewMask, const Frustum& frustum)
    {
        staticModelProcessor_.GetZonesAndOccluders(result, viewMask, frustum);
        defaultProcessor_.GetZonesAndOccluders(result, viewMask, frustum);
    }

    void GetAllDrawables(Vector<Drawable*>& drawables)
    {
        staticModelProcessor_.GetAllDrawables(drawables);
        defaultProcessor_.GetAllDrawables(drawables);
    }

    void UpdateDirtyDrawables()
    {
        staticModelProcessor_.UpdateDirtyDrawables();
        defaultProcessor_.UpdateDirtyDrawables();
    }

    StaticModelProcessor staticModelProcessor_;
//     LightModelProcessor lightProcessor_;
    DefaultDrawableProcessor defaultProcessor_;
private:
};

}
