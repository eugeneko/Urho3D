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

// #include "../Core/Context.h"
// #include "../Core/CoreEvents.h"
// #include "../Core/Profiler.h"
// #include "../Core/Thread.h"
// #include "../Core/WorkQueue.h"
// #include "../Graphics/DebugRenderer.h"
// #include "../Graphics/Graphics.h"
#include "../Graphics/DrawableProcessor.h"
// #include "../IO/Log.h"
// #include "../Scene/Scene.h"
// #include "../Scene/SceneEvents.h"

#include "../DebugNew.h"

namespace Urho3D
{

void LightProcessor::SetupDirLightShadowCamera(Camera* cullCamera, const SceneQueryGeometriesAndLightsResult& visibleGeometries,
    Camera* shadowCamera, Light* light, float nearSplit, float farSplit)
{
    Node* shadowCameraNode = shadowCamera->GetNode();
    Node* lightNode = light->GetNode();
    const float extrusionDistance = Min(cullCamera->GetFarClip(), light->GetShadowMaxExtrusion());
    const FocusParameters& parameters = light->GetShadowFocus();

    // Calculate initial position & rotation
    const Vector3 pos = cullCamera->GetNode()->GetWorldPosition() - extrusionDistance * lightNode->GetWorldDirection();
    shadowCameraNode->SetTransform(pos, lightNode->GetWorldRotation());

    // Calculate main camera shadowed frustum in light's view space
    farSplit = Min(farSplit, cullCamera->GetFarClip());
    // Use the scene Z bounds to limit frustum size if applicable
    if (parameters.focus_)
    {
        nearSplit = Max(visibleGeometries.minZ_, nearSplit);
        farSplit = Min(visibleGeometries.maxZ_, farSplit);
    }

    Frustum splitFrustum = cullCamera->GetSplitFrustum(nearSplit, farSplit);
    Polyhedron frustumVolume;
    frustumVolume.Define(splitFrustum);
    // If focusing enabled, clip the frustum volume by the combined bounding box of the lit geometries within the frustum
    if (parameters.focus_)
    {
        BoundingBox litGeometriesBox;
        const unsigned lightMask = light->GetLightMask();

        const unsigned numGeometries = visibleGeometries.geometries_.Size();
        for (unsigned i = 0; i < numGeometries; ++i)
        {
            const Vector2& minMaxZ = visibleGeometries.minmaxZ_[i];
            if (minMaxZ.x_ <= farSplit && minMaxZ.y_ >= nearSplit && (visibleGeometries.lightMasks_[i] & lightMask))
                litGeometriesBox.Merge(visibleGeometries.boundingBoxes_[i]);
        }

        if (litGeometriesBox.Defined())
        {
            frustumVolume.Clip(litGeometriesBox);
            // If volume became empty, restore it to avoid zero size
            if (frustumVolume.Empty())
                frustumVolume.Define(splitFrustum);
        }
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

void LightProcessor::QuantizeDirLightShadowCamera(Camera* shadowCamera, Light* light, const IntRect& shadowViewport,
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

void LightProcessor::SetupShadowCameras(Renderer* renderer, Camera* cullCamera,
    const SceneQueryGeometriesAndLightsResult& visibleGeometries, LightProcessingResult& query)
{
    Light* light = query.light_;

    unsigned splits = 0;

    switch (light->GetLightType())
    {
    case LIGHT_DIRECTIONAL:
        {
            const CascadeParameters& cascade = light->GetShadowCascade();

            float nearSplit = cullCamera->GetNearClip();
            float farSplit;
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
                Camera* shadowCamera = renderer->GetShadowCamera();
                query.shadowCameras_[splits] = shadowCamera;
                query.shadowNearSplits_[splits] = nearSplit;
                query.shadowFarSplits_[splits] = farSplit;
                SetupDirLightShadowCamera(cullCamera, visibleGeometries, shadowCamera, light, nearSplit, farSplit);

                nearSplit = farSplit;
                ++splits;
            }
        }
        break;

    case LIGHT_SPOT:
        {
            Camera* shadowCamera = renderer->GetShadowCamera();
            query.shadowCameras_[0] = shadowCamera;
            Node* cameraNode = shadowCamera->GetNode();
            Node* lightNode = light->GetNode();

            cameraNode->SetTransform(lightNode->GetWorldPosition(), lightNode->GetWorldRotation());
            shadowCamera->SetNearClip(light->GetShadowNearFarRatio() * light->GetRange());
            shadowCamera->SetFarClip(light->GetRange());
            shadowCamera->SetFov(light->GetFov());
            shadowCamera->SetAspectRatio(light->GetAspectRatio());

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
                Camera* shadowCamera = renderer->GetShadowCamera();
                query.shadowCameras_[i] = shadowCamera;
                Node* cameraNode = shadowCamera->GetNode();

                // When making a shadowed point light, align the splits along X, Y and Z axes regardless of light rotation
                cameraNode->SetPosition(light->GetNode()->GetWorldPosition());
                cameraNode->SetDirection(*directions[i]);
                shadowCamera->SetNearClip(light->GetShadowNearFarRatio() * light->GetRange());
                shadowCamera->SetFarClip(light->GetRange());
                shadowCamera->SetFov(90.0f);
                shadowCamera->SetAspectRatio(1.0f);
            }

            splits = MAX_CUBEMAP_FACES;
        }
        break;
    }

    query.numSplits_ = splits;
}

}
