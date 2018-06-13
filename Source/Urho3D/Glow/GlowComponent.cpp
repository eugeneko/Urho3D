//
// Copyright (c) 2014-2017 THUNDERBEAST GAMES LLC
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

#include "../IO/Log.h"
#include "../Core/Context.h"
#include "../Core/WorkQueue.h"
#include "../Graphics/StaticModel.h"
#include "../Resource/XMLFile.h"

// #include <ToolCore/ToolSystem.h>
// #include <ToolCore/Project/Project.h>

#include "Common/GlowSettings.h"
#include "Common/GlowEvents.h"
#include "Kernel/BakeModel.h"
#include "Kernel/BakeMaterial.h"
#include "Kernel/SceneBaker.h"
#include "GlowComponent.h"

// using namespace ToolCore;
using namespace AtomicGlow;

namespace Urho3D
{

GlowComponent::GlowComponent(Context *context) : Component(context)
{
    GlowSettings glowSettings;
    glowSettings.SetDefaults();
    SetFromGlowSettings(glowSettings);
}

GlowComponent::~GlowComponent()
{

}

void GlowComponent::SetFromGlowSettings(const GlowSettings& settings)
{

    lexelDensity_ = settings.lexelDensity_;

    giEnabled_ = settings.giEnabled_;
    giGranularity_ = settings.giGranularity_;
    giMaxBounces_ = settings.giMaxBounces_;

    aoEnabled_ = settings.aoEnabled_;
    aoDepth_ = settings.aoDepth_;
    nsamples_ = settings.nsamples_;
    aoMin_ = settings.aoMin_;
    aoMultiply_ = settings.aoMultiply_;

}

void GlowComponent::CopyToGlowSettings(GlowSettings& settings) const
{
    settings.lexelDensity_ = lexelDensity_ ;
    settings.sceneLexelDensityScale_ = 1.0f;

    settings.giEnabled_= giEnabled_;
    settings.giGranularity_ = giGranularity_;
    settings.giMaxBounces_ = giMaxBounces_;

    settings.aoEnabled_ = aoEnabled_;
    settings.aoDepth_ = aoDepth_;
    settings.nsamples_ = nsamples_;
    settings.aoMin_ = aoMin_;
    settings.aoMultiply_ = aoMultiply_;
}

bool GlowComponent::stub_ = false;

bool GlowComponent::Bake()
{
    if (!GetSubsystem<BakeModelCache>())
    {
        context_->RegisterSubsystem<BakeModelCache>();
    }
    if (!GetSubsystem<BakeMaterialCache>())
    {
        context_->RegisterSubsystem<BakeMaterialCache>();
    }

    const String projectPath = "D:/";
    auto sceneCopy = MakeShared<XMLFile>(context_);
    GetScene()->SaveXML(sceneCopy->GetOrCreateRoot("scene"));

    auto sceneBaker = MakeShared<SceneBaker>(context_, projectPath);
    sceneBaker->SetStandaloneMode(false);
    if (!sceneBaker->LoadScene(sceneCopy->GetRoot()))
    {
        assert(0);
        return false;
    }

    while (sceneBaker->GetCurrentLightMode() != GLOW_LIGHTMODE_COMPLETE)
    {
        GetSubsystem<WorkQueue>()->Complete(M_MAX_UNSIGNED);

        if (sceneBaker->GetCurrentLightMode() == GLOW_LIGHTMODE_UNDEFINED)
        {

            // light mode will either move to direct or complete, depending on work to do
            sceneBaker->Light(GLOW_LIGHTMODE_DIRECT);
            continue;
        }

        if (sceneBaker->GetCurrentLightMode() == GLOW_LIGHTMODE_DIRECT)
        {
            sceneBaker->LightFinishCycle();

            // light mode will either move to indirect or complete, depending on work to do
            sceneBaker->Light(GLOW_LIGHTMODE_INDIRECT);
            continue;
        }

        if (sceneBaker->GetCurrentLightMode() == GLOW_LIGHTMODE_INDIRECT)
        {
            sceneBaker->LightFinishCycle();

            // light mode will either move to indirect or complete, depending on work to do
            sceneBaker->Light(GLOW_LIGHTMODE_INDIRECT);
            continue;
        }
    }

    // if we're done, exit in standalone mode, or send IPC event if
    // running off Glow service
    sceneBaker->GenerateLightmaps();

    VectorBuffer bakeData = sceneBaker->GetBakeData();
    bakeData.Seek(0);

    unsigned count = bakeData.ReadUInt();

    PODVector<Node*> children;

    GetScene()->GetChildrenWithComponent<StaticModel>(children, true);

    for (unsigned i = 0; i < count; i++)
    {
        Node* node = 0;
        StaticModel* staticModel = 0;

        unsigned nodeID = bakeData.ReadUInt();
        unsigned staticModelID = bakeData.ReadUInt();
        unsigned lightMask = bakeData.ReadUInt();
        unsigned lightmapIndex = bakeData.ReadUInt();
        Vector4 tilingOffset = bakeData.ReadVector4();

        for (unsigned j = 0; j < children.Size(); j++)
        {
            if (children[j]->GetID() == nodeID)
            {
                node = children[j];
                staticModel = node->GetComponent<StaticModel>();

                if (!staticModel || staticModel->GetID() != staticModelID)
                {
                    URHO3D_LOGERROR("GlowService::ProcessBakeData() - mismatched node <-> static model ID");
                    return false;
                }

                break;
            }

        }

        if (!node)
        {
            URHO3D_LOGERROR("GlowService::ProcessBakeData() - unable to find node ID");
            return false;
        }

        staticModel->SetLightMask(lightMask);
        staticModel->SetLightmapIndex(lightmapIndex);
        staticModel->SetLightmapTilingOffset(tilingOffset);
    }

    return true;
}

void GlowComponent::HandleAtomicGlowBakeCancel(StringHash eventType, VariantMap& eventData)
{
//     GetSubsystem<GlowService>()->CancelBake();
}

void GlowComponent::HandleAtomicGlowServiceBakeResult(StringHash eventType, VariantMap& eventData)
{
    using namespace AtomicGlowServiceBakeResult;

    // convert to a glow component event, which contains the same fields
    SendEvent(E_ATOMICGLOWBAKERESULT, eventData);

}

void GlowComponent::HandleAtomicGlowServiceLogEvent(StringHash eventType, VariantMap& eventData)
{
    using namespace AtomicGlowServiceLogEvent;

    // convert to a glow component event, which contains the same fields
    SendEvent(E_ATOMICGLOWLOGEVENT, eventData);

}

void GlowComponent::RegisterObject(Context* context)
{
    context->RegisterFactory<GlowComponent>();

    URHO3D_ATTRIBUTE("Lexel Density", float, lexelDensity_, 0.1f, AM_FILE);

    URHO3D_ATTRIBUTE("GI Enabled", bool, giEnabled_, false, AM_FILE);
    URHO3D_ATTRIBUTE("GI Granularity", int, giGranularity_, 16, AM_FILE);
    URHO3D_ATTRIBUTE("GI Max Cycles", int, giMaxBounces_, 3, AM_FILE);

    URHO3D_ATTRIBUTE("AO Enabled", bool, aoEnabled_, false, AM_FILE);
    URHO3D_ATTRIBUTE("AO Depth", float, aoDepth_, 0.25f, AM_FILE);
    URHO3D_ATTRIBUTE("AO Samples", int, nsamples_, 64, AM_FILE);
    URHO3D_ATTRIBUTE("AO Min", float, aoMin_, 0.45f, AM_FILE);
    URHO3D_ATTRIBUTE("AO Multiply", float, aoMultiply_, 1.0f, AM_FILE);

    URHO3D_ATTRIBUTE("Bake", bool, stub_, false, AM_EDIT);

}


}
