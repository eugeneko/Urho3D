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

#include <Urho3D/Core/CoreEvents.h>
#include <Urho3D/Core/Profiler.h>
#include <Urho3D/Engine/Engine.h>
#include <Urho3D/Graphics/Camera.h>
#include <Urho3D/Graphics/Graphics.h>
#include <Urho3D/Graphics/Material.h>
#include <Urho3D/Graphics/Model.h>
#include <Urho3D/Graphics/Octree.h>
#include <Urho3D/Graphics/Renderer.h>
#include <Urho3D/Graphics/StaticModelGroup.h>
#include <Urho3D/Graphics/Zone.h>
#include <Urho3D/Input/Input.h>
#include <Urho3D/Resource/ResourceCache.h>
#include <Urho3D/Scene/Scene.h>
#include <Urho3D/UI/Font.h>
#include <Urho3D/UI/Text.h>
#include <Urho3D/UI/UI.h>

#include "HugeObjectCount.h"

#include <Urho3D/DebugNew.h>

URHO3D_DEFINE_APPLICATION_MAIN(HugeObjectCount)

HugeObjectCount::HugeObjectCount(Context* context) :
    Sample(context),
    animate_(false),
    singleMaterial_(false)
{
}

void HugeObjectCount::Start()
{
    // Execute base class startup
    Sample::Start();

    DebugHud* debugHud = GetSubsystem<DebugHud>();
    debugHud->SetMode(DEBUGHUD_SHOW_ALL);

    // Create the scene content
    CreateScene();

    // Create the UI content
    CreateInstructions();

    // Setup the viewport for displaying the scene
    SetupViewport();

    // Hook up to the frame update events
    SubscribeToEvents();

    // Set the mouse mode to use in the sample
    Sample::InitMouseMode(MM_RELATIVE);
}

static const bool USE_CLASSIC_SCENE = 1;
static const float XZ_GEOMETRY_RANGE = 800.0f;
static const float XZ_LIGHT_RANGE = 400.0f;
static const int NUM_OBJECTS = 125;
static const int NUM_LIGHTS = 32;
static const int NUM_MATERIALS = 256;

void HugeObjectCount::CreateScene()
{
    auto* cache = GetSubsystem<ResourceCache>();

    if (!scene_)
        scene_ = new Scene(context_);
    else
    {
        scene_->Clear();
        boxNodes_.Clear();
    }

    // Create the Octree component to the scene so that drawable objects can be rendered. Use default volume
    // (-1000, -1000, -1000) to (1000, 1000, 1000)
    scene_->CreateComponent<Octree>();

    // Create a Zone for ambient light & fog control
    Node* zoneNode = scene_->CreateChild("Zone");
    auto* zone = zoneNode->CreateComponent<Zone>();
    zone->SetBoundingBox(BoundingBox(-1000.0f, 1000.0f));
    zone->SetFogColor(Color(0.4f, 0.4f, 0.4f));
    zone->SetFogStart(200.0f);
    zone->SetFogEnd(300.0f);

    // Create a directional light
    Node* lightNode = scene_->CreateChild("DirectionalLight");
    lightNode->SetDirection(Vector3(-0.6f, -1.0f, -0.8f)); // The direction vector does not need to be normalized
    auto* light = lightNode->CreateComponent<Light>();
    light->SetLightType(LIGHT_DIRECTIONAL);
    light->SetCastShadows(true);
    light->SetShadowBias(BiasParameters(0.00025f, 0.5f));
    light->SetShadowCascade(CascadeParameters(10.0f, 50.0f, 200.0f, 0.0f, 0.8f));

    if (0)
    {
        Node* boxNode = scene_->CreateChild("BigBox");
        boxNode->SetPosition(Vector3(0.0f, 20.0f, 0.0f));
        boxNode->SetScale(20.0f);
        auto* boxObject = boxNode->CreateComponent<StaticModel>();
        boxObject->SetModel(cache->GetResource<Model>("Models/Box.mdl"));
        boxObject->SetMaterial(cache->GetResource<Material>("Materials/GreenTransparent.xml"));
        boxObject->SetOccluder(true);
        boxNodes_.Push(SharedPtr<Node>(boxNode));
    }

    {
        light->SetColor(Color(0.5f, 0.5f, 0.5f));

        // Create geometries
//         auto model = cache->GetResource<Model>("Models/Sphere.mdl");
        auto model = cache->GetResource<Model>("Models/Box.mdl");
        auto baseMaterial = cache->GetResource<Material>("Materials/StoneTiled.xml");
        Vector<SharedPtr<Material>> materials;
        for (unsigned i = 0; i < NUM_MATERIALS; ++i)
        {
            Color color;
            color.r_ = Random(0.1f, 1.0f);
            color.g_ = Random(0.1f, 1.0f);
            color.b_ = Random(0.1f, 1.0f);

            auto material = baseMaterial->Clone();
            material->SetShaderParameter("MatSpecColor", color);
            materials.Push(material);
        }
        for (int y = 0; y < NUM_OBJECTS; ++y)
        {
            for (int x = 0; x < NUM_OBJECTS; ++x)
            {
                unsigned materialIndex = (unsigned)Random((int)materials.Size());
                float factor = Vector2(x, y).Length() / Vector2(NUM_OBJECTS, NUM_OBJECTS).Length();
                float scale = Pow(factor, 1.5f);

                Vector3 position;
                position.x_ = scale * x * XZ_GEOMETRY_RANGE / NUM_OBJECTS;
                position.y_ = Random(-10.0f, 10.0f);
                position.z_ = scale * y * XZ_GEOMETRY_RANGE / NUM_OBJECTS;

                Node* boxNode = scene_->CreateChild("Box");
                boxNode->SetPosition(position);
                boxNode->SetScale(Max(1.0f, scale * 10.0f));
                auto* boxObject = boxNode->CreateComponent<StaticModel>();
                boxObject->SetModel(model);
                boxObject->SetMaterial(singleMaterial_ ? baseMaterial : materials[materialIndex]);
                boxObject->SetCastShadows(true);
                boxNodes_.Push(SharedPtr<Node>(boxNode));
            }
        }

        // Create lights
        {
            for (unsigned i = 0; i < NUM_LIGHTS; ++i)
            {
                Vector3 position;
                position.x_ = Random(0.0f, XZ_LIGHT_RANGE);
                position.y_ = 5.0f;
                position.z_ = Random(0.0f, XZ_LIGHT_RANGE);
                float radius = Random(15.0f, 30.0f);

                Node* lightNode = scene_->CreateChild("PointLight");
                lightNode->SetPosition(position);
                auto* light = lightNode->CreateComponent<Light>();
                light->SetLightType(LIGHT_POINT);
                light->SetRange(radius);
            }
        }
    }

    // Create the camera. Create it outside the scene so that we can clear the whole scene without affecting it
    if (!cameraNode_)
    {
        cameraNode_ = new Node(context_);
        cameraNode_->SetPosition(Vector3(-10.0f, 15.0f, -10.0f));
        cameraNode_->SetDirection(Vector3(4, -1, 4));
        auto* camera = cameraNode_->CreateComponent<Camera>();
        camera->SetFarClip(600.0f);
        pitch_ = cameraNode_->GetWorldRotation().EulerAngles().x_;
        yaw_ = cameraNode_->GetWorldRotation().EulerAngles().y_;
    }
}

void HugeObjectCount::CreateInstructions()
{
    auto* cache = GetSubsystem<ResourceCache>();
    auto* ui = GetSubsystem<UI>();

    // Construct new Text object, set string to display and font to use
    auto* instructionText = ui->GetRoot()->CreateChild<Text>();
    instructionText->SetText(
        "Use WASD keys and mouse/touch to move\n"
        "Space to toggle animation\n"
        "G to toggle object group optimization"
    );
    instructionText->SetFont(cache->GetResource<Font>("Fonts/Anonymous Pro.ttf"), 15);
    // The text has multiple rows. Center them in relation to each other
    instructionText->SetTextAlignment(HA_CENTER);

    // Position the text relative to the screen center
    instructionText->SetHorizontalAlignment(HA_CENTER);
    instructionText->SetVerticalAlignment(VA_CENTER);
    instructionText->SetPosition(0, ui->GetRoot()->GetHeight() / 4);
}

void HugeObjectCount::SetupViewport()
{
    auto* renderer = GetSubsystem<Renderer>();

    // Set up a viewport to the Renderer subsystem so that the 3D scene can be seen
    SharedPtr<Viewport> viewport(new Viewport(context_, scene_, cameraNode_->GetComponent<Camera>()));
    renderer->SetViewport(0, viewport);
}

void HugeObjectCount::SubscribeToEvents()
{
    // Subscribe HandleUpdate() function for processing update events
    SubscribeToEvent(E_UPDATE, URHO3D_HANDLER(HugeObjectCount, HandleUpdate));
}

void HugeObjectCount::MoveCamera(float timeStep)
{
    // Do not move if the UI has a focused element (the console)
    if (GetSubsystem<UI>()->GetFocusElement())
        return;

    auto* input = GetSubsystem<Input>();

    // Movement speed as world units per second
    const float MOVE_SPEED = 20.0f;
    // Mouse sensitivity as degrees per pixel
    const float MOUSE_SENSITIVITY = 0.1f;

    // Use this frame's mouse motion to adjust camera node yaw and pitch. Clamp the pitch between -90 and 90 degrees
    IntVector2 mouseMove = input->GetMouseMove();
    yaw_ += MOUSE_SENSITIVITY * mouseMove.x_;
    pitch_ += MOUSE_SENSITIVITY * mouseMove.y_;
    pitch_ = Clamp(pitch_, -90.0f, 90.0f);

    // Construct new orientation for the camera scene node from yaw and pitch. Roll is fixed to zero
    cameraNode_->SetRotation(Quaternion(pitch_, yaw_, 0.0f));

    // Read WASD keys and move the camera scene node to the corresponding direction if they are pressed
    if (input->GetKeyDown(KEY_W))
        cameraNode_->Translate(Vector3::FORWARD * MOVE_SPEED * timeStep);
    if (input->GetKeyDown(KEY_S))
        cameraNode_->Translate(Vector3::BACK * MOVE_SPEED * timeStep);
    if (input->GetKeyDown(KEY_A))
        cameraNode_->Translate(Vector3::LEFT * MOVE_SPEED * timeStep);
    if (input->GetKeyDown(KEY_D))
        cameraNode_->Translate(Vector3::RIGHT * MOVE_SPEED * timeStep);
}

void HugeObjectCount::AnimateObjects(float timeStep)
{
    URHO3D_PROFILE(AnimateObjects);

    const float ROTATE_SPEED = 15.0f;
    // Rotate about the Z axis (roll)
    Quaternion rotateQuat(ROTATE_SPEED * timeStep, Vector3::FORWARD);

    for (unsigned i = 0; i < boxNodes_.Size(); ++i)
        boxNodes_[i]->Rotate(rotateQuat);
}

void HugeObjectCount::HandleUpdate(StringHash eventType, VariantMap& eventData)
{
    using namespace Update;

    // Take the frame time step, which is stored as a float
    float timeStep = eventData[P_TIMESTEP].GetFloat();

    // Toggle animation with space
    auto* input = GetSubsystem<Input>();
    if (input->GetKeyPress(KEY_SPACE))
        animate_ = !animate_;

    // Toggle grouped / ungrouped mode
    if (input->GetKeyPress(KEY_G))
    {
        singleMaterial_ = !singleMaterial_;
        CreateScene();
    }

    // Move the camera, scale movement with time step
    MoveCamera(timeStep);

    // Animate scene if enabled
    if (animate_)
        AnimateObjects(timeStep);
}
