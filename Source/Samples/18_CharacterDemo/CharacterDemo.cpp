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
#include <Urho3D/Core/ProcessUtils.h>
#include <Urho3D/Engine/Engine.h>
#include <Urho3D/Graphics/AnimatedModel.h>
#include <Urho3D/Graphics/AnimationController.h>
#include <Urho3D/Graphics/Camera.h>
#include <Urho3D/Graphics/Geometry.h>
#include <Urho3D/Graphics/IndexBuffer.h>
#include <Urho3D/Graphics/Light.h>
#include <Urho3D/Graphics/Material.h>
#include <Urho3D/Graphics/Octree.h>
#include <Urho3D/Graphics/Renderer.h>
#include <Urho3D/Graphics/Terrain.h>
#include <Urho3D/Graphics/VertexBuffer.h>
#include <Urho3D/Graphics/Zone.h>
#include <Urho3D/Input/Controls.h>
#include <Urho3D/Input/Input.h>
#include <Urho3D/IO/FileSystem.h>
#include <Urho3D/Physics/CollisionShape.h>
#include <Urho3D/Physics/PhysicsWorld.h>
#include <Urho3D/Physics/RigidBody.h>
#include <Urho3D/Resource/ResourceCache.h>
#include <Urho3D/Scene/Scene.h>
#include <Urho3D/UI/Font.h>
#include <Urho3D/UI/Text.h>
#include <Urho3D/UI/UI.h>

#include "Character.h"
#include "CharacterDemo.h"
#include "Touch.h"

#include <random>

#include <Urho3D/DebugNew.h>


URHO3D_DEFINE_APPLICATION_MAIN(CharacterDemo)

//////////////////////////////////////////////////////////////////////////
using PointCloud2D = PODVector<Vector2>;
using PointCloud2DNorm = PODVector<Vector2>;
class PoissonRandom
{
public:
    PoissonRandom(unsigned seed);
    ~PoissonRandom();
    PointCloud2DNorm generate(float minDist, unsigned newPointsCount, unsigned numPoints);
private:
    struct Cell;
    struct Grid;
    float randomFloat();
    static IntVector2 imageToGrid(const Vector2& p, float cellSize);
    Vector2 popRandom(PointCloud2DNorm& points);
    Vector2 generateRandomPointAround(const Vector2& p, float minDist);
private:
    struct Core;
    UniquePtr<Core> impl_;
};

struct PoissonRandom::Core
{
    Core(unsigned seed) : m_generator(seed), m_distr(0.0, 1.0) {}
    std::mt19937 m_generator;
    std::uniform_real_distribution<> m_distr;
};

// Point Cloud
PointCloud2D samplePointCloud(const PointCloud2DNorm& cloud,
    const Vector2& begin, const Vector2& end,
    float scale)
{
    PointCloud2D dest;
    const Vector2 from = VectorFloor(begin / scale);
    const Vector2 to = VectorCeil(end / scale);
    for (float nx = from.x_; nx <= to.x_; ++nx)
    {
        for (float ny = from.y_; ny <= to.y_; ++ny)
        {
            const Vector2 tileBegin = Vector2(nx, ny);
            const Vector2 tileEnd = Vector2(nx + 1, ny + 1);
            const Vector2 clipBegin = VectorMax(begin / scale, VectorMin(end / scale, tileBegin));
            const Vector2 clipEnd = VectorMax(begin / scale, VectorMin(end / scale, tileEnd));
            for (const Vector2& sourcePoint : cloud)
            {
                const Vector2 point = VectorLerp(tileBegin, tileEnd, static_cast<Vector2>(sourcePoint));
                if (point.x_ < clipBegin.x_ || point.y_ < clipBegin.y_ || point.x_ > clipEnd.x_ || point.y_ > clipEnd.y_)
                    continue;

                dest.Push(point * scale);
            }
        }
    }
    return dest;
}


// Ctor
PoissonRandom::PoissonRandom(unsigned seed)
    : impl_(MakeUnique<Core>(seed))
{
}
PoissonRandom::~PoissonRandom()
{
}


// Generator
struct PoissonRandom::Cell
    : public Vector2
{
    bool isValid = false;
};
float PoissonRandom::randomFloat()
{
    return static_cast<float>(impl_->m_distr(impl_->m_generator));
}
IntVector2 PoissonRandom::imageToGrid(const Vector2& p, float cellSize)
{
    return IntVector2(static_cast<int>(p.x_ / cellSize), static_cast<int>(p.y_ / cellSize));

}
struct PoissonRandom::Grid
{
    Grid(const int w, const int h, const float cellSize)
        : m_width(w)
        , m_height(h)
        , m_cellSize(cellSize)
    {
        m_grid.Resize(m_height);

        for (auto i = m_grid.Begin(); i != m_grid.End(); i++)
        {
            i->Resize(m_width);
        }
    }
    void insert(const Vector2& p)
    {
        IntVector2 g = imageToGrid(p, m_cellSize);
        Cell& c = m_grid[g.x_][g.y_];

        c.x_ = p.x_;
        c.y_ = p.y_;
        c.isValid = true;
    }
    bool isInNeighbourhood(const Vector2& point, const float minDist, const float cellSize)
    {
        IntVector2 g = imageToGrid(point, cellSize);

        // Number of adjucent cells to look for neighbour points
        const int d = 5;

        // Scan the neighbourhood of the point in the grid
        for (int i = g.x_ - d; i < g.x_ + d; i++)
        {
            for (int j = g.y_ - d; j < g.y_ + d; j++)
            {
                // Wrap cells
                int wi = i;
                int wj = j;
                while (wi < 0) wi += m_width;
                while (wj < 0) wj += m_height;
                wi %= m_width;
                wj %= m_height;

                // Test wrapped distances
                Cell p = m_grid[wi][wj];
                const float dist = Min(Min((p - point).Length(),
                    (p + Vector2(1, 0) - point).Length()),
                    Min((p + Vector2(0, 1) - point).Length(),
                    (p + Vector2(1, 1) - point).Length()));

                if (p.isValid && dist < minDist)
                    return true;
            }
        }
        return false;
    }

private:
    int m_width;
    int m_height;
    float m_cellSize;

    Vector<Vector<Cell>> m_grid;
};
Vector2 PoissonRandom::popRandom(PointCloud2DNorm& points)
{
    std::uniform_int_distribution<> dis(0, points.Size() - 1);
    const int idx = dis(impl_->m_generator);
    const Vector2 p = points[idx];
    points.Erase(points.Begin() + idx);
    return p;
}
Vector2 PoissonRandom::generateRandomPointAround(const Vector2& p, float minDist)
{
    // Start with non-uniform distribution
    const float r1 = randomFloat();
    const float r2 = randomFloat();

    // Radius should be between MinDist and 2 * MinDist
    const float radius = minDist * (r1 + 1.0f);

    // Random angle
    const float angle = 2 * 3.141592653589f * r2;

    // The new point is generated around the point (x, y)
    const float x = p.x_ + radius * cos(angle);
    const float y = p.y_ + radius * sin(angle);

    return Vector2(x, y);
}
PointCloud2DNorm PoissonRandom::generate(const float minDist, const unsigned newPointsCount, const unsigned numPoints)
{
    PointCloud2DNorm samplePoints;
    PointCloud2DNorm processList;

    // Create the grid
    const float cellSize = minDist / sqrt(2.0f);

    const int gridW = static_cast<int>(ceil(1.0f / cellSize));
    const int gridH = static_cast<int>(ceil(1.0f / cellSize));

    Grid grid(gridW, gridH, cellSize);

    Vector2 firstPoint;
    firstPoint.x_ = randomFloat();
    firstPoint.y_ = randomFloat();

    // Update containers
    processList.Push(firstPoint);
    samplePoints.Push(firstPoint);
    grid.insert(firstPoint);

    // Generate new points for each point in the queue
    while (!processList.Empty() && samplePoints.Size() < numPoints)
    {
        const Vector2 point = popRandom(processList);

        for (unsigned i = 0; i < newPointsCount; i++)
        {
            const Vector2 newPoint = generateRandomPointAround(point, minDist);

            // Test
            const bool fits = newPoint.x_ >= 0 && newPoint.y_ >= 0 && newPoint.x_ <= 1 && newPoint.y_ <= 1;

            if (fits && !grid.isInNeighbourhood(newPoint, minDist, cellSize))
            {
                processList.Push(newPoint);
                samplePoints.Push(newPoint);
                grid.insert(newPoint);
                continue;
            }
        }
    }
    return samplePoints;
}
//////////////////////////////////////////////////////////////////////////
CharacterDemo::CharacterDemo(Context* context) :
    Sample(context),
    firstPerson_(false)
{
    // Register factory and attributes for the Character component so it can be created via CreateComponent, and loaded / saved
    Character::RegisterObject(context);
}

CharacterDemo::~CharacterDemo() = default;

void CharacterDemo::Start()
{
    // Execute base class startup
    Sample::Start();
    if (touchEnabled_)
        touch_ = new Touch(context_, TOUCH_SENSITIVITY);

    // Create static scene content
    CreateScene();

    // Create the controllable character
    CreateCharacter();

    // Create the UI content
    CreateInstructions();

    // Subscribe to necessary events
    SubscribeToEvents();

    // Set the mouse mode to use in the sample
    Sample::InitMouseMode(MM_RELATIVE);
}

void CharacterDemo::CreateScene()
{
    auto* cache = GetSubsystem<ResourceCache>();

    scene_ = new Scene(context_);

    auto sceneData = cache->GetResource<XMLFile>("ForestScene/Scene.xml");
    scene_->LoadXML(sceneData->GetRoot());

    // Create scene subsystem components
    scene_->CreateComponent<PhysicsWorld>();

    // Create camera and define viewport. We will be doing load / save, so it's convenient to create the camera outside the scene,
    // so that it won't be destroyed and recreated, and we don't have to redefine the viewport on load
    cameraNode_ = new Node(context_);
    auto* camera = cameraNode_->CreateComponent<Camera>();
    camera->SetFarClip(300.0f);
    GetSubsystem<Renderer>()->SetViewport(0, new Viewport(context_, scene_, camera));

    // Generate grass geometries
    Node* terrainNode = scene_->GetChild("Terrain", true);
    auto terrain = terrainNode->GetComponent<Terrain>();

    unsigned totalGrassInstances = 0;
    const int NUM_GRASS_CHUNKS = 8;
    const float MAX_GRASS_ANGLE = 25.0f;
    if (Node* zoneNode = scene_->GetChild("GrassRegion", true))
    {
        PODVector<float> vertexData;
        PODVector<unsigned short> indexData;

        PoissonRandom random(Rand());
        auto points = random.generate(0.02f, 10, 20000);

        auto zone = zoneNode->GetComponent<Zone>();
        const BoundingBox boundingBox = zone->GetWorldBoundingBox();
        IntVector2 chunkIndex;
        for (chunkIndex.x_ = 0; chunkIndex.x_ < NUM_GRASS_CHUNKS; ++chunkIndex.x_)
        {
            const float chunkSize = (boundingBox.max_.x_ - boundingBox.min_.x_) / static_cast<float>(NUM_GRASS_CHUNKS);
            for (chunkIndex.y_ = 0; chunkIndex.y_ < NUM_GRASS_CHUNKS; ++chunkIndex.y_)
            {
                const Vector2 chunkBegin(boundingBox.min_.x_ + chunkIndex.x_ * chunkSize,
                    boundingBox.min_.z_ + chunkIndex.y_ * chunkSize);

                // Prepare buffers
                auto vertexBuffer = MakeShared<VertexBuffer>(context_);
                auto indexBuffer = MakeShared<IndexBuffer>(context_);
                auto geometry = MakeShared<Geometry>(context_);
                geometry->SetVertexBuffer(0, vertexBuffer);
                geometry->SetIndexBuffer(indexBuffer);

                auto model = MakeShared<Model>(context_);
                model->SetNumGeometries(1);
                model->SetNumGeometryLodLevels(0, 1);
                model->SetGeometry(0, 0, geometry);

                // Generate grass instances
                PODVector<Vector3> positions;
                PODVector<Vector3> normals;
                PODVector<Quaternion> rotations;
                for (unsigned i = 0; i < points.Size(); ++i)
                {
                    Vector3 position{ points[i].x_ * chunkSize + chunkBegin.x_ , 0.0f, points[i].y_ * chunkSize + chunkBegin.y_ };
                    position.y_ = terrain->GetHeight(position);
                    const Vector3 normal = terrain->GetNormal(position);
                    const Quaternion rotation = Quaternion(Vector3::UP, normal)
                        * Quaternion(Random(0.0f, 360.0f), Vector3::UP)
                        * Quaternion(30.0f, Vector3::RIGHT);

                    if (normal.DotProduct(Vector3::UP) < Cos(MAX_GRASS_ANGLE))
                        continue;

                    positions.Push(position);
                    normals.Push(normal);
                    rotations.Push(rotation);
                }

                const unsigned numBillboards = positions.Size();
                totalGrassInstances += numBillboards;

                // Fill buffers
                vertexData.Resize(numBillboards * 4 * 8);
                indexData.Resize(numBillboards * 6);

                unsigned short* indexPtr = indexData.Buffer();
                unsigned vertexIndex = 0;
                for (unsigned i = 0; i < numBillboards; ++i)
                {
                    indexPtr[0] = (unsigned short)vertexIndex;
                    indexPtr[1] = (unsigned short)(vertexIndex + 1);
                    indexPtr[2] = (unsigned short)(vertexIndex + 2);
                    indexPtr[3] = (unsigned short)(vertexIndex + 2);
                    indexPtr[4] = (unsigned short)(vertexIndex + 3);
                    indexPtr[5] = (unsigned short)vertexIndex;

                    indexPtr += 6;
                    vertexIndex += 4;
                }

                BoundingBox boundingBox;
                float* vertexPtr = vertexData.Buffer();
                for (unsigned i = 0; i < numBillboards; ++i)
                {
                    const Matrix3 rotationMatrix = rotations[i].RotationMatrix();
                    const Vector3 xAxis = rotationMatrix.Column(0);
                    const Vector3 yAxis = rotationMatrix.Column(1);
                    static const Vector2 uvs[4] = { Vector2::UP, Vector2::ONE, Vector2::RIGHT, Vector2::ZERO };
                    for (unsigned j = 0; j < 4; ++j)
                    {
                        const Vector3 pos = positions[i] + xAxis * (uvs[j].x_ - 0.5f) + yAxis * (1.0f - uvs[j].y_);
                        boundingBox.Merge(pos);
                        vertexPtr[0] = pos.x_;
                        vertexPtr[1] = pos.y_;
                        vertexPtr[2] = pos.z_;
                        vertexPtr[3] = normals[i].x_;
                        vertexPtr[4] = normals[i].y_;
                        vertexPtr[5] = normals[i].z_;
                        vertexPtr[6] = uvs[j].x_;
                        vertexPtr[7] = uvs[j].y_;
                        vertexPtr += 8;
                    }
                }

                // Update GPU
                vertexBuffer->SetSize(static_cast<unsigned>(vertexData.Size() / 8),
                    MASK_POSITION | MASK_NORMAL /*| MASK_COLOR*/ | MASK_TEXCOORD1 /*| MASK_TEXCOORD2*/, true);
                vertexBuffer->SetData(vertexData.Buffer());
                indexBuffer->SetSize(static_cast<unsigned>(indexData.Size()), false, true);
                indexBuffer->SetData(indexData.Buffer());
                geometry->SetDrawRange(TRIANGLE_LIST, 0, indexData.Size(), false);
                model->SetBoundingBox(boundingBox);

                // Make node
                Node* grassChunk = zoneNode->CreateChild("GrassChunk");
                auto grassStaticModel = grassChunk->CreateComponent<StaticModel>();
                grassStaticModel->SetModel(model);
                grassStaticModel->SetMaterial(cache->GetResource<Material>("ForestScene/Grass/Grass_mat.xml"));
                grassStaticModel->SetCastShadows(true);
            }
        }
    }
    URHO3D_LOGINFOF("Num grass instances: %d", totalGrassInstances);
}

void CharacterDemo::CreateCharacter()
{
    auto* cache = GetSubsystem<ResourceCache>();

    Node* objectNode = scene_->CreateChild("Jack");
    objectNode->SetPosition(Vector3(0.0f, 1.0f, 0.0f));

    // Spin node
    Node* adjustNode = objectNode->CreateChild("AdjNode");
    adjustNode->SetRotation( Quaternion(180, Vector3(0,1,0) ) );

    // Create rigidbody, and set non-zero mass so that the body becomes dynamic
    auto* body = objectNode->CreateComponent<RigidBody>();
    body->SetCollisionLayer(1);
    body->SetMass(1.0f);

    // Set zero angular factor so that physics doesn't turn the character on its own.
    // Instead we will control the character yaw manually
    body->SetAngularFactor(Vector3::ZERO);

    // Set the rigidbody to signal collision also when in rest, so that we get ground collisions properly
    body->SetCollisionEventMode(COLLISION_ALWAYS);

    // Set a capsule shape for collision
    auto* shape = objectNode->CreateComponent<CollisionShape>();
    shape->SetCapsule(0.7f, 1.8f, Vector3(0.0f, 0.9f, 0.0f));

    // Create the character logic component, which takes care of steering the rigidbody
    // Remember it so that we can set the controls. Use a WeakPtr because the scene hierarchy already owns it
    // and keeps it alive as long as it's not removed from the hierarchy
    character_ = objectNode->CreateComponent<Character>();
}

void CharacterDemo::CreateInstructions()
{
    auto* cache = GetSubsystem<ResourceCache>();
    auto* ui = GetSubsystem<UI>();

    // Construct new Text object, set string to display and font to use
    auto* instructionText = ui->GetRoot()->CreateChild<Text>();
    instructionText->SetText(
        "Use WASD keys and mouse/touch to move\n"
        "Space to jump, F to toggle 1st/3rd person\n"
        "F5 to save scene, F7 to load"
    );
    instructionText->SetFont(cache->GetResource<Font>("Fonts/Anonymous Pro.ttf"), 15);
    // The text has multiple rows. Center them in relation to each other
    instructionText->SetTextAlignment(HA_CENTER);

    // Position the text relative to the screen center
    instructionText->SetHorizontalAlignment(HA_CENTER);
    instructionText->SetVerticalAlignment(VA_CENTER);
    instructionText->SetPosition(0, ui->GetRoot()->GetHeight() / 4);
}

void CharacterDemo::SubscribeToEvents()
{
    // Subscribe to Update event for setting the character controls before physics simulation
    SubscribeToEvent(E_UPDATE, URHO3D_HANDLER(CharacterDemo, HandleUpdate));

    // Subscribe to PostUpdate event for updating the camera position after physics simulation
    SubscribeToEvent(E_POSTUPDATE, URHO3D_HANDLER(CharacterDemo, HandlePostUpdate));

    // Unsubscribe the SceneUpdate event from base class as the camera node is being controlled in HandlePostUpdate() in this sample
    UnsubscribeFromEvent(E_SCENEUPDATE);
}

void CharacterDemo::HandleUpdate(StringHash eventType, VariantMap& eventData)
{
    using namespace Update;

    auto* input = GetSubsystem<Input>();

    if (character_)
    {
        // Clear previous controls
        character_->controls_.Set(CTRL_FORWARD | CTRL_BACK | CTRL_LEFT | CTRL_RIGHT | CTRL_JUMP, false);

        // Update controls using touch utility class
        if (touch_)
            touch_->UpdateTouches(character_->controls_);

        // Update controls using keys
        auto* ui = GetSubsystem<UI>();
        if (!ui->GetFocusElement())
        {
            if (!touch_ || !touch_->useGyroscope_)
            {
                character_->controls_.Set(CTRL_FORWARD, input->GetKeyDown(KEY_W));
                character_->controls_.Set(CTRL_BACK, input->GetKeyDown(KEY_S));
                character_->controls_.Set(CTRL_LEFT, input->GetKeyDown(KEY_A));
                character_->controls_.Set(CTRL_RIGHT, input->GetKeyDown(KEY_D));
            }
            character_->controls_.Set(CTRL_JUMP, input->GetKeyDown(KEY_SPACE));

            // Add character yaw & pitch from the mouse motion or touch input
            if (touchEnabled_)
            {
                for (unsigned i = 0; i < input->GetNumTouches(); ++i)
                {
                    TouchState* state = input->GetTouch(i);
                    if (!state->touchedElement_)    // Touch on empty space
                    {
                        auto* camera = cameraNode_->GetComponent<Camera>();
                        if (!camera)
                            return;

                        auto* graphics = GetSubsystem<Graphics>();
                        character_->controls_.yaw_ += TOUCH_SENSITIVITY * camera->GetFov() / graphics->GetHeight() * state->delta_.x_;
                        character_->controls_.pitch_ += TOUCH_SENSITIVITY * camera->GetFov() / graphics->GetHeight() * state->delta_.y_;
                    }
                }
            }
            else
            {
                character_->controls_.yaw_ += (float)input->GetMouseMoveX() * YAW_SENSITIVITY;
                character_->controls_.pitch_ += (float)input->GetMouseMoveY() * YAW_SENSITIVITY;
            }
            // Limit pitch
            character_->controls_.pitch_ = Clamp(character_->controls_.pitch_, -80.0f, 80.0f);
            // Set rotation already here so that it's updated every rendering frame instead of every physics frame
            character_->GetNode()->SetRotation(Quaternion(character_->controls_.yaw_, Vector3::UP));

            // Switch between 1st and 3rd person
            if (input->GetKeyPress(KEY_F))
                firstPerson_ = !firstPerson_;

            // Turn on/off gyroscope on mobile platform
            if (touch_ && input->GetKeyPress(KEY_G))
                touch_->useGyroscope_ = !touch_->useGyroscope_;
        }
    }
}

void CharacterDemo::HandlePostUpdate(StringHash eventType, VariantMap& eventData)
{
    if (!character_)
        return;

    Node* characterNode = character_->GetNode();

    // Get camera lookat dir from character yaw + pitch
    const Quaternion& rot = characterNode->GetRotation();
    Quaternion dir = rot * Quaternion(character_->controls_.pitch_, Vector3::RIGHT);
    cameraNode_->SetPosition(characterNode->GetWorldPosition() + Vector3(0.0f, 1.7f, 0.0f));
    cameraNode_->SetRotation(dir);
}
