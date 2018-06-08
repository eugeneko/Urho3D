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
#include <Urho3D/Physics/PhysicsEvents.h>
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
static const unsigned LAYER_MAIN = 1 << 0;
static const unsigned LAYER_GRASSPUSHER = 1 << 1;
static const unsigned LAYER_OBSTACLE = 1 << 2;
const float GRASS_BILLBOARD_SIZE = 0.7f;

//////////////////////////////////////////////////////////////////////////
float CubicSmooth(float x)
{
    return x * x * (3.0 - 2.0 * x);
}

float TriangleWave(float x)
{
    return Abs((Fract(x + 0.5) * 2.0) - 1.0);
}

float SmoothTriangleWave(float x)
{
    return CubicSmooth(TriangleWave(x));
}

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

PoissonRandom::PoissonRandom(unsigned seed)
    : impl_(MakeUnique<Core>(seed))
{
}

PoissonRandom::~PoissonRandom()
{
}

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
    Sample(context)
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

    // Collect grass pushers
    PODVector<RigidBody*> rigidBodies;
    scene_->GetComponents<RigidBody>(rigidBodies, true);
    for (RigidBody* body : rigidBodies)
    {
        if (body->GetCollisionLayer() & LAYER_GRASSPUSHER)
        {
            for (Component* component : body->GetNode()->GetComponents())
            {
                if (component->IsInstanceOf<CollisionShape>())
                {
                    grassPushers_.Push(static_cast<CollisionShape*>(component));
                }
            }
        }
    }

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
    Node* grassRegionNode = scene_->GetChild("GrassRegion", true);
    if (terrainNode && grassRegionNode)
    {
        grassLightMaterial_ = cache->GetResource<Material>("ForestScene/Grass/DynamicGrass_mat.xml");
        grassDarkMaterial_ = cache->GetResource<Material>("ForestScene/Grass/DynamicGrass_mat2.xml");
        CreateGrass(terrainNode, grassRegionNode);
    }

    darkTheme_ = true;
    SetupTheme(false);
}

void CharacterDemo::CreateGrass(Node* terrainNode, Node* grassRegionNode)
{
    auto* cache = GetSubsystem<ResourceCache>();

    auto terrain = terrainNode->GetComponent<Terrain>();

    unsigned totalGrassInstances = 0;
    const int NUM_GRASS_CHUNKS = 8;
    const float MAX_GRASS_ANGLE = 25.0f;
    PODVector<float> vertexData;
    PODVector<unsigned short> indexData;

    PoissonRandom random(Rand());
    auto points = random.generate(0.02f, 10, 20000);

    auto zone = grassRegionNode->GetComponent<Zone>();
    grassBoundingBox_ = zone->GetWorldBoundingBox();

    IntVector2 chunkIndex;
    for (chunkIndex.x_ = 0; chunkIndex.x_ < NUM_GRASS_CHUNKS; ++chunkIndex.x_)
    {
        const float chunkSize = (grassBoundingBox_.max_.x_ - grassBoundingBox_.min_.x_) / static_cast<float>(NUM_GRASS_CHUNKS);
        for (chunkIndex.y_ = 0; chunkIndex.y_ < NUM_GRASS_CHUNKS; ++chunkIndex.y_)
        {
            const Vector2 chunkBegin(grassBoundingBox_.min_.x_ + chunkIndex.x_ * chunkSize,
                grassBoundingBox_.min_.z_ + chunkIndex.y_ * chunkSize);

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
            static const unsigned VERTEX_STRIDE_FLOATS = 14;
            vertexData.Resize(numBillboards * 4 * VERTEX_STRIDE_FLOATS);
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

            BoundingBox modelBoundingBox;
            float* vertexPtr = vertexData.Buffer();
            for (unsigned i = 0; i < numBillboards; ++i)
            {
                const Vector3 basePos = positions[i];
                const float wave = (SmoothTriangleWave(basePos.x_ * 0.07f) + SmoothTriangleWave(basePos.z_ * 0.03f)) / 2;
                const float size = GRASS_BILLBOARD_SIZE * (Random(0.6f, 1.0f) + wave * 0.3);
                Matrix3 scaleMatrix;
                scaleMatrix.SetScale(size);
                const Matrix3 rotationScaleMatrix = rotations[i].RotationMatrix() * scaleMatrix;
                const Vector3 xAxis = rotationScaleMatrix.Column(0);
                const Vector3 yAxis = rotationScaleMatrix.Column(1);
                static const Vector2 uvs[4] = { Vector2::UP, Vector2::ONE, Vector2::RIGHT, Vector2::ZERO };
                for (unsigned j = 0; j < 4; ++j)
                {
                    const Vector3 pos = basePos + xAxis * (uvs[j].x_ - 0.5f) + yAxis * (1.0f - uvs[j].y_);
                    modelBoundingBox.Merge(pos);
                    vertexPtr[0] = pos.x_;
                    vertexPtr[1] = pos.y_;
                    vertexPtr[2] = pos.z_;
                    vertexPtr[3] = normals[i].x_;
                    vertexPtr[4] = normals[i].y_;
                    vertexPtr[5] = normals[i].z_;
                    vertexPtr[6] = uvs[j].x_;
                    vertexPtr[7] = uvs[j].y_;
                    vertexPtr[8] = (pos.x_ - grassBoundingBox_.min_.x_) / grassBoundingBox_.Size().x_;
                    vertexPtr[9] = (pos.z_ - grassBoundingBox_.min_.z_) / grassBoundingBox_.Size().z_;
                    vertexPtr[10] = basePos.x_;
                    vertexPtr[11] = basePos.y_;
                    vertexPtr[12] = basePos.z_;
                    vertexPtr[13] = (1.0f - uvs[j].y_) * GRASS_BILLBOARD_SIZE;
                    vertexPtr += VERTEX_STRIDE_FLOATS;
                }
            }

            // Update GPU
            vertexBuffer->SetSize(static_cast<unsigned>(vertexData.Size() / VERTEX_STRIDE_FLOATS),
                MASK_POSITION | MASK_NORMAL /*| MASK_COLOR*/ | MASK_TEXCOORD1 | MASK_TEXCOORD2 | MASK_TANGENT, true);
            vertexBuffer->SetData(vertexData.Buffer());
            indexBuffer->SetSize(static_cast<unsigned>(indexData.Size()), false, true);
            indexBuffer->SetData(indexData.Buffer());
            geometry->SetDrawRange(TRIANGLE_LIST, 0, indexData.Size(), false);
            model->SetBoundingBox(modelBoundingBox);

            // Make node
            Node* grassChunk = grassRegionNode->CreateChild("GrassChunk");
            auto grassStaticModel = grassChunk->CreateComponent<StaticModel>();
            grassStaticModel->SetModel(model);
            grassStaticModel->SetMaterial(grassLightMaterial_);
            grassStaticModel->SetCastShadows(true);
            grassStaticModel->SetLightMask(0xff & ~(1 << 0));
        }
    }
    URHO3D_LOGINFOF("Num grass instances: %d", totalGrassInstances);

    grassTexture_ = MakeShared<Texture2D>(context_);
    grassTexture_->SetNumLevels(1);
    grassTexture_->SetSize(grassTextureSize_, grassTextureSize_, Graphics::GetRGBAFormat(), TEXTURE_DYNAMIC);
    grassLightMaterial_->SetTexture(TU_NORMAL, grassTexture_);
    grassDarkMaterial_->SetTexture(TU_NORMAL, grassTexture_);
    grassTextureData_.Resize(grassTextureSize_ * grassTextureSize_ * 4);
    grassPushiness_.Resize(grassTextureSize_ * grassTextureSize_);
    grassBaseHeight_.Resize(grassTextureSize_ * grassTextureSize_);
    grassScale_.Resize(grassTextureSize_ * grassTextureSize_, 1.0f);

    auto physicsWorld = scene_->GetComponent<PhysicsWorld>();
    for (unsigned y = 0; y < grassTextureSize_; ++y)
    {
        for (unsigned x = 0; x < grassTextureSize_; ++x)
        {
            const unsigned i = y * grassTextureSize_ + x;

            const float kx = (x + 0.5f) / grassTextureSize_;
            const float ky = (y + 0.5f) / grassTextureSize_;
            const Vector3 samplePosition = VectorLerp(grassBoundingBox_.min_, grassBoundingBox_.max_, Vector3(kx, 0.0f, ky));
            grassBaseHeight_[i] = terrain->GetHeight(samplePosition);
        }
    }

    const float treeRadiusScale = 1.2f;
    Node* treesNode = scene_->GetChild("Trees");
    for (Node* treeNode : treesNode->GetChildren())
    {
        const Vector3 position = treeNode->GetWorldPosition();
        const float radius = treeNode->GetComponent<StaticModel>()->GetWorldBoundingBox().HalfSize().x_ * treeRadiusScale;
        for (unsigned y = 0; y < grassTextureSize_; ++y)
        {
            for (unsigned x = 0; x < grassTextureSize_; ++x)
            {
                const unsigned i = y * grassTextureSize_ + x;
                const float kx = (x + 0.5f) / grassTextureSize_;
                const float ky = (y + 0.5f) / grassTextureSize_;
                Vector3 samplePosition = VectorLerp(grassBoundingBox_.min_, grassBoundingBox_.max_, Vector3(kx, 0.0f, ky));
                samplePosition.y_ = grassBaseHeight_[i];
                const float distance = (samplePosition - position).Length();
                grassScale_[i] -= Pow(std::exp(-distance / radius), 2.0f);
            }
        }
    }
}

void CharacterDemo::SetupTheme(bool dark)
{
    if (darkTheme_ == dark)
        return;

    auto* cache = GetSubsystem<ResourceCache>();

    // Replace skull
    Node* ballNode = scene_->GetChild("Ball", true);
    Node* skullNode = scene_->GetChild("Skull", true);
    if (!darkTheme_)
    {
        const Vector3 position = ballNode->GetWorldPosition();
        ballNode->SetDeepEnabled(false);
        skullNode->SetDeepEnabled(true);
        skullNode->SetWorldPosition(position);
    }
    else
    {
        const Vector3 position = skullNode->GetWorldPosition();
        skullNode->SetDeepEnabled(false);
        ballNode->SetDeepEnabled(true);
        ballNode->SetWorldPosition(position + Vector3::UP);
    }

    // Replace environment
    Node* envLightNode = scene_->GetChild("EnvironmentLight", true);
    Node* envDarkNode = scene_->GetChild("EnvironmentDark", true);
    if (dark)
    {
        envLightNode->SetDeepEnabled(false);
        envDarkNode->SetDeepEnabled(true);
    }
    else
    {
        envLightNode->SetDeepEnabled(true);
        envDarkNode->SetDeepEnabled(false);
    }

    // Replace car
    Node* carNode = scene_->GetChild("Car", true);
    Node* ghostCarNode = carNode->GetChild("GhostCar", true);
    ghostCarNode->SetDeepEnabled(dark);
    auto carModel = carNode->GetComponent<StaticModel>();
    carModel->SetMaterial(0, dark
        ? cache->GetResource<Material>("ForestScene/Car/Body_mat2.xml")
        : cache->GetResource<Material>("ForestScene/Car/Body_mat.xml"));

    // Replace grass
    Node* grassRegionNode = scene_->GetChild("GrassRegion", true);
    for (Node* grassNode : grassRegionNode->GetChildren())
    {
        auto grassModel = grassNode->GetComponent<StaticModel>();
        grassModel->SetMaterial(0, dark ? grassDarkMaterial_ : grassLightMaterial_);
    }

    // Replace terrain
    Terrain* terrain = scene_->GetComponent<Terrain>(true);
    terrain->SetMaterial(dark
        ? cache->GetResource<Material>("ForestScene/Terrain/Terrain_mat2.xml")
        : cache->GetResource<Material>("ForestScene/Terrain/Terrain_mat.xml"));

    // Replace trees
    for (Node* treeNode : scene_->GetChild("Trees")->GetChildren())
    {
        auto treeModel = treeNode->GetComponent<StaticModel>();
        treeModel->SetMaterial(0, dark
            ? cache->GetResource<Material>("ForestScene/FirTree/Leaf_mat2.xml")
            : cache->GetResource<Material>("ForestScene/FirTree/Leaf_mat.xml"));
    }

    // Finalize
    darkTheme_ = dark;
}

void CharacterDemo::UpdateGrassTexture(float timeStep)
{
    const float MAX_WIND_PUSHINESS = 0.1f;
    const float WIND_SCALE_T = 0.17f;
    const float WIND_SCALE_X = 0.13f;
    const float WIND_SCALE_Z = 0.07f;
    const float GRASS_RESTORATION_SPEED = 0.2f;
    const float time = scene_->GetElapsedTime();
    // Reset grass texture
    for (unsigned y = 0; y < grassTextureSize_; ++y)
    {
        for (unsigned x = 0; x < grassTextureSize_; ++x)
        {
            const unsigned i = y * grassTextureSize_ + x;
            const float restorationSpeed = GRASS_RESTORATION_SPEED * Pow(Max(0.1f, 1.0f - grassPushiness_[i]),
                darkTheme_ ? 1.3f : 2.5f);
            grassPushiness_[i] = Max(0.0f, grassPushiness_[i] - restorationSpeed * timeStep);
        }
    }

    // Render physical objects
    auto v2i = [&](const Vector3& v)
    {
        IntVector2 i;
        i.x_ = Clamp((v.x_ - grassBoundingBox_.min_.x_) / grassBoundingBox_.Size().x_, 0.0f, 1.0f) * (grassTextureSize_ - 1);
        i.y_ = Clamp((v.z_ - grassBoundingBox_.min_.z_) / grassBoundingBox_.Size().z_, 0.0f, 1.0f) * (grassTextureSize_ - 1);
        return i;
    };
    auto i2v = [&](const IntVector2& i)
    {
        Vector3 v;
        v.x_ = (i.x_ + 0.5f) / (grassTextureSize_ - 1) * grassBoundingBox_.Size().x_ + grassBoundingBox_.min_.x_;
        v.y_ = grassBaseHeight_[i.y_ * grassTextureSize_ + i.x_];
        v.z_ = (i.y_ + 0.5f) / (grassTextureSize_ - 1) * grassBoundingBox_.Size().z_ + grassBoundingBox_.min_.z_;
        return v;
    };
    auto physicsWorld = scene_->GetComponent<PhysicsWorld>();
    PODVector<PhysicsRaycastResult> raycastResults;
    auto v4b = [&](const Vector3& v, const Matrix3x4& invTransform)
    {
        const Vector3 vt = invTransform * v;
        if (Abs(vt.x_) < 1.0 && Abs(vt.y_) < 1.0 && Abs(vt.z_) < 1.0)
            return 1.0f;

        raycastResults.Clear();

        physicsWorld->Raycast(raycastResults, Ray(v, Vector3::UP), GRASS_BILLBOARD_SIZE, LAYER_GRASSPUSHER);

        if (!raycastResults.Empty())
            return GRASS_BILLBOARD_SIZE - raycastResults[0].distance_;

        return 0.0f;
    };
    auto v2s = [&](const Vector3& v, const Matrix3x4& invTransform)
    {
        const Vector3 vt = invTransform * v;
        if (vt.Length() < 1.0f)
            return 1.0f;

        raycastResults.Clear();

        physicsWorld->Raycast(raycastResults, Ray(v, Vector3::UP), GRASS_BILLBOARD_SIZE, LAYER_GRASSPUSHER);

        if (!raycastResults.Empty())
            return GRASS_BILLBOARD_SIZE - raycastResults[0].distance_;

        return 0.0f;
    };
    const float texelSize = grassBoundingBox_.Size().x_ / grassTextureSize_;
    for (CollisionShape* shape : grassPushers_)
    {
        Node* node = shape->GetNode();
        const Matrix3x4 invTransform = node->GetWorldTransform().Inverse();
        const BoundingBox bbox = shape->GetWorldBoundingBox();
        const IntVector2 fromXY = v2i(bbox.min_ - 2 * texelSize * Vector3::ONE);
        const IntVector2 toXY = v2i(bbox.max_ + 2 * texelSize * Vector3::ONE);
        IntVector2 xy;
        switch (shape->GetShapeType())
        {
        case SHAPE_BOX:
            for (xy.y_ = fromXY.y_; xy.y_ <= toXY.y_; ++xy.y_)
            {
                for (xy.x_ = fromXY.x_; xy.x_ <= toXY.x_; ++xy.x_)
                {
                    const unsigned i = xy.y_ * grassTextureSize_ + xy.x_;
                    grassPushiness_[i] = Max(grassPushiness_[i], v4b(i2v(xy), invTransform));
                }
            }
            break;
        case SHAPE_SPHERE:
            for (xy.y_ = fromXY.y_; xy.y_ <= toXY.y_; ++xy.y_)
            {
                for (xy.x_ = fromXY.x_; xy.x_ <= toXY.x_; ++xy.x_)
                {
                    const unsigned i = xy.y_ * grassTextureSize_ + xy.x_;
                    const float radius = Max(texelSize * 0.5f, node->GetWorldScale().x_ * shape->GetSize().x_);
                    const Vector3 center = node->GetWorldPosition();
                    const Vector3 grassPosition = i2v(xy);
                    const float supression = Min(1.0f, (Vector3(1, 0, 1) * (center - grassPosition)).Length() / radius);
                    const float dh = Max(0.0f, center.y_ - grassPosition.y_ - radius);
                    grassPushiness_[i] = Max(grassPushiness_[i], Max(0.0f, Sqrt(1.0f - supression * supression) - dh));
                }
            }
            break;
        default:
            break;
        }
    }

    // Apply wind
    for (unsigned y = 0; y < grassTextureSize_; ++y)
    {
        for (unsigned x = 0; x < grassTextureSize_; ++x)
        {
            const unsigned i = y * grassTextureSize_ + x;
            const float kx = (x + 0.5f) / grassTextureSize_ * grassBoundingBox_.Size().x_;
            const float ky = (y + 0.5f) / grassTextureSize_ * grassBoundingBox_.Size().z_;
            const float wavePhase = (SmoothTriangleWave(kx * WIND_SCALE_X) + SmoothTriangleWave(ky * WIND_SCALE_Z)) / 2;
            const float minPushiness = SmoothTriangleWave(time * WIND_SCALE_T + wavePhase);
            grassPushiness_[i] = Min(maxGrassPushiness_, Max(MAX_WIND_PUSHINESS * minPushiness, grassPushiness_[i]));
        }
    }
    // Fill texture
    auto getPushiness = [&](unsigned x, unsigned y)
    {
        return grassPushiness_[y * grassTextureSize_ + x];
    };
    const float deltaScale = 0.5f * GRASS_BILLBOARD_SIZE / texelSize;
    for (unsigned y = 0; y < grassTextureSize_; ++y)
    {
        for (unsigned x = 0; x < grassTextureSize_; ++x)
        {
            const unsigned i = y * grassTextureSize_ + x;
#if 0
            const float hx0 = x > 0 ? getPushiness(x - 1, y) : 0.0f;
            const float hx1 = x + 1 < grassTextureSize_ ? getPushiness(x + 1, y) : 0.0f;
            const float hy0 = y > 0 ? getPushiness(x, y - 1) : 0.0f;
            const float hy1 = y + 1 < grassTextureSize_ ? getPushiness(x, y + 1) : 0.0f;
            const float dhx = deltaScale * (hx1 - hx0);
            const float dhy = deltaScale * (hy1 - hy0);
            grassTextureData_[i * 4 + 0] = static_cast<unsigned char>(Clamp((dhx * 0.5f + 0.5f) * 255, 0.0f, 255.0f));
            grassTextureData_[i * 4 + 1] = static_cast<unsigned char>(Clamp((dhy * 0.5f + 0.5f) * 255, 0.0f, 255.0f));
#else
            grassTextureData_[i * 4 + 0] = 0;
            grassTextureData_[i * 4 + 1] = 0;
#endif
            grassTextureData_[i * 4 + 2] = static_cast<unsigned char>(Clamp(grassScale_[i] * 255, 0.0f, 255.0f));
            grassTextureData_[i * 4 + 3] = static_cast<unsigned char>(Clamp(grassPushiness_[i] * 255, 0.0f, 255.0f));
        }
    }
    grassTexture_->SetData(0, 0, 0, grassTextureSize_, grassTextureSize_, grassTextureData_.Buffer());
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
    body->SetCollisionLayer(LAYER_MAIN);
    body->SetCollisionMask(LAYER_OBSTACLE);
    body->SetMass(1.0f);

    // Set zero angular factor so that physics doesn't turn the character on its own.
    // Instead we will control the character yaw manually
    body->SetAngularFactor(Vector3::ZERO);

    // Set the rigidbody to signal collision also when in rest, so that we get ground collisions properly
    body->SetCollisionEventMode(COLLISION_ALWAYS);

    const float KICK_VELOCITY = 10.0f;

    SubscribeToEvent(objectNode, E_NODECOLLISION,
        [=](StringHash eventType, VariantMap& eventData)
    {
        Node* otherNode = static_cast<Node*>(eventData[NodeCollision::P_OTHERNODE].GetPtr());

        if (kickCooldownTimer_ <= 0.0f && otherNode->HasTag("KICKABLE"))
        {
            kickCooldownTimer_ = kickCooldown_;
            RigidBody* otherBody = static_cast<RigidBody*>(eventData[NodeCollision::P_OTHERBODY].GetPtr());
            const Vector3 kickDirection = (cameraNode_->GetWorldDirection() * Vector3(1, 0, 1)).Normalized() + Vector3::UP * 0.3f;
            otherBody->ApplyImpulse(kickDirection.Normalized() * otherBody->GetMass() * KICK_VELOCITY);
            otherBody->SetAngularVelocity(cameraNode_->GetRight() * 10);
        }
    });

    // Set a capsule shape for collision
    auto* shape = objectNode->CreateComponent<CollisionShape>();
    shape->SetCapsule(0.7f, 1.8f, Vector3(0.0f, 0.9f, 0.0f));

    // Create the character logic component, which takes care of steering the rigidbody
    // Remember it so that we can set the controls. Use a WeakPtr because the scene hierarchy already owns it
    // and keeps it alive as long as it's not removed from the hierarchy
    character_ = objectNode->CreateComponent<Character>();

    // Create grass pusher
    Node* grassPusher = objectNode->CreateChild("CharacterGrassPusher");
    {
        auto* body = objectNode->CreateComponent<RigidBody>();
        body->SetCollisionLayer(LAYER_GRASSPUSHER);
        body->SetCollisionMask(0);

        auto* shape = objectNode->CreateComponent<CollisionShape>();
        shape->SetSphere(0.4f, Vector3(0.0f, 0.35f, 0.0f));
    }
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

    const float timeStep = eventData[Update::P_TIMESTEP].GetFloat();
    auto* input = GetSubsystem<Input>();

    kickCooldownTimer_ -= timeStep;

    if (input->GetKeyPress(KEY_TAB))
    {
        SetupTheme(!darkTheme_);
    }

    // Update grass
    if (grassTextureSize_ > 0)
    {
        UpdateGrassTexture(timeStep);
    }

    // Update skull light
//     Node* skullNode = scene_->GetChild("Skull");
//     Node* skullLightNode = skullNode->GetChild("SkullLight");
//     skullLightNode->SetWorldPosition(skullNode->GetWorldPosition() + Vector3::UP * 1.0);

    Node* handGrassPusher = scene_->GetChild("HandGrassPusher", true);
    if (input->GetMouseButtonDown(MOUSEB_LEFT))
    {
        auto physicsWorld = scene_->GetComponent<PhysicsWorld>();
        const Vector3 camPosition = cameraNode_->GetWorldPosition();
        const Vector3 camDirection = cameraNode_->GetDirection();
        const float maxDistance = 5.0f;
        PhysicsRaycastResult raycast;
        physicsWorld->RaycastSingle(raycast, Ray(camPosition, camDirection), maxDistance, LAYER_MAIN);
        const float distance = raycast.body_ ? Min(maxDistance, raycast.distance_) : maxDistance;
        handGrassPusher->SetWorldPosition(camPosition + camDirection * distance);
    }
    else
    {
        handGrassPusher->SetWorldPosition(Vector3(0, 100, 0));
    }

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
