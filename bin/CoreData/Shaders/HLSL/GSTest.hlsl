#include "Uniforms.hlsl"
#include "Samplers.hlsl"
#include "Transform.hlsl"
#include "ScreenPos.hlsl"

struct VS_INPUT
{
    float4 iPos : POSITION;
};

struct VS_OUTPUT
{
    float4 oPos : OUTPOSITION;
    float4 oWorldPos : TEXCOORD0;
};

struct GS_OUTPUT
{
    float4 oPos : SV_POSITION;
    float4 oColor : TEXCOORD0;
};

#ifdef COMPILEVS
VS_OUTPUT VS(in VS_INPUT input)
{
    VS_OUTPUT outData = (VS_OUTPUT)0;

    float4x3 modelMatrix = cModel;
    float3 worldPos = mul(input.iPos, cModel);
    outData.oPos = GetClipPos(worldPos);
    outData.oWorldPos = float4(worldPos, GetDepth(outData.oPos));
    //float3 worldPos = GetWorldPos(modelMatrix);
    //outData.oPos = GetClipPos(worldPos);
    //outData.oWorldPos = float4(worldPos, GetDepth(oPos));
    //outData.oPos = GetClipPos(input.iPos);
    //outData.oPos = worldPos;
    //outData.oPos = input.iPos;

    //float4x3 modelMatrix = iModelMatrix;
    //float3 worldPos = GetWorldPos(modelMatrix);
    //outData.oPos = GetClipPos(worldPos);

    return outData;
}
#endif

void CreateVertex(inout TriangleStream<GS_OUTPUT> triStream, float4 pos, float4 col)
{
    GS_OUTPUT temp = (GS_OUTPUT)0;
    temp.oPos = pos;
    temp.oColor = col;
    triStream.Append(temp);
}

#if defined(COMPILEGS)
[maxvertexcount(3 * 3)]
void GS(triangle in VS_OUTPUT vertices[3], inout TriangleStream<GS_OUTPUT> triStream)
{
    const float gLayers = 3.0;

    VS_OUTPUT v1 = vertices[0],
        v2 = vertices[1],
        v3 = vertices[2];
    float length = 0.1;

    float3 offsets[] = {
        float3(3, 3, 0),
        float3(-3, -3, 0),
        float3(3, 0, -3)
    };
    float4 colors[] = {
        float4(0.8, 0.6, 0.3, 0.5),
        float4(0.2, 0.9, 0.5, 0.5),
        float4(0.1, 0.33, 0.7, 0.5),
    };

    float offset = length / gLayers;
    CreateVertex(triStream, v1.oPos, colors[0]);
    CreateVertex(triStream, v2.oPos, colors[0]);
    CreateVertex(triStream, v3.oPos, colors[0]);
    triStream.RestartStrip();

    for (int i = 0; i < 3; ++i)
    {
        float3 norm = offsets[i];
        CreateVertex(triStream, v1.oPos + float4(norm*0.03, 1) + v1.oPos*0.1, colors[i]);
        CreateVertex(triStream, v2.oPos, colors[i]);
        CreateVertex(triStream, v3.oPos + float4(norm*0.06, 1) + v1.oPos*0.1, colors[i]);
        triStream.RestartStrip();
    }
}
#endif

float4 PS(in GS_OUTPUT gsData) : SV_TARGET
{
    return gsData.oColor;
}

