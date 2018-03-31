#include "Constants.hlsl"
#include "Uniforms.hlsl"
#include "Samplers.hlsl"
#include "Transform.hlsl"
#include "Lighting.hlsl"
#include "ScreenPos.hlsl"
#include "Fog.hlsl"

#if defined(COMPILEPS) && defined(SOFTPARTICLES)
#ifndef D3D11
// D3D9 uniform
uniform float cSoftParticleFadeScale;
#else
// D3D11 constant buffer
cbuffer CustomPS : register(b6)
{
    float cSoftParticleFadeScale;
}
#endif
#endif

#ifdef COMPILEVS

void VS(float4 iPos : POSITION,
    #if !defined(BILLBOARD) && !defined(TRAILFACECAM) && !defined(POINTBILLBOARD)
        float3 iNormal : NORMAL,
    #endif
    #ifndef NOUV
        #ifdef POINTEXPAND
            float4 iTexCoord : TEXCOORD0,
        #else
            float2 iTexCoord : TEXCOORD0,
        #endif
    #endif
    #ifdef VERTEXCOLOR
        float4 iColor : COLOR0,
    #endif
    #if defined(BILLBOARD) || defined(DIRBILLBOARD)
        float2 iSize : TEXCOORD1,
    #endif
    #if defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD)
        float4 iSize : TEXCOORD1,
    #endif
    #if defined(TRAILFACECAM) || defined(TRAILBONE)
        float4 iTangent : TANGENT,
    #endif
    #if defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD)
        out float4 oTexCoord : TEXCOORD0,
        out float4 oSize : TEXCOORD1,
        #if defined(POINTDIRBILLBOARD)
            out float3 oNormal : NORMAL,
        #endif
    #else
        out float2 oTexCoord : TEXCOORD0,
    #endif
    #if defined(SOFTPARTICLES) && !defined(POINTEXPAND)
        out float4 oScreenPos : TEXCOORD1,
    #endif
    out float4 oWorldPos : TEXCOORD3,
    #if defined(PERPIXEL) && !defined(POINTEXPAND)
        #ifdef SHADOW
            out float4 oShadowPos[NUMCASCADES] : TEXCOORD4,
        #endif
        #ifdef SPOTLIGHT
            out float4 oSpotPos : TEXCOORD5,
        #endif
        #ifdef POINTLIGHT
            out float3 oCubeMaskVec : TEXCOORD5,
        #endif
    #elif !defined(POINTEXPAND)
        out float3 oVertexLight : TEXCOORD4,
    #endif
    #ifdef VERTEXCOLOR
        out float4 oColor : COLOR0,
    #endif
    #if defined(D3D11) && defined(CLIPPLANE)
        out float oClip : SV_CLIPDISTANCE0,
    #endif
    out float4 oPos : OUTPOSITION)
{
    // Define a 0,0 UV coord if not expected from the vertex data
    #ifdef NOUV
    float2 iTexCoord = float2(0.0, 0.0);
    #endif
    
    #if defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD)
        float3 worldPos = iPos.xyz;
        oPos = GetClipPos(worldPos);
        oTexCoord = iTexCoord;
        oSize = iSize;
        oWorldPos = float4(worldPos, GetDepth(oPos));
        #ifdef POINTDIRBILLBOARD
            oNormal = iNormal;
        #endif
    #else
        float4x3 modelMatrix = iModelMatrix;
        float3 worldPos = GetWorldPos(modelMatrix);
        oPos = GetClipPos(worldPos);
        oTexCoord = GetTexCoord(iTexCoord);
        oWorldPos = float4(worldPos, GetDepth(oPos));
    #endif

    #if defined(D3D11) && defined(CLIPPLANE)
        oClip = dot(oPos, cClipPlane);
    #endif

    #if defined(SOFTPARTICLES) && !defined(POINTEXPAND)
        oScreenPos = GetScreenPos(oPos);
    #endif

    #ifdef VERTEXCOLOR
        oColor = iColor;
    #endif

    #if defined(PERPIXEL) && !(defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD))
        // Per-pixel forward lighting
        float4 projWorldPos = float4(worldPos.xyz, 1.0);

        #ifdef SHADOW
            // Shadow projection: transform from world space to shadow space
            GetShadowPos(projWorldPos, float3(0, 0, 0), oShadowPos);
        #endif

        #ifdef SPOTLIGHT
            // Spotlight projection: transform from world space to projector texture coordinates
            oSpotPos = mul(projWorldPos, cLightMatrices[0]);
        #endif

        #ifdef POINTLIGHT
            oCubeMaskVec = mul(worldPos - cLightPos.xyz, (float3x3)cLightMatrices[0]);
        #endif
    #elif !defined(POINTEXPAND)
        // Ambient & per-vertex lighting
        oVertexLight = GetAmbient(GetZonePos(worldPos));

        #ifdef NUMVERTEXLIGHTS
            for (int i = 0; i < NUMVERTEXLIGHTS; ++i)
                oVertexLight += GetVertexLightVolumetric(i, worldPos) * cVertexLights[i * 3].rgb;
        #endif
    #endif
}

#endif


#ifdef COMPILEGS

struct VSOutput
{
    float4 iTexCoord : TEXCOORD0;
    float4 iSize : TEXCOORD1;
    #ifdef POINTDIRBILLBOARD
        float3 iNormal : NORMAL;
    #endif
    float4 iWorldPos: TEXCOORD3;
    #ifdef VERTEXCOLOR
        float4 iColor : COLOR0;
    #endif
};

struct GSOutput
{
    float2 oTexCoord : TEXCOORD0;
    #ifdef SOFTPARTICLES
        float4 oScreenPos: TEXCOORD1;
    #endif
    float4 oWorldPos: TEXCOORD3;
    
    #ifdef PERPIXEL
        #ifdef SHADOW
            float4 oShadowPos[NUMCASCADES] : TEXCOORD4;
        #endif
        #ifdef SPOTLIGHT
            float4 oSpotPos : TEXCOORD5;
        #endif
        #ifdef POINTLIGHT
            float3 oCubeMaskVec : TEXCOORD5;
        #endif
    #else
        float3 oVertexLight : TEXCOORD4;
    #endif
    
    #ifdef VERTEXCOLOR
        float4 oColor : COLOR0;
    #endif
    
    float4 oPos : SV_POSITION;
};

[maxvertexcount(4)]
void GS(point in VSOutput vertexData[1], inout TriangleStream<GSOutput> triStream)
{
    VSOutput vertData = vertexData[0];

    float2 minUV = vertData.iTexCoord.xy;
    float2 maxUV = vertData.iTexCoord.zw;
    float2 partSize = vertData.iSize.xy;

    float2 partUV[] = {
        float2(minUV.x, maxUV.y),
        minUV,
        maxUV,
        float2(maxUV.x, minUV.y),
    };

    float s;
    float c;
    sincos(vertData.iSize.z * M_DEGTORAD, s, c);
    float3 vUpNew    = c * float3(1, 0, 0) - s * float3(0, 1, 0);
    float3 vRightNew = s * float3(1, 0, 0) + c * float3(0, 1, 0);
    vUpNew *= partSize.y;
    vRightNew *= partSize.x;
    
    float3 partOffset[] = {
        -vUpNew + -vRightNew,
        -vUpNew + vRightNew,
        vUpNew + -vRightNew,
        vUpNew + vRightNew,
    };
    
    float3 worldLoc = mul(vertData.iWorldPos.xyz, iModelMatrix);
    
    [unroll]
    for (int i = 0; i < 4; ++i)
    {
        GSOutput outVert = (GSOutput)0;
        
        #ifdef POINTBILLBOARD       
            float3 vPos = mul(partOffset[i], cBillboardRot);
        #endif
        #ifdef POINTDIRBILLBOARD
            float3 vPos = mul(partOffset[i].xzy, GetFaceCameraRotation(worldLoc, vertData.iNormal));
        #endif
        
        outVert.oWorldPos = float4(worldLoc + vPos, 0);
        outVert.oPos = GetClipPos(outVert.oWorldPos);
        outVert.oTexCoord = GetTexCoord(partUV[i]);
        outVert.oColor = vertData.iColor;
        
        #ifdef SOFTPARTICLES
            outVert.oScreenPos = GetScreenPos(outVert.oPos);
        #endif
        
        #ifdef PERPIXEL
            // Per-pixel forward lighting
            float4 projWorldPos = float4(outVert.oPos.xyz, 1.0);

            #ifdef SHADOW
                // Shadow projection: transform from world space to shadow space
                GetShadowPos(projWorldPos, float3(0, 0, 0), outVert.oShadowPos);
            #endif

            #ifdef SPOTLIGHT
                // Spotlight projection: transform from world space to projector texture coordinates
                outVert.oSpotPos = mul(projWorldPos, cLightMatrices[0]);
            #endif

            #ifdef POINTLIGHT
                outVert.oCubeMaskVec = mul(projWorldPos - cLightPos.xyz, (float3x3)cLightMatrices[0]);
            #endif
        #else
            // Ambient & per-vertex lighting
            outVert.oVertexLight = GetAmbient(GetZonePos(outVert.oWorldPos));

            #ifdef NUMVERTEXLIGHTS
                for (int i = 0; i < NUMVERTEXLIGHTS; ++i)
                    outVert.oVertexLight += GetVertexLightVolumetric(i, outVert.oWorldPos) * cVertexLights[i * 3].rgb;
            #endif
        #endif
        
        triStream.Append(outVert);
    }
}


#endif

#ifdef COMPILEPS

void PS(float2 iTexCoord : TEXCOORD0,
    #ifdef SOFTPARTICLES
        float4 iScreenPos: TEXCOORD1,
    #endif
    float4 iWorldPos : TEXCOORD3,
    #ifdef PERPIXEL
        #ifdef SHADOW
            float4 iShadowPos[NUMCASCADES] : TEXCOORD4,
        #endif
        #ifdef SPOTLIGHT
            float4 iSpotPos : TEXCOORD5,
        #endif
        #ifdef POINTLIGHT
            float3 iCubeMaskVec : TEXCOORD5,
        #endif
    #else
        float3 iVertexLight : TEXCOORD4,
    #endif
    #ifdef VERTEXCOLOR
        float4 iColor : COLOR0,
    #endif
    #if defined(D3D11) && defined(CLIPPLANE)
        float iClip : SV_CLIPDISTANCE0,
    #endif
    out float4 oColor : OUTCOLOR0)
{
    // Get material diffuse albedo
    #ifdef DIFFMAP
        float4 diffInput = Sample2D(DiffMap, iTexCoord);
        #ifdef ALPHAMASK
            if (diffInput.a < 0.5)
                discard;
        #endif
        float4 diffColor = cMatDiffColor * diffInput;
    #else
        float4 diffColor = cMatDiffColor;
    #endif

    #ifdef VERTEXCOLOR
        diffColor *= iColor;
    #endif

    // Get fog factor
    #ifdef HEIGHTFOG
        float fogFactor = GetHeightFogFactor(iWorldPos.w, iWorldPos.y);
    #else
        float fogFactor = GetFogFactor(iWorldPos.w);
    #endif

    // Soft particle fade
    // In expand mode depth test should be off. In that case do manual alpha discard test first to reduce fill rate
    #ifdef SOFTPARTICLES
        #ifdef EXPAND
            if (diffColor.a < 0.01)
                discard;
        #endif

        float particleDepth = iWorldPos.w;
        float depth = Sample2DProj(DepthBuffer, iScreenPos).r;
        #ifdef HWDEPTH
            depth = ReconstructDepth(depth);
        #endif

        #ifdef EXPAND
            float diffZ = max(particleDepth - depth, 0.0) * (cFarClipPS - cNearClipPS);
            float fade = saturate(diffZ * cSoftParticleFadeScale);
        #else
            float diffZ = (depth - particleDepth) * (cFarClipPS - cNearClipPS);
            float fade = saturate(1.0 - diffZ * cSoftParticleFadeScale);
        #endif

        diffColor.a = max(diffColor.a - fade, 0.0);
    #endif

    #ifdef PERPIXEL
        // Per-pixel forward lighting
        float3 lightColor;
        float3 finalColor;
        
        float diff = GetDiffuseVolumetric(iWorldPos.xyz);

        #ifdef SHADOW
            diff *= GetShadow(iShadowPos, iWorldPos.w);
        #endif

        #if defined(SPOTLIGHT)
            lightColor = iSpotPos.w > 0.0 ? Sample2DProj(LightSpotMap, iSpotPos).rgb * cLightColor.rgb : 0.0;
        #elif defined(CUBEMASK)
            lightColor = texCUBE(sLightCubeMap, iCubeMaskVec).rgb * cLightColor.rgb;
        #else
            lightColor = cLightColor.rgb;
        #endif

        finalColor = diff * lightColor * diffColor.rgb;
        oColor = float4(GetLitFog(finalColor, fogFactor), diffColor.a);
    #else
        // Ambient & per-vertex lighting
        float3 finalColor = iVertexLight * diffColor.rgb;

        oColor = float4(GetFog(finalColor, fogFactor), diffColor.a);
    #endif
}

#endif
