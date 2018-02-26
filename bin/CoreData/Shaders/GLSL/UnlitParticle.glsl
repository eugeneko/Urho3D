#include "Uniforms.glsl"
#include "Samplers.glsl"
#include "Transform.glsl"
#include "ScreenPos.glsl"
#include "Fog.glsl"

#ifdef SOFTPARTICLES
    uniform float cSoftParticleFadeScale;
#endif

#ifndef GL3

    varying vec2 vTexCoord;
    varying vec4 vWorldPos;
    #ifdef VERTEXCOLOR
        varying vec4 vColor;
    #endif
    #ifdef SOFTPARTICLES
        varying vec4 vScreenPos;
    #endif

#else

#ifdef COMPILEGS
    #define VSDATA_TAIL []
#else
    #define VSDATA_TAIL
#endif

#ifdef COMPILEVS
    out 
#else
    in
#endif
    VertexData
    {
        #if (defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD)) && !defined(COMPILEPS)
            vec4 vTexCoord;
            vec4 vSize;
        #else
            vec2 vTexCoord;
        #endif
        #if defined(POINTDIRBILLBOARD) && !defined(COMPILEPS)
            vec3 vNormal;
        #endif
        vec4 vWorldPos;
        #ifdef VERTEXCOLOR
            vec4 vColor;
        #endif
        #ifdef SOFTPARTICLES
            vec4 vScreenPos;
        #endif
    } vs_out VSDATA_TAIL;
    
#endif

#ifdef COMPILEVS

    #define vTexCoord vs_out.vTexCoord
    #define vWorldPos vs_out.vWorldPos
    #define vColor vs_out.vColor
    #define vScreenPos vs_out.vScreenPos
    #define vSize vs_out.vSize
    #define vNormal vs_out.vNormal

void VS()
{
    #if defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD)
        mat4 modelMatrix = iModelMatrix;
        vec3 worldPos = GetWorldPos(modelMatrix);
        gl_Position = GetClipPos(worldPos);
        vTexCoord = iTexCoord;
        vSize = iTexCoord1;
        vWorldPos = vec4(worldPos, GetDepth(gl_Position));
        #if defined(POINTDIRBILLBOARD)
            vNormal = iNormal;
        #endif
    #else
        mat4 modelMatrix = iModelMatrix;
        vec3 worldPos = GetWorldPos(modelMatrix);
        gl_Position = GetClipPos(worldPos);
        vTexCoord = GetTexCoord(iTexCoord);
        vWorldPos = vec4(worldPos, GetDepth(gl_Position));
    #endif

    #ifdef SOFTPARTICLES
        vScreenPos = GetScreenPos(gl_Position);
    #endif

    #ifdef VERTEXCOLOR
        vColor = iColor;
    #endif

}

#endif

#if defined(COMPILEGS) && (defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD))
#line 100
    
    out VertexData
    {
        vec2 vTexCoord;
        vec4 vWorldPos;
        #ifdef POINTDIRBILLBOARD
            vec3 vNormal;
        #endif
        #ifdef VERTEXCOLOR
            vec4 vColor;
        #endif
        #ifdef SOFTPARTICLES
            vec4 vScreenPos;
        #endif
    } gs_out;
    
    #define oTexCoord gs_out.vTexCoord
    #define oWorldPos gs_out.vWorldPos
    #define oColor gs_out.vColor
    #define oScreenPos gs_out.vScreenPos

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

mat3 GetFaceCameraRotation(vec3 position, vec3 direction)
{
    vec3 cameraDir = normalize(position - cCameraPos);
    vec3 front = normalize(direction);
    vec3 right = normalize(cross(front, cameraDir));
    vec3 up = normalize(cross(front, right));

    return mat3(
        right.x, up.x, front.x,
        right.y, up.y, front.y,
        right.z, up.z, front.z
    );
}

void GS()
{
    vec2 minUV = vs_out[0].vTexCoord.xy;
    vec2 maxUV = vs_out[0].vTexCoord.zw;
    vec2 minSize = vs_out[0].vSize.xy;
    vec2 maxSize = vs_out[0].vSize.zw;
    
    vec2 partUV[] = {
        minUV,
        vec2(minUV.x, maxUV.y),
        vec2(maxUV.x, minUV.y),
        maxUV,
    };
    
    vec3 partOffset[] = {
        vec3(minSize.x, minSize.y, 0),
        vec3(minSize.x, maxSize.y, 0),
        vec3(maxSize.x, minSize.y, 0),
        vec3(maxSize.x, maxSize.y, 0),
    };

    vec4 wPos = vs_out[0].vWorldPos;
    for (int i = 0; i < 4; ++i)
    {        
        #if defined(POINTBILLBOARD)
            vec3 vPos = partOffset[i] * cBillboardRot;
        #elif defined(POINTDIRBILLBOARD)
            vec3 vPos = vec3(partOffset[i].x, 0, partOffset[i].y) * GetFaceCameraRotation(wPos.xyz, vs_out[0].vNormal);
        #endif
        
        oWorldPos = vec4((wPos + vPos).xyz, 0);
        oTexCoord = partUV[i];
        oColor = vs_out[0].vColor;
        gl_Position = vec4(oWorldPos.xyz, 1.0) * cViewProj;
        EmitVertex();
    }
    EndPrimitive();
}

#endif

#ifdef COMPILEPS

    #define vTexCoord vs_out.vTexCoord
    #define vWorldPos vs_out.vWorldPos
    #define vColor vs_out.vColor
    #define vScreenPos vs_out.vScreenPos

void PS()
{
    // Get material diffuse albedo
    #ifdef DIFFMAP
        vec4 diffColor = cMatDiffColor * texture2D(sDiffMap, vTexCoord);
        #ifdef ALPHAMASK
            if (diffColor.a < 0.5)
                discard;
        #endif
    #else
        vec4 diffColor = cMatDiffColor;
    #endif

    #ifdef VERTEXCOLOR
        diffColor *= vColor;
    #endif

    // Get fog factor
    #ifdef HEIGHTFOG
        float fogFactor = GetHeightFogFactor(vWorldPos.w, vWorldPos.y);
    #else
        float fogFactor = GetFogFactor(vWorldPos.w);
    #endif

    // Soft particle fade
    // In expand mode depth test should be off. In that case do manual alpha discard test first to reduce fill rate
    #ifdef SOFTPARTICLES
        #ifdef EXPAND
            if (diffColor.a < 0.01)
                discard;
        #endif

        float particleDepth = vWorldPos.w;
        #ifdef HWDEPTH
            float depth = ReconstructDepth(texture2DProj(sDepthBuffer, vScreenPos).r);
        #else
            float depth = DecodeDepth(texture2DProj(sDepthBuffer, vScreenPos).rgb);
        #endif

        #ifdef EXPAND
            float diffZ = max(particleDepth - depth, 0.0) * (cFarClipPS - cNearClipPS);
            float fade = clamp(diffZ * cSoftParticleFadeScale, 0.0, 1.0);
        #else
            float diffZ = (depth - particleDepth) * (cFarClipPS - cNearClipPS);
            float fade = clamp(1.0 - diffZ * cSoftParticleFadeScale, 0.0, 1.0);
        #endif

        #ifndef ADDITIVE
            diffColor.a = max(diffColor.a - fade, 0.0);
        #else
            diffColor.rgb = max(diffColor.rgb - fade, vec3(0.0, 0.0, 0.0));
        #endif
    #endif

    gl_FragColor = vec4(GetFog(diffColor.rgb, fogFactor), diffColor.a);
}

#endif
