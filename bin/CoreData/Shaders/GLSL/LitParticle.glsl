#include "Constants.glsl"
#include "Uniforms.glsl"
#include "Samplers.glsl"
#include "Transform.glsl"
#include "ScreenPos.glsl"
#include "Lighting.glsl"
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
        uniform float cSoftParticleFadeScale;
    #endif
    #ifdef PERPIXEL
        #ifdef SHADOW
            #ifndef GL_ES
                varying vec4 vShadowPos[NUMCASCADES];
            #else
                varying highp vec4 vShadowPos[NUMCASCADES];
            #endif
        #endif
        #ifdef SPOTLIGHT
            varying vec4 vSpotPos;
        #endif
        #ifdef POINTLIGHT
            varying vec3 vCubeMaskVec;
        #endif
    #else
        varying vec3 vVertexLight;
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
        #if (defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD)) && defined(POINTEXPAND) && !defined(COMPILEPS)
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
        #if defined(SOFTPARTICLES) && (defined(COMPILEPS) || !defined(POINTEXPAND))
            vec4 vScreenPos;
        #endif

        #if defined(COMPILEPS) || (defined(COMPILEVS) && !defined(POINTEXPAND))
            #ifdef PERPIXEL
                #ifdef SHADOW
                    #ifndef GL_ES
                        vec4 vShadowPos[NUMCASCADES];
                    #else
                        highp vec4 vShadowPos[NUMCASCADES];
                    #endif
                #endif
                #ifdef SPOTLIGHT
                    vec4 vSpotPos;
                #endif
                #ifdef POINTLIGHT
                    vec3 vCubeMaskVec;
                #endif
            #else
                vec3 vVertexLight;
            #endif
        #endif
        
    } vs_out VSDATA_TAIL;

#endif

#ifdef COMPILEVS

#ifdef GL3
    #define vTexCoord vs_out.vTexCoord
    #define vWorldPos vs_out.vWorldPos
    #define vColor vs_out.vColor
    #define vScreenPos vs_out.vScreenPos
    #define vSize vs_out.vSize
    #define vNormal vs_out.vNormal
    #define vCubeMaskVec vs_out.vCubeMaskVec
    #define vSpotPos vs_out.vSpotPos
    #define vShadowPos vs_out.vShadowPos
    #define vVertexLight vs_out.vVertexLight
#endif

void VS()
{
    #if defined(POINTBILLBOARD) || defined(POINTDIRBILLBOARD)
        mat4 modelMatrix = iModelMatrix;
        vec3 worldPos = GetWorldPos(modelMatrix);
        gl_Position = GetClipPos(worldPos);
        vTexCoord = iTexCoord;
        vSize = iTexCoord1;
        vWorldPos = vec4(worldPos, GetDepth(gl_Position));
        #ifdef POINTDIRBILLBOARD
            vNormal = iNormal;
        #endif
    #else
        mat4 modelMatrix = iModelMatrix;
        vec3 worldPos = GetWorldPos(modelMatrix);
        gl_Position = GetClipPos(worldPos);
        vTexCoord = GetTexCoord(iTexCoord);
        vWorldPos = vec4(worldPos, GetDepth(gl_Position));
    #endif

    #if defined(SOFTPARTICLES) && !defined(POINTEXPAND)
        vScreenPos = GetScreenPos(gl_Position);
    #endif

    #ifdef VERTEXCOLOR
        vColor = iColor;
    #endif

    #if !defined(POINTEXPAND)
        #ifdef PERPIXEL
            // Per-pixel forward lighting
            vec4 projWorldPos = vec4(worldPos, 1.0);

            #ifdef SHADOW
                // Shadow projection: transform from world space to shadow space
                for (int i = 0; i < NUMCASCADES; i++)
                    vShadowPos[i] = GetShadowPos(i, vec3(0, 0, 0), projWorldPos);
            #endif

            #ifdef SPOTLIGHT
                // Spotlight projection: transform from world space to projector texture coordinates
                vSpotPos = projWorldPos * cLightMatrices[0];
            #endif
        
            #ifdef POINTLIGHT
                vCubeMaskVec = (worldPos - cLightPos.xyz) * mat3(cLightMatrices[0][0].xyz, cLightMatrices[0][1].xyz, cLightMatrices[0][2].xyz);
            #endif
        #else
            // Ambient & per-vertex lighting
            vVertexLight = GetAmbient(GetZonePos(worldPos));

            #ifdef NUMVERTEXLIGHTS
                for (int i = 0; i < NUMVERTEXLIGHTS; ++i)
                    vVertexLight += GetVertexLightVolumetric(i, worldPos) * cVertexLights[i * 3].rgb;
            #endif
        #endif
    #endif
}

#endif

#ifdef COMPILEGS

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
        #ifdef PERPIXEL
            #ifdef SHADOW
                #ifndef GL_ES
                    vec4 vShadowPos[NUMCASCADES];
                #else
                    highp vec4 vShadowPos[NUMCASCADES];
                #endif
            #endif
            #ifdef SPOTLIGHT
                vec4 vSpotPos;
            #endif
            #ifdef POINTLIGHT
                vec3 vCubeMaskVec;
            #endif
        #else
            vec3 vVertexLight;
        #endif
    } gs_out;
    
    #define oTexCoord gs_out.vTexCoord
    #define oWorldPos gs_out.vWorldPos
    #define oColor gs_out.vColor
    #define oScreenPos gs_out.vScreenPos
    #define oShadowPos gs_out.vShadowPos
    #define oSpotPos gs_out.vSpotPos
    #define oCubeMaskVec gs_out.vCubeMaskVec
    #define oVertexLight gs_out.vVertexLight
    
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
    
layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

void GS()
{
    vec2 minUV = vs_out[0].vTexCoord.xy;
    vec2 maxUV = vs_out[0].vTexCoord.zw;
    vec2 partSize = vs_out[0].vSize.xy;
    
    vec2 partUV[] = {
        minUV,
        vec2(minUV.x, maxUV.y),
        vec2(maxUV.x, minUV.y),
        maxUV,
    };
    
    float s = sin(vs_out[0].vSize.z * M_DEGTORAD);
    float c = cos(vs_out[0].vSize.z * M_DEGTORAD);
    vec3 vUpNew    = c * vec3(1, 0, 0) - s * vec3(0, 1, 0);
    vec3 vRightNew = s * vec3(1, 0, 0) + c * vec3(0, 1, 0);
    vUpNew *= partSize.y;
    vRightNew *= partSize.x;
    
    vec3 partOffset[] = {
        -vUpNew + -vRightNew,
        -vUpNew + vRightNew,
        vUpNew + -vRightNew,
        vUpNew + vRightNew,
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
        
        #ifdef SOFTPARTICLES
            oScreenPos = GetScreenPos(gl_Position);
        #endif
        
        #ifdef PERPIXEL
            // Per-pixel forward lighting
            vec4 projWorldPos = vec4(oWorldPos.xyz, 1.0);

            #ifdef SHADOW
                // Shadow projection: transform from world space to shadow space
                for (int i = 0; i < NUMCASCADES; i++)
                    oShadowPos[i] = GetShadowPos(i, vec3(0, 0, 0), projWorldPos);
            #endif

            #ifdef SPOTLIGHT
                // Spotlight projection: transform from world space to projector texture coordinates
                oSpotPos = projWorldPos * cLightMatrices[0];
            #endif

            #ifdef POINTLIGHT
                oCubeMaskVec = (oWorldPos.xyz - cLightPos.xyz) * mat3(cLightMatrices[0][0].xyz, cLightMatrices[0][1].xyz, cLightMatrices[0][2].xyz);
            #endif
        #else
            // Ambient & per-vertex lighting
            oVertexLight = GetAmbient(GetZonePos(oWorldPos.xyz));

            #ifdef NUMVERTEXLIGHTS
                for (int i = 0; i < NUMVERTEXLIGHTS; ++i)
                    oVertexLight += GetVertexLightVolumetric(i, oWorldPos) * cVertexLights[i * 3].rgb;
            #endif
        #endif
        
        EmitVertex();
    }
    EndPrimitive();
}

#endif

#ifdef COMPILEPS

#ifdef GL3
    #define vTexCoord vs_out.vTexCoord
    #define vWorldPos vs_out.vWorldPos
    #define vColor vs_out.vColor
    #define vScreenPos vs_out.vScreenPos
    #define vVertexLight vs_out.vVertexLight
#endif

void PS()
{
    // Get material diffuse albedo
    #ifdef DIFFMAP
        vec4 diffInput = texture2D(sDiffMap, vTexCoord);
        #ifdef ALPHAMASK
            if (diffInput.a < 0.5)
                discard;
        #endif
        vec4 diffColor = cMatDiffColor * diffInput;
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

        diffColor.a = max(diffColor.a - fade, 0.0);
    #endif

    #ifdef PERPIXEL
        // Per-pixel forward lighting
        vec3 lightColor;
        vec3 lightDir;
        vec3 finalColor;

        float diff = GetDiffuseVolumetric(vWorldPos.xyz);

        #ifdef SHADOW
            diff *= GetShadow(vShadowPos, vWorldPos.w);
        #endif
    
        #if defined(SPOTLIGHT)
            lightColor = vSpotPos.w > 0.0 ? texture2DProj(sLightSpotMap, vSpotPos).rgb * cLightColor.rgb : vec3(0.0, 0.0, 0.0);
        #elif defined(CUBEMASK)
            lightColor = textureCube(sLightCubeMap, vCubeMaskVec).rgb * cLightColor.rgb;
        #else
            lightColor = cLightColor.rgb;
        #endif

        finalColor = diff * lightColor * diffColor.rgb;
        gl_FragColor = vec4(GetLitFog(finalColor, fogFactor), diffColor.a);
    #else
        // Ambient & per-vertex lighting
        vec3 finalColor = vVertexLight * diffColor.rgb;

        gl_FragColor = vec4(GetFog(finalColor, fogFactor), diffColor.a);
    #endif
}

#endif
