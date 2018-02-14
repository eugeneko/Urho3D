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

#include "../Core/Context.h"
#include "../Graphics/Graphics.h"
#include "../Graphics/Shader.h"
#include "../Graphics/ShaderVariation.h"
#include "../IO/Deserializer.h"
#include "../IO/FileSystem.h"
#include "../IO/Log.h"
#include "../Resource/ResourceCache.h"

#include "../DebugNew.h"

namespace Urho3D
{

void CommentOutFunction(String& code, const String& signature)
{
    unsigned startPos = code.Find(signature);
    unsigned braceLevel = 0;
    if (startPos == String::NPOS)
        return;

    code.Insert(startPos, "/*");

    for (unsigned i = startPos + 2 + signature.Length(); i < code.Length(); ++i)
    {
        if (code[i] == '{')
            ++braceLevel;
        else if (code[i] == '}')
        {
            --braceLevel;
            if (braceLevel == 0)
            {
                code.Insert(i + 1, "*/");
                return;
            }
        }
    }
}

Shader::Shader(Context* context) :
    Resource(context),
    timeStamp_(0),
    numVariations_(0)
{
    RefreshMemoryUse();
}

Shader::~Shader()
{
    auto* cache = GetSubsystem<ResourceCache>();
    if (cache)
        cache->ResetDependencies(this);
}

void Shader::RegisterObject(Context* context)
{
    context->RegisterFactory<Shader>();
}

bool Shader::BeginLoad(Deserializer& source)
{
    auto* graphics = GetSubsystem<Graphics>();
    if (!graphics)
        return false;

    // Load the shader source code and resolve any includes
    timeStamp_ = 0;
    String shaderCode;
    if (!ProcessSource(shaderCode, source))
        return false;

    // Comment out the unneeded shader function
    vertexShader_.sourceCode_ = shaderCode;
    pixelShader_.sourceCode_ = shaderCode;

    CommentOutFunction(vertexShader_.sourceCode_, "void PS(");
    CommentOutFunction(vertexShader_.sourceCode_, "void GS(");
    CommentOutFunction(vertexShader_.sourceCode_, "void HS(");
    CommentOutFunction(vertexShader_.sourceCode_, "void DS(");
    
    CommentOutFunction(pixelShader_.sourceCode_, "void VS(");
    CommentOutFunction(pixelShader_.sourceCode_, "void GS(");
    CommentOutFunction(pixelShader_.sourceCode_, "void HS(");
    CommentOutFunction(pixelShader_.sourceCode_, "void DS(");

#if !defined(GL_ES_VERSION_2_0) && !defined(URHO3D_D3D9)
    geometryShader_.sourceCode_ = shaderCode;
    tcsShader_.sourceCode_ = shaderCode;
    tesShader_.sourceCode_ = shaderCode;

    CommentOutFunction(geometryShader_.sourceCode_, "void VS(");
    CommentOutFunction(geometryShader_.sourceCode_, "void PS(");
    CommentOutFunction(geometryShader_.sourceCode_, "void HS(");
    CommentOutFunction(geometryShader_.sourceCode_, "void DS(");
    
    CommentOutFunction(tcsShader_.sourceCode_, "void VS(");
    CommentOutFunction(tcsShader_.sourceCode_, "void PS(");
    CommentOutFunction(tcsShader_.sourceCode_, "void GS(");
    CommentOutFunction(tcsShader_.sourceCode_, "void DS(");

    CommentOutFunction(tesShader_.sourceCode_, "void VS(");
    CommentOutFunction(tesShader_.sourceCode_, "void PS(");
    CommentOutFunction(tesShader_.sourceCode_, "void GS(");
    CommentOutFunction(tesShader_.sourceCode_, "void HS(");
#endif


    // OpenGL: rename either VS() or PS() or GS() to main()
#ifdef URHO3D_OPENGL
    vertexShader_.sourceCode_.Replace("void VS(", "void main(");
    pixelShader_.sourceCode_.Replace("void PS(", "void main(");
#ifndef GL_ES_VERSION_2_0
    geometryShader_.sourceCode_.Replace("void GS(", "void main(");
    tcsShader_.sourceCode_.Replace("void HS(", "void main(");
    tesShader_.sourceCode_.Replace("void DS(", "void main(");
#endif
#endif

    RefreshMemoryUse();
    return true;
}

bool Shader::EndLoad()
{
    // If variations had already been created, release them and require recompile
    for (HashMap<StringHash, SharedPtr<ShaderVariation> >::Iterator i = vertexShader_.variations_.Begin(); i != vertexShader_.variations_.End(); ++i)
        i->second_->Release();
    for (HashMap<StringHash, SharedPtr<ShaderVariation> >::Iterator i = pixelShader_.variations_.Begin(); i != pixelShader_.variations_.End(); ++i)
        i->second_->Release();
    
#if !defined(GL_ES_VERSION_2_0) && !defined(URHO3D_D3D9)
    for (HashMap<StringHash, SharedPtr<ShaderVariation> >::Iterator i = geometryShader_.variations_.Begin(); i != geometryShader_.variations_.End(); ++i)
        i->second_->Release();
    for (HashMap<StringHash, SharedPtr<ShaderVariation> >::Iterator i = tcsShader_.variations_.Begin(); i != tcsShader_.variations_.End(); ++i)
        i->second_->Release();
    for (HashMap<StringHash, SharedPtr<ShaderVariation> >::Iterator i = tesShader_.variations_.Begin(); i != tesShader_.variations_.End(); ++i)
        i->second_->Release();
#endif

    return true;
}

ShaderVariation* Shader::GetVariation(ShaderType type, const String& defines)
{
    return GetVariation(type, defines.CString());
}

const String& Shader::GetSourceCode(ShaderType type) const
{
    switch (type)
    {
    case VS:
        return vertexShader_.sourceCode_;
    case PS:
        return pixelShader_.sourceCode_;
#if !defined(GL_ES_VERSION_2_0) && !defined(URHO3D_D3D9)
    case GS:
        return geometryShader_.sourceCode_;
    case TCS:
        return tcsShader_.sourceCode_;
    case TES:
        return tesShader_.sourceCode_;
#endif
    }
    return vertexShader_.sourceCode_;
}

ShaderVariation* Shader::GetVariation(ShaderType type, const char* defines)
{
    StringHash definesHash(defines);
    HashMap<StringHash, SharedPtr<ShaderVariation> >& variations = GetVariations(type);
    HashMap<StringHash, SharedPtr<ShaderVariation> >::Iterator i = variations.Find(definesHash);
    if (i == variations.End())
    {
        // If shader not found, normalize the defines (to prevent duplicates) and check again. In that case make an alias
        // so that further queries are faster
        String normalizedDefines = NormalizeDefines(defines);
        StringHash normalizedHash(normalizedDefines);

        i = variations.Find(normalizedHash);
        if (i != variations.End())
            variations.Insert(MakePair(definesHash, i->second_));
        else
        {
            // No shader variation found. Create new
            i = variations.Insert(MakePair(normalizedHash, SharedPtr<ShaderVariation>(new ShaderVariation(this, type))));
            if (definesHash != normalizedHash)
                variations.Insert(MakePair(definesHash, i->second_));

            i->second_->SetName(GetFileName(GetName()));
            i->second_->SetDefines(normalizedDefines);
            ++numVariations_;
            RefreshMemoryUse();
        }
    }

    return i->second_;
}

HashMap<StringHash, SharedPtr<ShaderVariation> >& Shader::GetVariations(ShaderType type)
{
    switch (type)
    {
    case VS:
        return vertexShader_.variations_;
    case PS:
        return pixelShader_.variations_;
#if !defined(GL_ES_VERSION_2_0) && !defined(URHO3D_D3D9)
    case GS:
        return geometryShader_.variations_;
    case TCS:
        return tcsShader_.variations_;
    case TES:
        return tesShader_.variations_;
#endif
    }
    return vertexShader_.variations_;
}

bool Shader::ProcessSource(String& code, Deserializer& source)
{
    auto* cache = GetSubsystem<ResourceCache>();

    // If the source if a non-packaged file, store the timestamp
    auto* file = dynamic_cast<File*>(&source);
    if (file && !file->IsPackaged())
    {
        auto* fileSystem = GetSubsystem<FileSystem>();
        String fullName = cache->GetResourceFileName(file->GetName());
        unsigned fileTimeStamp = fileSystem->GetLastModifiedTime(fullName);
        if (fileTimeStamp > timeStamp_)
            timeStamp_ = fileTimeStamp;
    }

    // Store resource dependencies for includes so that we know to reload if any of them changes
    if (source.GetName() != GetName())
        cache->StoreResourceDependency(this, source.GetName());

    while (!source.IsEof())
    {
        String line = source.ReadLine();

        if (line.StartsWith("#include"))
        {
            String includeFileName = GetPath(source.GetName()) + line.Substring(9).Replaced("\"", "").Trimmed();

            SharedPtr<File> includeFile = cache->GetFile(includeFileName);
            if (!includeFile)
                return false;

            // Add the include file into the current code recursively
            if (!ProcessSource(code, *includeFile))
                return false;
        }
        else
        {
            code += line;
            code += "\n";
        }
    }

    // Finally insert an empty line to mark the space between files
    code += "\n";

    return true;
}

String Shader::NormalizeDefines(const String& defines)
{
    Vector<String> definesVec = defines.ToUpper().Split(' ');
    Sort(definesVec.Begin(), definesVec.End());
    return String::Joined(definesVec, " ");
}

void Shader::RefreshMemoryUse()
{
    unsigned sourcesSize = vertexShader_.sourceCode_.Length() + pixelShader_.sourceCode_.Length();
#if !defined(GL_ES_VERSION_2_0) && !defined(URHO3D_D3D9)
    sourcesSize += geometryShader_.sourceCode_.Length();
    sourcesSize += tcsShader_.sourceCode_.Length() + tesShader_.sourceCode_.Length();
#endif
    SetMemoryUse(
        (unsigned)(sizeof(Shader) + sourcesSize + numVariations_ * sizeof(ShaderVariation)));
}

}
