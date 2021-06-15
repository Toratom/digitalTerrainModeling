#version 330 core            // Minimal GL version support expected from the GPU

layout(location=0) in vec3 vPosition; // The 1st input attribute is the position (CPU side: glVertexAttrib 0)
layout(location=1) in vec3 vNormal;
layout(location=2) in vec4 vColor;

out vec3 fNormal;
out vec3 fPos;
out vec4 fColor;

uniform mat4 viewMat, projMat, modelMat;



void main() {
    gl_Position = projMat * viewMat * modelMat* vec4(vPosition, 1.0); // mandatory to rasterize properly
    fNormal = transpose(inverse(mat3(modelMat))) * vNormal;
    fPos = vec3(modelMat * vec4(vPosition,1.0));
    fColor = vColor;
    //fColor = vec4(fNormal, 1.0);
    //fTexCoord = vTexCoord;
}
