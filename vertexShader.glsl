#version 330 core            // Minimal GL version support expected from the GPU

layout(location = 0) in vec2 vGridPos; //I, J donc Z, X
layout(location = 1) in float vHeight;
layout(location = 2) in vec4 vNormal; //Attention ici normal en vec4 car buffer de vec3 mal accepté dans shader storage buffer object
layout(location = 3) in vec4 vColor;


out vec3 fPos;
out vec3 fNormal;
out vec3 fColor;

uniform mat4 viewMat, projMat, modelMat;
uniform uint nbOfLayers;


void main() {
    float height = vHeight;
    vec3 vPosition = vec3(vGridPos.y, height, vGridPos.x);

    gl_Position = projMat * viewMat * modelMat* vec4(vPosition, 1.0); // mandatory to rasterize properly
    fNormal = transpose(inverse(mat3(modelMat))) * vNormal.xyz;
    fPos = vec3(modelMat*vec4(vPosition, 1.0));

    //fColor = vec3(fNormal.x, fNormal.y, fNormal.z);
    fColor = vColor.xyz;

    //fTexCoord = vTexCoord;
}
