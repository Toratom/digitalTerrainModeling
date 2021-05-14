#version 330 core            // Minimal GL version support expected from the GPU

layout(location=0) in vec2 vGridPos; //I, J donc Z, X
layout(location=1) in float vLayersHeight; //Toujours vec4 mais nbOfLayers change

out vec3 fPos;
//out vec3 fNormal;
out vec3 fColor;

uniform mat4 viewMat, projMat, modelMat;
uniform uint nbOfLayers;


void main() {
    //C'est le vertex shader qui calcule la hauteur resultante, la couleur à partir de vertext attributs et list layerColor... (la passer en uniforme?)
    float height = vLayersHeight; //Sommer si vecI
    vec3 vPosition = vec3(vGridPos.y, height, vGridPos.x);

    gl_Position = projMat * viewMat * modelMat* vec4(vPosition, 1.0); // mandatory to rasterize properly
    //fNormal = transpose(inverse(mat3(modelMat)))*vNormal;
    fPos = vec3(modelMat*vec4(vPosition, 1.0));

    fColor = vec3(1.0, 0.0, 0.0);

    //fTexCoord = vTexCoord;
}
