#version 330 core            // Minimal GL version support expected from the GPU

layout(location = 0) in vec2 vGridCoords; //I, J donc Z, X
layout(location = 1) in float vHeight;
//layout(location = 2) in vec4 vNormal; //Pas besoin ici, mais laisse pour n'avoir qu'un seul vao ?
//layout(location = 3) in vec4 vColor; //Pas besoin ici
layout(location = 4) in uvec2 vGridPos; //Le (i,j), pas de probleme dans conversion en float lors de la division ?

out vec2 fColor;

uniform mat4 viewMat, projMat, modelMat;
uniform int gridHeight, gridWidth;

void main() {
    float height = vHeight;
    vec3 vPosition = vec3(vGridCoords.y, height, vGridCoords.x);
    gl_Position = projMat * viewMat * modelMat * vec4(vPosition, 1.0); // mandatory to rasterize properly

    fColor.x = vGridPos.x / float(gridHeight); //X correspond à I et Y a J
    fColor.y = vGridPos.y / float(gridWidth);
}
