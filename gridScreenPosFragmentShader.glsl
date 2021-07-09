#version 330 core	     // Minimal GL version support expected from the GPU

in vec2 fColor;

out vec4 color;

void main() {
	color = vec4(fColor, 0.0, 1.0);
}
