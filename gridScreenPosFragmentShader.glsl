#version 330 core	     // Minimal GL version support expected from the GPU

in vec2 fColor;

out vec2 color;

void main() {
	color = fColor;
	//color = vec2(1,1);
}
