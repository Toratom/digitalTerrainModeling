#version 330 core	     // Minimal GL version support expected from the GPU

struct Material {
	sampler2D albedoTex;
};
uniform Material material;
uniform vec3 camPos;
uniform vec2 cursorPos; //Position normaliser du cursor (i,j) dans [0,1] attention different des texCoords (u, v)
uniform bool displayEdition;
uniform sampler2D gridScreenPosTex;

//in vec2 fTexCoord;
in vec3 fCoords;
in vec3 fNormal;
in vec4 fColor;
in vec2 fGridPos;

out vec4 color;	  // Shader output: the color response attached to this fragment

vec3 lightSun = vec3(-5.0, 5.0, -5.0);
vec4 ka = vec4(1.);
vec4 kd = vec4(1.0);
vec4 ks = vec4(0.8, 0.8, 0.8, 1.);
float alpha = 20;

//TODO :
//Faire le test si cursorGridPos en dehors du terrain cestadire si negatif en les deux coords
//PB avec la texture et displayEdition...

void main() {
	vec2 cursorGridPos = texture(gridScreenPosTex, vec2(cursorPos.y, 1.0 - cursorPos.x)).xy;

	vec3 normal = normalize(fNormal);
	vec3 lightDirection = normalize(lightSun - fCoords);
	vec4 ambient = ka;

	float ln = max(dot(lightDirection, normal), 0.0);
	vec4 diffuse = ln * kd;

	vec4 specular = ks * vec4(0.0);

	vec4 basedColor = fColor;
	if (displayEdition && cursorGridPos.x >= 0.0 && cursorGridPos.y >= 0.0 && (distance(cursorGridPos, fGridPos) < 0.005)) {
		basedColor = vec4(0.5, 0.0, 0.0, 1.0);
	}

	color = basedColor * (diffuse + ambient + specular);// build an RGBA from an RGB
	//color = vec4(fColor, 1.0);
}
