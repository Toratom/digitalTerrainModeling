#version 330 core	     // Minimal GL version support expected from the GPU

struct Material {
	sampler2D albedoTex;
};
uniform Material material;
uniform vec3 camPos;

in vec2 fTexCoord;
in vec3 fPos;
in vec3 fNormal;
in vec4 fColor;

out vec4 color;	  // Shader output: the color response attached to this fragment

vec3 lightSun = vec3(-5.0, 5.0, -5.0);
vec4 ka = vec4(1.0);
vec4 kd = vec4(1.0);
vec4 ks = vec4(0.8, 0.8, 0.8, 1.0);
float alpha = 20;

void main() {
	
	vec3 normal = normalize(fNormal);
	vec3 lightDirection = normalize(lightSun - fPos);
	//vec3 view = normalize(camPos-fPos);
	//vec3 reflection = 2*dot(normal,lightDirection)*normal - lightDirection;
	//vec3 reflection = normalize(-reflect(lightDirection, lightDirection))
	vec4 ambient = ka;

	float ln = max(dot(lightDirection, normal), 0.0);
	vec4 diffuse = vec4(ln) * kd;

	//float rv = max(dot(reflection,view),0.0);
	//vec3 specular = vec3(pow(rv,alpha))*ks;
	vec4 specular = vec4(0.0);

	//vec3 texColor = texture(material.albedoTex, fTexCoord).rgb;

	color = vec4(fColor * (diffuse + ambient + specular));// build an RGBA from an RGB
	//color = fColor;
	//color = vec4(normalize(fNormal).x, 0.f, normalize(fNormal).z, 1.0);
	//color = vec4(normalize(fNormal),1.0);
}
