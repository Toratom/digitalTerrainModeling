#version 330 core	     // Minimal GL version support expected from the GPU

struct Material {
	sampler2D albedoTex;
};
uniform Material material;
uniform vec3 camPos;

in vec2 fTexCoord;
in vec3 fPos;
in vec3 fNormal;
in vec3 fColor;

out vec4 color;	  // Shader output: the color response attached to this fragment

vec3 lightSun = vec3(1);
vec3 ka = vec3(0.8);
vec3 kd = vec3(0.5);
vec3 ks = vec3(0.8);
float alpha = 20;

void main() {
	
	/*vec3 normal = normalize(fNormal);
	vec3 lightDirection = normalize(lightSun-fPos);
	vec3 view = normalize(camPos-fPos);
	//vec3 reflection = 2*dot(normal,lightDirection)*normal - lightDirection;
	vec3 reflection = normalize(-reflect(lightDirection, lightDirection))
	vec3 ambient = ka;

	float ln = max(dot(lightDirection,normal),0.0);
	vec3 diffuse = vec3(ln)*kd;

	float rv = max(dot(reflection,view),0.0);
	vec3 specular = vec3(pow(rv,alpha))*ks;

	vec3 texColor = texture(material.albedoTex, fTexCoord).rgb;*/

	//color = vec4(texColor * (diffuse + ambient + specular), 1.0);// build an RGBA from an RGB
	//color = vec4(normalize(fNormal),1.0);
	color = vec4(fColor, 1.0);

	//color = vec4(normalize(fNormal).x, 0.f, normalize(fNormal).z, 1.0);
}
