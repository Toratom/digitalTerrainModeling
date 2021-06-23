//Permet de calculer les normales du terrain et la hauteur du terrain en chaque i,j avant le rendering
#define PI 3.14159265

layout(std430, binding = 0) readonly buffer ThickR {
	float ThicknessR[][NB_OF_LAYERS];
};

layout(std140, binding = 1) writeonly buffer NormW {
	vec4 NormalsW[];
};

layout(std430, binding = 2) writeonly buffer HW {
	float HeightsW[];
};

layout(std140, binding = 3) writeonly buffer ColW {
	vec4 ColorsW[];
};

layout(std430, binding = 4) readonly buffer Vel {
	float Velocity[][2];
};

layout(std430, binding = 5) readonly buffer Sed {
	float Sediment[];
};

uniform int gridHeight, gridWidth; 
uniform float cellHeight, cellWidth;
uniform vec4 layersColor[NB_OF_LAYERS];
uniform bool renderWater, displayVel, displaySed;

layout(local_size_x = PATCH_HEIGHT, local_size_y = PATCH_WIDTH, local_size_z = 1) in;
//x correspond à i et y correspond à j



uint getIndex(int i, int j) {
	//i, j peuvent être en dehors de la grille, pas un probleme car clamp
	int iC = clamp(i, 0, gridHeight - 1);
	int jC = clamp(j, 0, gridWidth - 1);
	return gridWidth * iC + jC;
}

float getHeight(int i, int j, bool renderWater) {
	uint kMax = NB_OF_LAYERS - 1;
	if (renderWater) {
		kMax = NB_OF_LAYERS;
	}
	
	float h = 0;
	for (uint k = 0; k < kMax; k += 1) {
		h += ThicknessR[getIndex(i , j)][k];
	}
	return h;
}

float getIDerivate(int i, int j, bool renderWater) {
	//On fait attention au bord
	if (i == 0) {
		return (getHeight(i + 1, j, renderWater) - getHeight(i, j, renderWater)) / cellHeight;
	}
	if (i == gridHeight - 1) {
		return (getHeight(i, j, renderWater) - getHeight(i - 1, j, renderWater)) / cellHeight;
	}

	return (getHeight(i + 1, j, renderWater) - getHeight(i - 1, j, renderWater)) / (2.f * cellHeight);
}

float getJDerivate(int i, int j, bool renderWater) {
	//On fait attention au bord
	if (j == 0) {
		return (getHeight(i, j + 1, renderWater) - getHeight(i, j, renderWater)) / cellWidth;
	}
	if (j == gridWidth - 1) {
		return (getHeight(i, j, renderWater) - getHeight(i, j - 1, renderWater)) / cellWidth;
	}

	//return (getTerrainH(i, j) - getTerrainH(i, j - 1)) / m_cellWidth;
	return (getHeight(i, j + 1, renderWater) - getHeight(i, j - 1, renderWater)) / (2.f * cellWidth);
}

vec2 getGradient(int i, int j, bool renderWater) {
	return vec2(getIDerivate(i, j, renderWater), getJDerivate(i, j, renderWater));
}

uint getTopLayerId(int i, int j, bool renderWater) {
	uint id = NB_OF_LAYERS - 2;
	if (renderWater) {
		id = NB_OF_LAYERS - 1;
	}

	while (ThicknessR[getIndex(i, j)][id] < 0.0001 && id > 0) { //Si le layer pas assez epais on ne le concidere pas comme afleurant resout pb de clignotage
		id = id - 1;
	}

	return id;
}

vec3 hueToRGB(float H) { //Avec H dans 0, 360
	//Conversion de HSV en RGB
	float V = 0.5;
	float S = 1;
	float C = V * S;
	H = H / 60.;
	float X = C * (1 - abs(mod(H, 2) - 1));

	vec3 RGB;
	if ((0. <= H) && (H <= 1.)) {
		RGB.x = C;
		RGB.y = X;
		RGB.z = 0;
	}
	else if ((1. < H) && (H <= 2.)) {
		RGB.x = X;
		RGB.y = C;
		RGB.z = 0;
	}
	else if ((2. < H) && (H <= 3.)) {
		RGB.x = 0;
		RGB.y = C;
		RGB.z = X;
	}
	else if ((3. < H) && (H <= 4.)) {
		RGB.x = 0;
		RGB.y = X;
		RGB.z = C;
	}
	else if ((4. < H) && (H <= 5.)) {
		RGB.x = X;
		RGB.y = 0;
		RGB.z = C;
	}
	else if ((5. < H) && (H <= 6.)) {
		RGB.x = C;
		RGB.y = 0;
		RGB.z = X;
	}

	float m = V - C;
	RGB += vec3(m, m, m);

	return RGB;
}

vec4 velocityToColor(int i, int j) {
	vec2 vel = vec2(Velocity[getIndex(i, j)][0], Velocity[getIndex(i, j)][1]);
	vel = normalize(vel);
	float theta = PI / 2. + atan(vel.x / vel.y); //Dans 0, 2PI
	if (vel.y < 0) {
		theta += PI;
	}
	theta = (theta / (2. * PI)) * 360.; //Dans 0, 360

	vec3 RGB = hueToRGB(theta);

	return vec4(RGB, 1.);
}

vec4 sedimentToColor(int i, int j) {
	float sedC = clamp(100. * Sediment[getIndex(i, j)], 0., 1.); //Paramètre 100 à tune à la main ajouter slider sur IMGUI ?
	vec4 pointColor = layersColor[getTopLayerId(i, j, renderWater)];
	return vec4((1. - sedC) * pointColor.xyz + sedC * vec3(1., 0., 0.), pointColor.w);
}


void main() {
	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);
	const int i = pointIJ.x;
	const int j = pointIJ.y;
	vec2 grad;
	float h = 0.;

	if (i < gridHeight && j < gridWidth) {
		//Update normales
		grad = getGradient(i, j, renderWater);
		NormalsW[getIndex(i ,j)] = vec4(normalize(vec3(-grad.y, 1, -grad.x)), 0);

		//Upadte heights
		h = getHeight(i, j, renderWater);
		if ((renderWater) && (getTopLayerId(i, j, renderWater) < NB_OF_LAYERS - 1)) { //Si c'est pas de l'eau et fait rendu de l'eau enleve un eps
			h -= 0.001;
		}
		HeightsW[getIndex(i, j)] = h;

		//Update color
		ColorsW[getIndex(i, j)] = layersColor[getTopLayerId(i, j, renderWater)];
		if (displaySed && renderWater) ColorsW[getIndex(i, j)] = sedimentToColor(i, j);
		if (displayVel && renderWater) ColorsW[getIndex(i, j)] = velocityToColor(i, j);
	}
}