//Permet de calculer les normales du terrain et la hauteur du terrain en chaque i,j avant le rendering

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

uniform int gridHeight, gridWidth; 
uniform float cellHeight, cellWidth;
uniform vec4 layersColor[NB_OF_LAYERS];
uniform bool renderWater;

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
	}
}