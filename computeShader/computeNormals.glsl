layout(std430, binding = 0) readonly buffer ThickR {
	float ThicknessR[][NB_OF_LAYERS]; //Faire la même chose pour les water props, utiliser des arrays d'array
};

layout(std140, binding = 1) writeonly buffer NormW {
	vec4 NormalsW[];
};

uniform uint gridHeight, gridWidth; 
uniform float cellHeight, cellWidth;

layout(local_size_x = PATCH_HEIGHT, local_size_y = PATCH_WIDTH, local_size_z = 1) in;
//x correspond à i et y correspond à j



uint getIndex(int i, int j) {
	//i, j peuvent être en dehors de la grille, pas un probleme car clamp
	uint iC = uint(clamp(i, 0, gridHeight));
	uint jC = uint(clamp(j, 0, gridWidth));
	//ivec2 xyC = clamp(ivec2(i, j), vec2(0, 0), vec2(gridHeight, gridWidth)); //Permet de ne pas se poser de question quand 
	return gridWidth * iC + jC;
}

float getHeight(int i, int j) {
	float h = 0;
	for (uint k = 0; k < NB_OF_LAYERS; k += 1) {
		h += ThicknessR[getIndex(i , j)][k];
	}
	return h;
}

float getIDerivate(int i, int j) {
	//On fait attention au bord
	if (i == 0) {
		return (getHeight(i + 1, j) - getHeight(i, j)) / cellHeight;
	}
	if (i == gridHeight - 1) {
		return (getHeight(i, j) - getHeight(i - 1, j)) / cellHeight;
	}

	return (getHeight(i + 1, j) - getHeight(i - 1, j)) / (2.f * cellHeight);
}

float getJDerivate(int i, int j) {
	//On fait attention au bord
	if (j == 0) {
		return (getHeight(i, j + 1) - getHeight(i, j)) / cellWidth;
	}
	if (j == gridWidth - 1) {
		return (getHeight(i, j) - getHeight(i, j - 1)) / cellWidth;
	}

	//return (getTerrainH(i, j) - getTerrainH(i, j - 1)) / m_cellWidth;
	return (getHeight(i, j + 1) - getHeight(i, j - 1)) / (2.f * cellWidth);
}

vec2 getGradient(int i, int j) {
	return vec2(getIDerivate(i, j), getJDerivate(i, j));
}



void main() {
	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);
	const int i = pointIJ.x;
	const int j = pointIJ.y;
	vec2 grad;

	if (i < gridHeight && j < gridWidth) {
		grad = getGradient(i, j);
		NormalsW[getIndex(i ,j)] = vec4(normalize(vec3(-grad.y, 1, -grad.x)), 0);
	}
}