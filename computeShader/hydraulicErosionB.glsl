//Realise l'étape 4
layout(std430, binding = 0) readonly buffer SedR { //avec std140 pb car stride de 4 donc avant on ecrivait dans le padding !?
	float SedimentR[]; //Faire la même chose pour les water props, utiliser des arrays d'array
};

layout(std430, binding = 1) writeonly buffer SedW {
	float SedimentW[];
};

layout(std430, binding = 2) readonly buffer Vel {
	float Velocity[][2]; //Vitesse en i, vitesse en j
};

uniform int gridHeight, gridWidth; 
uniform float dt;

layout(local_size_x = PATCH_HEIGHT, local_size_y = PATCH_WIDTH, local_size_z = 1) in;
//x correspond à i et y correspond à j

uint getIndex(int i, int j) {
	//i, j peuvent être en dehors de la grille, pas un probleme car clamp
	int iC = clamp(i, 0, gridHeight - 1);
	int jC = clamp(j, 0, gridWidth - 1);
	return gridWidth * iC + jC;
}

float getSediment(int i, int j) { //Si en dehors de la grille pas de sediment, car sinon risque de creer des sediments
	if (i < 0 || j < 0) return 0.;
	if (i > gridHeight - 1 || j > gridWidth - 1) return 0.;
	return SedimentR[getIndex(i, j)];
}

void main() {
	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);

	if (pointIJ.x < gridHeight && pointIJ.y < gridWidth) {
		float sourcePointAbsI = pointIJ.x - Velocity[getIndex(pointIJ.x, pointIJ.y)][0] * dt;
		float sourcePointAbsJ = pointIJ.y - Velocity[getIndex(pointIJ.x, pointIJ.y)][1] * dt;
		int sourcePointI = int(floor(sourcePointAbsI)); //Floor la position avant pour les nbs negatif
		int sourcePointJ = int(floor(sourcePointAbsJ));
		float coeff1 = sourcePointAbsI - sourcePointI;
		float coeff2 = sourcePointAbsJ - sourcePointJ;

		float outputS = (1. - coeff1) * (1. - coeff2) * getSediment(sourcePointI, sourcePointJ) 
			+ coeff2 * (1 - coeff1) * getSediment(sourcePointI, sourcePointJ + 1)
			+ coeff1 * (1 - coeff2) * getSediment(sourcePointI + 1, sourcePointJ)
			+ coeff1 * coeff2 * getSediment(sourcePointI + 1, sourcePointJ + 1);

		SedimentW[getIndex(pointIJ.x, pointIJ.y)] = outputS;
	}
}