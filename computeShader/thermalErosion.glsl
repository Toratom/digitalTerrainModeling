layout(std430, binding = 0) readonly buffer ThickR { //avec std140 pb car stride de 4 donc avant on ecrivait dans le padding !?
	float ThicknessR[][NB_OF_LAYERS]; //Faire la même chose pour les water props, utiliser des arrays d'array
};

layout(std430, binding = 1) writeonly buffer ThickW {
	float ThicknessW[][NB_OF_LAYERS];
};

uniform uint gridHeight, gridWidth; 
//Mettre aussi la taille d'une cellule

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

void main() {
	const int kSize = 5; //Mettre un nb impair !!
	const int halfKSize = (kSize - 1) / 2;

	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);
	const int i = pointIJ.x;
	const int j = pointIJ.y;


	if (i < gridHeight && j < gridWidth) {
		float newH = 0;
		float currentH = 0;
		for (int dI = -halfKSize; dI < halfKSize + 1; dI += 1) {
			for (int dJ = -halfKSize; dJ < halfKSize + 1; dJ += 1) {
				newH += getHeight(i + dI, j + dJ) / (kSize * kSize);
			}
		}

		ThicknessW[getIndex(i, j)][0] = newH;
	}
}