#version 430

//#extension GL_ARB_compute_shader : enable
//#extension GL_ARB_shader_storage_buffer_object : enable

//Shader dans le cas terrain à 1 layer

#define NB_OF_LAYERS 1 //Faire la même chose pour les water props
#define PATCH_HEIGHT 32
#define PATCH_WIDTH 32
//ATTENTION doit etre coherent avec les defines dans main.cpp

layout(std140, binding = 0) buffer readonly ThickR {
	float ThicknessR[][NB_OF_LAYERS]; //Faire la même chose pour les water props, utiliser des arrays d'array
};

layout(std140, binding = 1) buffer writeonly ThickW {
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

//Pour test flou gaussien
//float getGaussWeight(float h, float sig) {
//	return exp(- (h * h) / (2.0 * sig * sig));
//}

void main() {
	//Pour test flou gaussien
	//const float sig = 2.0;
	const int kSize = 5; //Mettre un nb impair !!
	const int halfKSize = (kSize - 1) / 2;

	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);
	const int i = pointIJ.x;
	const int j = pointIJ.y;

	if (i < gridHeight && j < gridWidth) {
		//float newH = 0;
		//float currentH = 0;
		//for (int dI = -halfKSize; dI < halfKSize + 1; dI += 1) {
		//	for (int dJ = -halfKSize; dJ < halfKSize + 1; dJ += 1) {
		//		newH += getHeight(i + dI, j + dJ) / (kSize * kSize);
		//	}
		//}

		//ThicknessW[getIndex(i, j)][0] = newH;

		ThicknessW[getIndex(i, j)][0] = getHeight(i, j)/1.5;
	}
}