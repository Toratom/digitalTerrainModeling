//L'entête est ajouté par le main lors avant la compilation du shader, l'entête est la suivante :
//#version 430
//#extension GL_ARB_compute_shader : enable
//#extension GL_ARB_shader_storage_buffer_object : enable
//#define NB_OF_LAYERS 1 //Faire la même chose pour les water props
//#define PATCH_HEIGHT 32
//#define PATCH_WIDTH 32

//Ce sont des interface block i.e. une interface qui permet de regrouper plusieurs variable GPU, par exemple on regroupe ici dans le tableau (buffer) ThicknessR des floats
//Il faut preciser comment sont rangees les donnes dans le block, le stride qui les separe (i.e. comment sont ranges les floats dans le buffer ou les uns à la suite des autres
//un float puis du padding pour avoir la taille d'un vec4 puis un autre float). Pour preciser la convention utiliser pour ranger les donnes dans le block on utilise layout(...)
//Le deuxieme exemple correspond à std140... ?
//On appelle membre les elements à l'interieure du block (dans notre cas des floats).
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
		for (int dI = -halfKSize; dI < halfKSize + 1; dI += 1) {
			for (int dJ = -halfKSize; dJ < halfKSize + 1; dJ += 1) {
				newH += getHeight(i + dI, j + dJ) / (kSize * kSize);
			}
		}

		for (uint k = 0; k < NB_OF_LAYERS; k += 1) {
			ThicknessW[getIndex(i, j)][k] = newH / NB_OF_LAYERS;
		}
	}
}