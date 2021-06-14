#define PATCH_W_NEIGHBORHOOD_HEIGHT (PATCH_HEIGHT + 2)
#define PATCH_W_NEIGHBORHOOD_WIDTH (PATCH_WIDTH + 2)

layout(std430, binding = 0) readonly buffer ThickR { //avec std140 pb car stride de 4 donc avant on ecrivait dans le padding !?
	float ThicknessR[][NB_OF_LAYERS]; //Faire la même chose pour les water props, utiliser des arrays d'array
};

layout(std430, binding = 1) writeonly buffer ThickW {
	float ThicknessW[][NB_OF_LAYERS];
};

uniform int gridHeight, gridWidth; 
uniform float cellHeight, cellWidth, dt;
uniform float erosionCoeffs[NB_OF_LAYERS];
uniform float thetasLimit[NB_OF_LAYERS];

const ivec2 offset = ivec2(1, 1);
const int indexOfSand = NB_OF_LAYERS - 2;
const ivec2 neighborTranslations[8] = ivec2[8](
	ivec2(-1, -1),
	ivec2(-1, 0),
	ivec2(-1, 1),
	ivec2(0, -1),
	ivec2(0, 1),
	ivec2(1, -1),
	ivec2(1, 0),
	ivec2(1, 1));

shared float dHOut[PATCH_W_NEIGHBORHOOD_HEIGHT][PATCH_W_NEIGHBORHOOD_WIDTH][8]; //8 connexite


layout(local_size_x = PATCH_HEIGHT, local_size_y = PATCH_WIDTH, local_size_z = 1) in;
//x correspond à i et y correspond à j


uint getIndex(int i, int j) {
	//i, j peuvent être en dehors de la grille, pas un probleme car clamp
	int iC = clamp(i, 0, gridHeight - 1);
	int jC = clamp(j, 0, gridWidth - 1);
	return gridWidth * iC + jC;
}

float getHeight(int i, int j) {
	//Thermal Erosion du hauteur du terrain sans eau d'ou NB_OF_LAYERS - 1
	float h = 0;
	for (uint k = 0; k < NB_OF_LAYERS - 1; k += 1) {
		h += ThicknessR[getIndex(i , j)][k];
	}
	return h;
}

uint getTopLayerId(int i, int j) {
	//Thermal Erosion du hauteur du terrain sans eau d'ou NB_OF_LAYERS - 2
	uint id = NB_OF_LAYERS - 2;
	while (ThicknessR[getIndex(i, j)][id] == 0. && id > 0) {
		id = id - 1;
	}

	return id;
}

void main() {
	//i, j indice global VS x, y indice locale
	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);
	const ivec2 pointXY = ivec2(gl_LocalInvocationID);

	ivec2 currentIJ;
	ivec2 neighborIJ;
	ivec2 currentXY;
	ivec2 neighborXY;
	float distanceToNextCell;
	uint currentTopLayerId;
	float dHOutTot = 0.;
	float K = 1.; 
	

	//Phase 1 : Calculs des dHOut pour PATCH_W_NEIGHBORHOOD (la double boucle permet d'explorer tous le patch w neigbour intelligemment)
	//On fait cette phase en dehors du if pas un probleme car on clamp les coord
	//Besoin de le sortir du if pour que la retire phase soit en dehors du if, car sinon il arrive que certain thread n'atteignent pas la barriere (car dans un if) ce qui bloque tout !!
	for (int dI = 0; dI < PATCH_W_NEIGHBORHOOD_HEIGHT; dI += PATCH_HEIGHT) {
		for (int dJ = 0; dJ < PATCH_W_NEIGHBORHOOD_WIDTH; dJ += PATCH_WIDTH) {
			if (pointXY.x + dI < PATCH_W_NEIGHBORHOOD_HEIGHT && pointXY.y + dJ < PATCH_W_NEIGHBORHOOD_WIDTH) {
				currentIJ = pointIJ + ivec2(dI, dJ) - offset;
				currentXY = pointXY + ivec2(dI, dJ);
				currentTopLayerId = getTopLayerId(currentIJ.x, currentIJ.y);

				//Cas bedrock gere par erosionCoeff = 0 grace à la liste erosionCoeffs
				for (uint k = 0; k < 8; k += 1) {
					neighborIJ = currentIJ + neighborTranslations[k];
					distanceToNextCell = sqrt(pow(cellHeight * neighborTranslations[k].x, 2) + pow(cellWidth * neighborTranslations[k].y, 2));
					dHOut[currentXY.x][currentXY.y][k] = max(0, erosionCoeffs[currentTopLayerId] * ((getHeight(currentIJ.x, currentIJ.y) - getHeight(neighborIJ.x, neighborIJ.y)) / distanceToNextCell - tan(thetasLimit[currentTopLayerId])) * dt);
					dHOutTot += dHOut[currentXY.x][currentXY.y][k];
				}

				//Calcule du coeff de normalisation K, pour eviter de perdre plus de matiere que la colonne courante du layer affleurant
				K = min(1, ThicknessR[getIndex(currentIJ.x, currentIJ.y)][currentTopLayerId] / dHOutTot);

				//Applique normalisation
				for (uint k = 0; k < 8; k += 1) {
					dHOut[currentXY.x][currentXY.y][k] = K * dHOut[currentXY.x][currentXY.y][k];
				}
			}
		}
	}

	//Retire Phase
	barrier();
	memoryBarrierShared();

	if (pointIJ.x < gridHeight && pointIJ.y < gridWidth) {
	//On update la hauteur que si le point est dans la grille d'ou le if
		//Phase 2 : calcul de ce qui entre, ce qui sort
		currentXY = pointXY + offset; //Localisation de pointXY dans dHOut
		dHOutTot = 0.;
		float dHInTot = 0.;
		uint kN = 0;
		for (uint k = 0; k < 8; k += 1) {
			dHOutTot += dHOut[currentXY.x][currentXY.y][k];
			
			neighborIJ = currentXY + neighborTranslations[k];
			kN = 7 - k;
			dHInTot += dHOut[neighborIJ.x][neighborIJ.y][kN];
		}

		//Phase 3 : ecriture avec convention que l'érodé devient sable, besoin de separer en R et W a cause des frontieres/bords des patchs qui se trouve dans plusieurs working group
		float thicknessOutput[NB_OF_LAYERS]; //On passe par cette variable car on n'a mis ThicknessW en W only donc on ne peut pas utiliser le += et -= + permet de regler le pb ou currentTopLayerId = sandId
		for (uint k = 0; k < NB_OF_LAYERS; k += 1) {
			thicknessOutput[k] = ThicknessR[getIndex(pointIJ.x, pointIJ.y)][k];
		}
		currentTopLayerId = getTopLayerId(pointIJ.x, pointIJ.y);
		thicknessOutput[currentTopLayerId] -= dHOutTot;
		thicknessOutput[indexOfSand] += dHInTot;
		for (uint k = 0; k < NB_OF_LAYERS; k += 1) {
			ThicknessW[getIndex(pointIJ.x, pointIJ.y)][k] = thicknessOutput[k];
		}
	}
}