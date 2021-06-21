//Realise les etapes 2 et 3 de l'article Fast Hydraulic Erosion on GPU

#define PATCH_W_NEIGHBORHOOD_HEIGHT (PATCH_HEIGHT + 2)
#define PATCH_W_NEIGHBORHOOD_WIDTH (PATCH_WIDTH + 2)

layout(std430, binding = 0) readonly buffer ThickR { //avec std140 pb car stride de 4 donc avant on ecrivait dans le padding !?
	float ThicknessR[][NB_OF_LAYERS]; //Faire la même chose pour les water props, utiliser des arrays d'array
};

layout(std430, binding = 1) writeonly buffer ThickW {
	float ThicknessW[][NB_OF_LAYERS];
};

layout(std430, binding = 2) readonly buffer FlowR {
	float FlowOutR[][4]; //B (i - 1, j) , L (i, j - 1), R (i, j + 1), T (i + 1, j)
};

layout(std430, binding = 3) writeonly buffer FlowW {
	float FlowOutW[][4]; //B (i - 1, j) , L (i, j - 1), R (i, j + 1), T (i + 1, j)
};

layout(std430, binding = 4) coherent buffer Vel {
	float Velocity[][2]; //Vitesse en i, vitesse en j
};

layout(std430, binding = 5) coherent buffer Sed {
	float Sediment[];
};

uniform int gridHeight, gridWidth; 
uniform float cellHeight, cellWidth, dt;

const ivec2 offset = ivec2(1, 1);
const int indexOfWater = NB_OF_LAYERS - 1;
const int indexOfSand = NB_OF_LAYERS - 2;
const ivec2 neighborTranslations[4] = ivec2[4](
	ivec2(-1, 0),
	ivec2(0, -1),
	ivec2(0, 1),
	ivec2(1, 0));
const float g = 9.81;
const float Kc = 0.5;
const float Ks = 0.5;
const float Kd = 0.5;

shared float flowOut[PATCH_W_NEIGHBORHOOD_HEIGHT][PATCH_W_NEIGHBORHOOD_WIDTH][4]; //4 connexite


layout(local_size_x = PATCH_HEIGHT, local_size_y = PATCH_WIDTH, local_size_z = 1) in;
//x correspond à i et y correspond à j

uint getIndex(int i, int j) {
	//i, j peuvent être en dehors de la grille, pas un probleme car clamp
	int iC = clamp(i, 0, gridHeight - 1);
	int jC = clamp(j, 0, gridWidth - 1);
	return gridWidth * iC + jC;
}

float getHeight(int i, int j) {
	//Hydraulic erosion dc hauteur du terrain avec eau d'ou NB_OF_LAYERS
	float h = 0;
	for (uint k = 0; k < NB_OF_LAYERS; k += 1) {
		h += ThicknessR[getIndex(i , j)][k];
	}
	return h;
}

float getFlow(int i, int j, unsigned int k) {
	if (i < 0 || j < 0) return 0.;
	if (i > gridHeight - 1 || j > gridWidth - 1) return 0.;
	return FlowOutR[getIndex(i, j)][k];
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
	float A = cellHeight * cellWidth; //La surface d'un pipe
	float deltaH = 0.;
	float distanceToNextCell = 0.;
	float flowOutTot = 0.;
	float K = 1.; 
	

	//Phase 1 : Calculs des outputFlow pour PATCH_W_NEIGHBORHOOD (la double boucle permet d'explorer tous le patch w neigbour intelligemment)
	for (int dI = 0; dI < PATCH_W_NEIGHBORHOOD_HEIGHT; dI += PATCH_HEIGHT) {
		for (int dJ = 0; dJ < PATCH_W_NEIGHBORHOOD_WIDTH; dJ += PATCH_WIDTH) {
			if (pointXY.x + dI < PATCH_W_NEIGHBORHOOD_HEIGHT && pointXY.y + dJ < PATCH_W_NEIGHBORHOOD_WIDTH) {
				currentIJ = pointIJ + ivec2(dI, dJ) - offset;
				currentXY = pointXY + ivec2(dI, dJ);
				flowOutTot = 0.;

				for (uint k = 0; k < 4; k += 1) {
					neighborIJ = currentIJ + neighborTranslations[k];
					distanceToNextCell = sqrt(pow(cellHeight * neighborTranslations[k].x, 2) + pow(cellWidth * neighborTranslations[k].y, 2));
					deltaH = getHeight(currentIJ.x, currentIJ.y) - getHeight(neighborIJ.x, neighborIJ.y);
					flowOut[currentXY.x][currentXY.y][k] = max(0, getFlow(currentIJ.x, currentIJ.y, k) + dt * A * g * deltaH / distanceToNextCell); //Atention au clampin dans flowOutR
					flowOutTot += flowOut[currentXY.x][currentXY.y][k];
				}

				//Calcule du coeff de normalisation K, pour eviter de perdre plus de matiere que la colonne courante du layer
				K = min(1, ThicknessR[getIndex(currentIJ.x, currentIJ.y)][indexOfWater] * A / (flowOutTot * dt));

				//Applique normalisation
				for (uint k = 0; k < 4; k += 1) {
					flowOut[currentXY.x][currentXY.y][k] = K * flowOut[currentXY.x][currentXY.y][k];
				}
			}
		}
	}

	//Retire Phase
	barrier();
	memoryBarrierShared();

	if (pointIJ.x < gridHeight && pointIJ.y < gridWidth) {
	//On update la hauteur et la veloicty que si le point est dans la grille d'ou le if
		//Phase 2 : calcul de ce qui entre, ce qui sort
		currentXY = pointXY + offset; //Localisation de pointXY dans flowOut
		float volumeOutTot = 0.;
		float volumeInTot = 0.;
		uint kN = 0;
		for (uint k = 0; k < 4; k += 1) {
			volumeOutTot += dt * flowOut[currentXY.x][currentXY.y][k];
			
			neighborXY = currentXY + neighborTranslations[k];
			kN = 3 - k;
			volumeInTot += dt * flowOut[neighborXY.x][neighborXY.y][kN];
		}

		//Phase 3 : ecriture avec convention que l'érodé devient sable, besoin de separer en R et W a cause des frontieres/bords des patchs qui se trouve dans plusieurs working group
		//La hauteur de l'eau update
		float d2 = ThicknessR[getIndex(pointIJ.x, pointIJ.y)][indexOfWater] + (volumeInTot - volumeOutTot) / A;
		float dMean = (d2 + ThicknessR[getIndex(pointIJ.x, pointIJ.y)][indexOfWater]) / 2.;
		//Vitesse de l'eau dans la direction de i, X = u
		float deltaWI = (flowOut[currentXY.x - 1][currentXY.y][3] - flowOut[currentXY.x][currentXY.y][0] + flowOut[currentXY.x][currentXY.y][3] - flowOut[currentXY.x + 1][currentXY.y][0]) / 2.;
		float u = deltaWI / (cellWidth * dMean);
		//Vitesse de l'eau dans la direction de j, Y = v
		float deltaWJ = (flowOut[currentXY.x][currentXY.y - 1][2] - flowOut[currentXY.x][currentXY.y][1] + flowOut[currentXY.x][currentXY.y][2] - flowOut[currentXY.x][currentXY.y + 1][1]) / 2.;
		float v = deltaWJ / (cellHeight * dMean);


		//Phase 4 erosion et deposition...
		//To Do


		//Phase 5, ecriture des outputs
		for (uint k = 0; k < 4; k += 1) {
			FlowOutW[getIndex(pointIJ.x, pointIJ.y)][k] = flowOut[currentXY.x][currentXY.y][k];
		}
		for (uint k = 0; k < NB_OF_LAYERS; k += 1) {
			ThicknessW[getIndex(pointIJ.x, pointIJ.y)][k] = ThicknessR[getIndex(pointIJ.x, pointIJ.y)][k];
			//ThicknessW[getIndex(pointIJ.x, pointIJ.y)][k] = 0.5;
		}
		ThicknessW[getIndex(pointIJ.x, pointIJ.y)][indexOfWater] = d2;
		Velocity[getIndex(pointIJ.x, pointIJ.y)][0] = u;
		Velocity[getIndex(pointIJ.x, pointIJ.y)][1] = v;
	}
}