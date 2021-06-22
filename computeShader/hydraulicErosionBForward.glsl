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
uniform float cellHeight, cellWidth, dt;

const ivec2 neighborTranslations[9] = ivec2[9](
	ivec2(-1, -1),
	ivec2(-1, 0),
	ivec2(-1, 1),
	ivec2(0, -1),
	ivec2(0, 0), //Ici on concidere le point dans les voisins qui il se transporte des sendiments sur lui aussi
	ivec2(0, 1),
	ivec2(1, -1),
	ivec2(1, 0),
	ivec2(1, 1));

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

vec2 getScaledVelocity(int i, int j) { //Si en dehors de la grille pas d'eau, donc vitesse nulle, aussi on scale la vitesse pour que dans [0, 1]^2 une fois mult par dt, car utilise que pour repartition dans voisins
	if (i < 0 || j < 0) return vec2(0., 0.);
	if (i > gridHeight - 1 || j > gridWidth - 1) return vec2(0., 0.);
	return vec2(Velocity[getIndex(i, j)][0] / cellHeight, Velocity[getIndex(i, j)][1] / cellWidth); //Vitesse en I, vitesse en J scaled, on prefere grille carre donc ???
}

void main() {
	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);

	if (pointIJ.x < gridHeight && pointIJ.y < gridWidth) {
		//On calcule la qqt de sediment qu'il recoit par ses 9 voisins, qui se deverse possiblement sur lui
		float outputS = 0;
		vec2 d;//Ce le vecteur qui part de pointIJ et qui va vers le point pointé par la vitesse * dt du voisin courant, auquel on applique la fonction abs coord par coord
		ivec2 neighborIJ;
		for (uint k = 0; k < 9; k += 1) {
			neighborIJ = pointIJ + neighborTranslations[k];
			d = abs(neighborIJ + dt * getScaledVelocity(neighborIJ.x, neighborIJ.y) - pointIJ);
			//Test si d suffisamment proche de pointIJ i.e. si dans son voisinage proche tel que pointIJ recoive de la matiere
			if (d.x < 1. && d.y < 1.) {
				outputS += (1 - d.x) * (1 - d.y)  * getSediment(neighborIJ.x, neighborIJ.y);
			}
		}

		SedimentW[getIndex(pointIJ.x, pointIJ.y)] = outputS;
	}
}