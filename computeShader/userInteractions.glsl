layout(std430, binding = 0) buffer ThickR {
	float ThicknessR[][NB_OF_LAYERS];
};

uniform float gridHeight, gridWidth; //Pour la division
uniform vec2 cursorPos;
uniform sampler2D gridScreenPosTex;
uniform bool addMaterial;
uniform bool removeMaterial;

const int indexOfSand = NB_OF_LAYERS - 2;

layout(local_size_x = PATCH_HEIGHT, local_size_y = PATCH_WIDTH, local_size_z = 1) in;

uint getIndex(int i, int j) {
	//Ne doivent pas être en dehors de la grille !!
	return int(gridWidth) * i + j; //Conversation Okay ??
}

float ker(float d, float r) {
	//On note r le rayon du disque sur lequel le kernet a de l'effet.
	if (d > r) return 0.0;

	float sig = r / 3;
	return 0.01 * exp(d * d / (-2 * sig * sig));
}

void main() {
	const ivec2 pointIJ = ivec2(gl_GlobalInvocationID);
	const float iN = float(pointIJ.x) / gridHeight;
	const float jN = float(pointIJ.y) / gridWidth;


	if (iN < 1 && jN < 1) {
		float newH = ThicknessR[getIndex(pointIJ.x, pointIJ.y)][indexOfSand];
		vec2 cursorGridPos = texture(gridScreenPosTex, vec2(cursorPos.y, 1.0 - cursorPos.x)).xy;

		if (cursorGridPos.x >= 0.0 && cursorGridPos.y >= 0.0) {
			float d = distance(cursorGridPos, vec2(iN, jN));

			if (addMaterial) {
				newH = newH + ker(d, 0.05);
			}
			if (removeMaterial) {
				newH = newH - ker(d, 0.05);
			}

		}

		ThicknessR[getIndex(pointIJ.x, pointIJ.y)][indexOfSand] = max(newH, 0.0);
	}
}