// ----------------------------------------------------------------------------
// main.cpp
//
//  Created on: 24 Jul 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: IGR201 Practical; OpenGL and Shaders (DO NOT distribute!)
//
// Copyright 2020 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/ext.hpp>

#include <imgui.h>
#include <imgui_impl_opengl3.h>
#include <imgui_impl_glfw.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>

#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define PI 3.141592


int g_layer = 1;
float g_h_min = 0.f;
float g_h_max = 2.f;
float g_dy = 0.5f;

// Window parameters
GLFWwindow *g_window = nullptr;
int g_windowWidth = 1024;
int g_windowHeight = 768;

//Simulation parameters

//Thermal erosion
float g_thetaLimit_t = 0.3;
float g_erosionCoeff_t = 0.3;
float g_dt_t = 0.0500;
int g_iter_t = 500;
int g_strategyErosion_t = 0; //stratégie A
int g_typeErosion_t = 1; //Layer breaks into sand
int g_connexity_t = 0; // 4 connexité
int g_descentDirection_t = 0; //Heighest gradient
int g_typeGradient_t = 0; //Heighest gradient
int g_neighbourReceiver_t = 0; //all neighbors
unsigned int g_nbOfIterations_t = 0;
unsigned int g_nbOfIterations_t_max = 1;

//Hydraulic erosion
float g_dt_h = 0.0500;
int g_iter_h = 500;
unsigned int g_nbOfIterations_h = 0;
unsigned int g_nbOfIterations_h_max = 1;
bool g_waterarriving = true;
float g_water_flow = 10.f;
int g_water_coords[2] = {50, 50};
float g_k_h[3] = {0.03f, 0.03f, 0.03f};

//Fault parameters
int g_fault_mode = 2;
int g_fault_niter = 1;
float g_fault_dy = 0.5f;
int g_fault_circlemode = 0;
unsigned int g_fault_nbOfIterations = 0;
unsigned int g_fault_nbOfIterations_max = 1;

//Perlin parameters
unsigned int g_perl_nbOfIterations = 0;
unsigned int g_perl_nbOfIterations_max = 1;
int g_perl_niter = 1;
int g_noise_type = 0; //0 : random, 1 : perlin
int g_perl_layer = 1;
int g_max_height_perlin = 2;
int g_min_height_perlin = 0;
float g_max_noise = g_h_max*3/2;
int g_resolutionX = 1;
int g_resolutionY = 1;
std::vector<float> noiseVector;
std::vector<glm::vec2> randomGradients; //liste de gradients random de norme 1
int g_number_octaves = 10;

GLFWwindow* g_window2 = nullptr;

// GPU objects
GLuint g_program = 0; // A GPU program contains at least a vertex shader and a fragment shader


// Basic camera model
class Camera {
public:
    const glm::vec3& getPosition() const { return m_pos; }
    void setPosition(const glm::vec3& t) { m_pos = t; }
    const glm::vec3& getRotation() const { return m_rotation; }
    void setRotation(const glm::vec3& r) { m_rotation = r; }
    inline float getFov() const { return m_fov; }
    inline void setFoV(const float f) { m_fov = f; }
    inline float getAspectRatio() const { return m_aspectRatio; }
    inline void setAspectRatio(const float a) { m_aspectRatio = a; }
    inline float getNear() const { return m_near; }
    inline void setNear(const float n) { m_near = n; }
    inline float getFar() const { return m_far; }
    inline void setFar(const float n) { m_far = n; }

    glm::mat4 computeViewMatrix() const {
        glm::mat4 rot = glm::rotate(glm::mat4(1.0), m_rotation[0], glm::vec3(1.0, 0.0, 0.0));
        rot = glm::rotate(rot, m_rotation[1], glm::vec3(0.0, 1.0, 0.0));
        rot = glm::rotate(rot, m_rotation[2], glm::vec3(0.0, 0.0, 1.0));
        const glm::mat4 trn = glm::translate(glm::mat4(1.0), m_pos);
        return glm::inverse(rot * trn);
    }

    // Returns the projection matrix stemming from the camera intrinsic parameter.
    glm::mat4 computeProjectionMatrix() const {
        return glm::perspective(glm::radians(m_fov), m_aspectRatio, m_near, m_far);
    }

private:
    glm::vec3 m_pos = glm::vec3(0, 0, -10);
    glm::vec3 m_rotation = glm::vec3(0, 0, 0);

    float m_fov = 45.f;        // Field of view, in degrees
    float m_aspectRatio = 1.f; // Ratio between the width and the height of the image
    float m_near = 0.1f; // Distance before which geometry is excluded fromt he rasterization process
    float m_far = 4.f; // Distance after which the geometry is excluded fromt he rasterization process
};
Camera g_camera;
// camera control variables
float g_meshScale = 1.0; // to update based on the mesh size, so that navigation runs at scale
bool g_rotatingP = false;
bool g_panningP = false;
bool g_zoomingP = false;
double g_baseX = 0.0, g_baseY = 0.0;
glm::vec3 g_baseTrans(0.0);
glm::vec3 g_baseRot(0.0);

//Mesh
class Mesh {
public:
    Mesh(const std::vector<std::string>& filenames, const std::vector<glm::vec4> layersColor, const glm::vec4& corners, const glm::vec2& e);
    Mesh(const int nbOfLayers, const int width, const int height, const int mode, const std::vector<glm::vec3> layersColor, const glm::vec4& corners, const glm::vec2& e);

    void init(bool withWater); //Génère la surface à partir des différentes épaisseur de niveau, ces normales, puis envoie l'info au GPU
    void render();
    unsigned int getIndex(unsigned int k, unsigned int i, unsigned int j);
    float getH(int i, int j) const; //Donne la hauteur du terrain eau comprise, autorise indice negatif grace au clamping
    float getTerrainH(int i, int j) const; //Pour avoir la hauteur issue des différentes epaisseurs des layers au point (i,j) de la grille
    float getLayerThickness(unsigned int k, int i, int j) const;
    float getLayerH(unsigned int k, int i, int j) const;
    void setLayerThickness(float value, unsigned int k, unsigned int i, unsigned int j);
    unsigned int getTopLayerId(unsigned int i, unsigned int j) const; //Id de 0 à m_nbOfLayers - 1 correspond à indice dans m_layersColor
    unsigned int getTerrainTopLayerId(unsigned int i, unsigned int j) const; //Id de 0 à m_nbOfLayers - 2 correspond on exclu l'eau
    glm::uvec2 getLowestNeighbor(int i, int j) const;
    float getIDerivate(int i, int j, bool withWater) const; //Dervive par rapport à i c'est à dire quand passe de ligne i à ligne i + 1 (Z)
    float getJDerivate(int i, int j, bool withWater) const; //Derive par rapport à j (X)
    glm::vec2 getGradient(int i, int j, bool withWater) const;
    void setLayersColors(int layer, float color[]);
    std::vector<glm::vec4> getLayersColors();
    std::vector<glm::uvec2> getAllLowNeighbors(int i, int j, bool connexity8) const;
    void thermalErosionA(float thetaLimit,float erosionCoeff,float dt, bool neighbourReceiver, bool descentDirection, bool typeErosion, bool connexity8);
    void thermalErosionB(float thetaLimit, float erosionCoeff, float dt, bool connexity8);
    void hydraulicErosion(unsigned int N, float dt, float k[3]);
    //float Mesh::getLayersFlux(unsigned int k, unsigned int i, unsigned int j);
    void setLayersFlux(float value, unsigned int k, unsigned int i, unsigned int j);
    glm::vec4 getLayersFlux(int i, int j);
    void setLayersFlux(glm::vec4 values, unsigned int i, unsigned int j);
    glm::vec2 getLayersVelocity(int i, int j); 
    void setLayersVelocity(glm::vec2 values, unsigned int i, unsigned int j);
    float getLayersSediment(int i, int j);
    void setLayersSediment(float values, unsigned int i, unsigned int j);
    float waterArriving(unsigned int i, unsigned int j) const;
    void applyNThermalErosion(unsigned int N, float thetaLimit, float erosionCoeff, float dt, bool neighbourReceiver, bool descentDirection, bool typeErosion, bool connexity8, bool strategyB);
    void applyFault(const int& mode, const int& n_iter, const float& _dy );
    std::vector<std::string> getLayersFileNames();
    void applyNoise(const int& n_iter, int layerID, float maxNoise, int typeNoise, int numberOfOctaves);
    void createGridPerlinNoise();
    float perlinNoise(float x, float y,int resolutionX, int resolutionY);
    float fractionalBrownianMotion(float x, float y, int numberOfOctaves);


private:
    unsigned int m_gridWidth = 0; //Nb de colonnes de la grille de discretisation (image)
    unsigned int m_gridHeight = 0; //Nb de lignes de la grille de discretisation (image)
    unsigned int m_nbOfLayers = 0;
    glm::vec2 m_gridTopLeftCorner;
    glm::vec2 m_gridBottomRightCorner;
    float m_cellWidth = 0;
    float m_cellHeight = 0;
    std::vector<std::string> m_layersFileNames;
    std::vector<float> m_layersThickness;
    std::vector<glm::vec4> m_layersColor;
    std::vector<glm::vec4> m_layersFlux; //Top, Left, Right, Bottom
    std::vector<glm::vec2> m_layersVelocity;
    std::vector<float> m_layersSediment;

    std::vector<float> m_vertexPositions;
    std::vector<float> m_vertexNormals;
    std::vector<float> m_vertexColors;
    std::vector<unsigned int> m_triangleIndices;
    std::vector<float> m_vertexTexCoords;
    GLuint m_texCoordsVbo = 0;
    GLuint m_vao = 0;
    GLuint m_posVbo = 0;
    GLuint m_normalVbo = 0;
    GLuint m_ibo = 0;
    GLuint m_colVbo = 0;
    //unsigned char* Mesh::createHeightMapFault(const int& width, const int& height, const int& mode, const int& n_iter);
    unsigned char* loadHeightMapFromFile(const std::string& filename, int& width, int& height, int& channels);
};

//Mesh CPU
Mesh* mesh;

void Mesh::init(bool withWater) {
    //Update cpu
    float ax = m_gridTopLeftCorner.x; //a coin en haut gauche
    float az = m_gridTopLeftCorner.y;
    float bx = m_gridBottomRightCorner.x; //b coin en bas droite
    float bz = m_gridBottomRightCorner.y;
    glm::vec2 grad(0.f, 0.f);
    glm::vec3 normal(0.f, 0.f, 0.f);
    glm::vec4 color(0.f, 0.f, 0.f, 0.f);
    unsigned int ind = 0;
    unsigned int ind_c = 0;
    float h = 0.f;
    for (int i = 0; i < m_gridHeight; i++) {
        for (int j = 0; j < m_gridWidth; j++) {//On calcule les vecteurs normaux, en utilisant le gradient de la fonction d'élévation
            grad = getGradient(i, j, withWater); //Met a true car affiche l'eau
            normal = glm::normalize(glm::vec3(-grad.y, 1, -grad.x)); //On fait attention bien mettre dans bon ordre i.e. derive par rapport à x, correspond à derive par rapport à j...
            color = m_layersColor[getTerrainTopLayerId(i, j)];
            if (withWater) {
                color = m_layersColor[getTopLayerId(i, j)];
            }
            //glm::vec2 gradN = (grad.x >0 && grad.y > 0) ? glm::normalize(grad) : grad;
            //color.x = 1 / 2.f * (gradN.x + 1);
            //color.y = 1 / 2.f * (gradN.y + 1);
            //color.z = 1 / 2.f;
            //color = glm::vec3(0.5, 0.5, 0.5);
            //color = m_layersColor[getTopLayerId(i, j)];

            m_vertexPositions[ind] = (ax + (bx - ax) * j / (m_gridWidth - 1)); //x
            m_vertexNormals[ind] = normal.x;
            m_vertexColors[ind_c] = color.x;
            ind += 1;
            ind_c += 1;

            // = getTerrainH(i, j); //y
            if (withWater) {
                h = getH(i, j);
                if (getTopLayerId(i, j) < m_nbOfLayers - 1) {
                    h -= 0.001f;
                }
            } else {
                h = getTerrainH(i, j);
            }

            m_vertexPositions[ind] = h;
            m_vertexNormals[ind] = normal.y;
            m_vertexColors[ind_c] = color.y;
            ind += 1;
            ind_c += 1;

            m_vertexPositions[ind] = (az + (bz - az) * i / (m_gridWidth - 1)); //z
            m_vertexNormals[ind] = normal.z;
            m_vertexColors[ind_c] = color.z;
            ind += 1;
            ind_c += 1;

            m_vertexColors[ind_c] = color.w;
            ind_c += 1;
        }
    }

    size_t vertexBufferSize = sizeof(float) * m_vertexPositions.size(); // Gather the size of the buffer from the CPU-side vector
    //Update buffer position
    glNamedBufferSubData(m_posVbo, 0, vertexBufferSize, m_vertexPositions.data());
    //Update buffer des normales
    glNamedBufferSubData(m_normalVbo, 0, vertexBufferSize, m_vertexNormals.data());
    //Update buffer des couleur
    vertexBufferSize = sizeof(float) * m_vertexColors.size();
    glNamedBufferSubData(m_colVbo, 0, vertexBufferSize, m_vertexColors.data());

    //Buffer for texture
    //size_t vertexBufferSize3 = sizeof(float) * m_vertexTexCoords.size(); // Gather the size of the buffer from the CPU-side vector
    //glCreateBuffers(1, &m_texCoordsVbo);
    //glNamedBufferStorage(m_texCoordsVbo, vertexBufferSize3, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    //glNamedBufferSubData(m_texCoordsVbo, 0, vertexBufferSize3, m_vertexTexCoords.data()); // Fill the data storage from a CPU array
    //glBindBuffer(GL_ARRAY_BUFFER, m_texCoordsVbo);
    //glEnableVertexAttribArray(3);
    //glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), 0);
}

void Mesh::setLayersColors(int layer, float color[]) {
    std::cout << color[0] << " " << color[1] << std::endl;
    m_layersColor[layer] = glm::vec4(color[0], color[1], color[2], color[3]);

}

std::vector<std::string> Mesh::getLayersFileNames() {
    return m_layersFileNames;
};

std::vector<glm::vec4> Mesh::getLayersColors() {
    return m_layersColor;
};


void Mesh::render() {
    const glm::vec3 camPosition = g_camera.getPosition();
    glUniform3f(glGetUniformLocation(g_program, "camPos"), camPosition[0], camPosition[1], camPosition[2]);

    glm::mat4 model = glm::mat4(1.f); //matrice modele
    GLuint g_textureID = 0; 

    //On donne la matrice du modèle au vertex shader et les couleurs au fragment shader
    glUniformMatrix4fv(glGetUniformLocation(g_program, "modelMat"), 1, GL_FALSE, glm::value_ptr(model)); 
    
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(glGetUniformLocation(g_program, "material.albedoTex"), 0);
    glBindTexture(GL_TEXTURE_2D, g_textureID);
    
    glBindVertexArray(m_vao);     // bind the VAO storing geometry data
    glDrawElements(GL_TRIANGLES, m_triangleIndices.size(), GL_UNSIGNED_INT, 0); // Call for rendering: stream the current GPU geometry through the current GPU program
}

unsigned char* Mesh::loadHeightMapFromFile(const std::string& filename, int& width, int& height, int& channels) {
    
    unsigned char* heightMap = stbi_load(
        filename.c_str(),
        &width, &height,
        &channels, // 1 for a 8 bit greyscale image, 3 for 24bits RGB image, 4 for 32bits RGBA image
        0);

    if (heightMap == NULL) {
        printf("Error in loading the height map image.\n");
        exit(1);
    }
    else {
        printf("Loaded image with a width of %dpx, a height of %dpx and %d channels\n", width, height, channels);
    }

    size_t img_size = width * height * channels;
    int gray_channels = channels == 4 ? 2 : 1;
    size_t gray_img_size = width * height * gray_channels;

    unsigned char * gray_img = (unsigned char * ) malloc(gray_img_size);

    if (gray_img == NULL) {
        printf("Unable to allocate memory for the gray image.\n");
        exit(1);
    }

    for (unsigned char* p = heightMap, *pg = gray_img; p != heightMap + img_size; p += channels, pg += gray_channels) {
        //To gray 
        *pg = (uint8_t)(*p * 0.39 + *(p + 1) * 0.50 + *(p + 2) * 0.11);

        if (channels == 4) {
            *(pg + 1) = *(p + 3);
        }
    }

    stbi_image_free(heightMap);

    return gray_img;
}


Mesh::Mesh(const std::vector<std::string>& filenames, const std::vector<glm::vec4> layersColor, const glm::vec4& corners, const glm::vec2& e) {
    int width = 0, height = 0, channels = 0;
    m_nbOfLayers = filenames.size();
    m_layersColor = layersColor;
    m_gridTopLeftCorner = glm::vec2(corners.x, corners.y);
    m_gridBottomRightCorner = glm::vec2(corners.z, corners.w);
    float emin = e.x;
    float emax = e.y;

    m_layersFileNames = filenames;

    for (unsigned int k = 0; k < m_nbOfLayers; k = k + 1) {

        unsigned char* gray_img = loadHeightMapFromFile(filenames[k], width, height, channels);

        //Met a jour la taille de la grille avec la valeur de la première map
        if (k == 0) {
            m_gridWidth = width;
            m_gridHeight = height;

            m_cellWidth = abs(m_gridBottomRightCorner.x - m_gridTopLeftCorner.x) / (width - 1.f);
            m_cellHeight = abs(m_gridBottomRightCorner.y - m_gridTopLeftCorner.y) / (height - 1.f);

            m_layersThickness.resize(m_nbOfLayers * m_gridWidth * m_gridHeight);
            m_layersFlux.resize(4 * m_gridWidth * m_gridHeight);
            m_layersVelocity.resize(2 * m_gridWidth * m_gridHeight);
            m_layersSediment.resize(m_gridWidth * m_gridHeight);
            m_vertexPositions.resize(3 * m_gridWidth * m_gridHeight);
            m_vertexNormals.resize(3 * m_gridWidth * m_gridHeight);
            m_vertexColors.resize(4 * m_gridWidth * m_gridHeight);

            for (int i = 0; i < m_gridWidth * m_gridHeight - m_gridWidth; i++) {
                if (i % m_gridWidth != m_gridWidth - 1) {
                    m_triangleIndices.push_back(i);
                    m_triangleIndices.push_back(i + m_gridWidth);
                    m_triangleIndices.push_back(i + 1);
                    m_triangleIndices.push_back(i + 1);
                    m_triangleIndices.push_back(i + m_gridWidth);
                    m_triangleIndices.push_back(i + m_gridWidth + 1);
                }
            }

            for (int i = 0; i < height; i++) {
                for (int j = 0; j < width; j++) {
                    m_layersFlux[j + width * i] = glm::vec4(0.f, 0.f, 0.f, 0.f);
                    m_layersVelocity[j + width * i] = glm::vec2(0.f, 0.f);
                }
            }
        }

        //Verfier que bonne dim
        if (width != m_gridWidth || height != m_gridHeight) {
            std::cout << "ERROR : wrong size" << std::endl;
            break;
        }

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int src_index = j + width * i; //Image déroule ligne par ligne numComponents * j + numComponents * width< * i (au cas ou)
                setLayerThickness(gray_img[src_index] / 255.f * (emax - emin) + emin, k, i, j);
                
            }
        }

        free(gray_img);
    }

    //Pour les coordonnées des textures, on avance de 1/resolution pour remplir avec la texture entre 0 et 1
    //m_vertexTexCoords = {};
    /*for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            //m_vertexTexCoords.push_back(float(j) / (width - 1)); //x
            //m_vertexTexCoords.push_back(float(i) / (height - 1)); //y
        }
    }*/

    //Creation des buffers et création du vao
    //Create a single handle that joins together attributes (vertex positions,
    //normals) and connectivity (triangles indices)
    glCreateVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    // Generate a GPU buffer to store the positions of the vertices
    size_t vertexBufferSize = sizeof(float) * m_vertexPositions.size(); // Gather the size of the buffer from the CPU-side vector
    glCreateBuffers(1, &m_posVbo);
    glNamedBufferStorage(m_posVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

    //Buffer for normal
    glCreateBuffers(1, &m_normalVbo);
    glNamedBufferStorage(m_normalVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    glBindBuffer(GL_ARRAY_BUFFER, m_normalVbo);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

    //Buffer for color
    vertexBufferSize = sizeof(float) * m_vertexColors.size();
    glCreateBuffers(1, &m_colVbo);
    glNamedBufferStorage(m_colVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    glBindBuffer(GL_ARRAY_BUFFER, m_colVbo);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    //Buffer for texture
    //size_t vertexBufferSize3 = sizeof(float) * m_vertexTexCoords.size(); // Gather the size of the buffer from the CPU-side vector
    //glCreateBuffers(1, &m_texCoordsVbo);
    //glNamedBufferStorage(m_texCoordsVbo, vertexBufferSize3, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    //glBindBuffer(GL_ARRAY_BUFFER, m_texCoordsVbo);
    //glEnableVertexAttribArray(3);
    //glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), 0);

    // Same for the index buffer that stores the list of indices of the
    // triangles forming the mesh
    size_t indexBufferSize = sizeof(unsigned int) * m_triangleIndices.size();
    glCreateBuffers(1, &m_ibo);
    glNamedBufferStorage(m_ibo, indexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT);
    glNamedBufferSubData(m_ibo, 0, indexBufferSize, m_triangleIndices.data());
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo); // bind the IBO storing geometry data

    glBindVertexArray(0);
    glBindVertexArray(0);
}

/*Mesh::Mesh(const int nbOfLayers, const int width, const int height, const int mode, const std::vector<glm::vec3> layersColor, const glm::vec4& corners, const glm::vec2& e) {
    m_nbOfLayers = nbOfLayers;
    m_layersColor = layersColor;
    m_gridTopLeftCorner = glm::vec2(corners.x, corners.y);
    m_gridBottomRightCorner = glm::vec2(corners.z, corners.w);
    float emin = e.x;
    float emax = e.y;

    for (unsigned int k = 0; k < m_nbOfLayers; k = k + 1) {
        unsigned char* gray_img = createHeightMapFault(width, height, g_fault_mode, g_fault_niter);

        //Met a jour la taille de la grille avec la valeur de la première map
        if (k == 0) {
            m_gridWidth = width;
            m_gridHeight = height;

            m_cellWidth = abs(m_gridBottomRightCorner.x - m_gridTopLeftCorner.x) / (width - 1.f);
            m_cellHeight = abs(m_gridBottomRightCorner.y - m_gridTopLeftCorner.y) / (height - 1.f);

            //std::cout << m_cellHeight << " " << m_cellHeight << std::endl;

            m_layersThickness.resize(m_nbOfLayers * m_gridWidth * m_gridHeight);
            m_vertexPositions.resize(3 * m_gridWidth * m_gridHeight);
            m_vertexNormals.resize(3 * m_gridWidth * m_gridHeight);
            m_vertexColors.resize(3 * m_gridWidth * m_gridHeight);

            for (int i = 0; i < m_gridWidth * m_gridHeight - m_gridWidth; i++) {
                if (i % m_gridWidth != m_gridWidth - 1) {
                    m_triangleIndices.push_back(i);
                    m_triangleIndices.push_back(i + m_gridWidth);
                    m_triangleIndices.push_back(i + 1);
                    m_triangleIndices.push_back(i + 1);
                    m_triangleIndices.push_back(i + m_gridWidth);
                    m_triangleIndices.push_back(i + m_gridWidth + 1);
                }
            }
        }

        //Verfier que bonne dim
        if (width != m_gridWidth || height != m_gridHeight) {
            std::cout << "ERROR : wrong size" << std::endl;
            break;
        }

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                int src_index = j + width * i; //Image déroule ligne par ligne numComponents * j + numComponents * width< * i (au cas ou)
                setLayerThickness(gray_img[src_index] / 255.f * (emax - emin) + emin, k, i, j);
            }
        }

        free(gray_img);
    }

    //Pour les coordonnées des textures, on avance de 1/resolution pour remplir avec la texture entre 0 et 1
    //m_vertexTexCoords = {};
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            m_vertexTexCoords.push_back(float(j)); //x
            m_vertexTexCoords.push_back(float(i)); //y
        }
    }
}*/

unsigned int Mesh::getTopLayerId(unsigned int i, unsigned int j) const {
    unsigned int id = m_nbOfLayers - 1;

    while (getLayerThickness(id, i, j) == 0.f && id > 0) {
        id = id - 1;
    }

    return id;
}

unsigned int Mesh::getTerrainTopLayerId(unsigned int i, unsigned int j) const {
    unsigned int id = m_nbOfLayers - 2;

    while (getLayerThickness(id, i, j) == 0.f && id > 0) {
        id = id - 1;
    }

    return id;
}

unsigned int Mesh::getIndex(unsigned int k, unsigned int i, unsigned int j) {
    return k * m_gridHeight * m_gridWidth + i * m_gridWidth + j;
}

float Mesh::getLayerH(unsigned int k, int i, int j) const {
    float h = 0;
    for (unsigned int l = 0; l <= k; l = l + 1) {
        h += getLayerThickness(l, i, j);
    }
    return h;
}

glm::uvec2 Mesh::getLowestNeighbor(int i, int j) const {
    //Prend en compte le pixel central donc pas besoin de vérifier qu'on a bien un voisin plus bas que ce pixel central

    //Test en connexité 8
    int lowestDI = -1;
    int lowestDJ = -1;

    float lowestH = getTerrainH(i + lowestDI, j + lowestDJ);
    float currentH = 0;

    for (int di = -1; di < 2; di += 1) {
        for (int dj = -1; dj < 2; dj += 1) {
            currentH = getTerrainH(i + di, j + dj);
            if (currentH < lowestH) {
                lowestH = currentH;
                lowestDI = di;
                lowestDJ = dj;
            }
        }
    }

    return glm::uvec2(i, j) + glm::uvec2(lowestDI, lowestDJ);
}

std::vector<glm::uvec2> Mesh::getAllLowNeighbors(int i, int j, bool connexity8) const {
    //Test en connexité 8, retourne le vecteur de positions des voisins dans l'ordre des plus bas au plus haut 
    // (jusqu'à la hauteur du pixel central)

    float currentH;
    float hCenter = getTerrainH(i, j);
    std::vector<glm::uvec2> vectorOfLowNeighbors;
    

    if (connexity8) {
        for (int di = -1; di < 2; di += 1) {
            for (int dj = -1; dj < 2; dj += 1) {
                int nextI = i + di;
                int nextJ = j + dj;
                //on vérifie qu'on ne sort pas de la grille
                if (nextI >= 0 && nextJ >= 0 && nextI < m_gridHeight && nextJ < m_gridWidth) {
                    currentH = getTerrainH(nextI, nextJ);
                    if (currentH < hCenter) {
                        vectorOfLowNeighbors.push_back(glm::uvec2(nextI, nextJ));
                    }
                }

            }
        }
    }
    // 4 connexité
    else {
        for (int di = -1; di < 2; di += 1) {
            for (int dj = -1; dj < 2; dj += 1) {
                //on est en connexité 4
                if (di != 1 || dj != 1) {

                    int nextI = i + di;
                    int nextJ = j + dj;
                    //on vérifie qu'on ne sort pas de la grille
                    if (nextI >= 0 && nextJ >= 0 && nextI < m_gridHeight && nextJ < m_gridWidth) {
                        currentH = getTerrainH(nextI, nextJ);
                        if (currentH < hCenter) {
                            vectorOfLowNeighbors.push_back(glm::uvec2(nextI, nextJ));
                        }
                    }

                }

            }
        }
    }

    if (vectorOfLowNeighbors.size() == 0) {
        //on compte le pixel central car on lui a enlevé de la matière donc s'il
        //n'a aucun voisin plus bas que lui, on lui remet sa matière
        vectorOfLowNeighbors.push_back(glm::uvec2(i, j));
    }
    
    return vectorOfLowNeighbors;
}

float Mesh::getLayerThickness(unsigned int k, int i, int j) const {
    //Clamp (i, j) si en dehors de la grille
    if (i < 0) i = 0;
    if (i >= m_gridHeight) i = m_gridHeight - 1;
    if (j < 0) j = 0;
    if (j >= m_gridWidth) j = m_gridWidth - 1;
    return m_layersThickness[k * m_gridHeight * m_gridWidth + i * m_gridWidth + j];
}

void Mesh::setLayerThickness(float value, unsigned int k, unsigned int i, unsigned int j) {
    if (value < 0) {
        value = 0.f;
    }

    m_layersThickness[k * m_gridHeight * m_gridWidth + i * m_gridWidth + j] = value;
}

float Mesh::getH(int i, int j) const {
    //Donne la hauteur en un point avec l'eau comprise (pour hydraulic simulation)
    float h = 0;
    for (unsigned int k = 0; k < m_nbOfLayers; k = k + 1) {
        h += getLayerThickness(k, i, j);
    }
    return h;
}

float Mesh::getTerrainH(int i, int j) const {
    //Donne la hauteur en un point du terrain, dont exclu l'eau qui est par convention dans l'index m_nbOfLayers - 1 (pour thermal erosion)
    float h = 0;
    for (unsigned int k = 0; k < m_nbOfLayers - 1; k = k + 1) {
        h += getLayerThickness(k, i, j);
    }
    return h;
}

float Mesh::getIDerivate(int i, int j, bool withWater) const {
    auto computeH = &Mesh::getTerrainH;
    if (withWater) {
        computeH = &Mesh::getH;
    }

    //On fait attention au bord
    if (i == 0) {
        return ((this->*computeH)(i + 1, j) - (this->*computeH)(i, j)) / m_cellHeight;
    }
    if (i == m_gridHeight - 1) {
        return ((this->*computeH)(i, j) - (this->*computeH)(i - 1, j)) / m_cellHeight;
    }
    
    //return (getTerrainH(i, j) - getTerrainH(i - 1, j)) / m_cellHeight;
    if (g_typeGradient_t == 0) return ((this->*computeH)(i + 1, j) - (this->*computeH)(i - 1, j)) / (2.f * m_cellHeight);

    else if (g_typeGradient_t == 1) return ((this->*computeH)(i, j) - (this->*computeH)(i-1, j)) / m_cellHeight;

    else return ((this->*computeH)(i+1, j) - (this->*computeH)(i, j)) / m_cellHeight;
}

float Mesh::getJDerivate(int i, int j, bool withWater) const {
    auto computeH = &Mesh::getTerrainH;
    if (withWater) {
        computeH = &Mesh::getH;
    }

    //On fait attention au bord
    if (j == 0) {
        return ((this->*computeH)(i, j + 1) - (this->*computeH)(i, j)) / m_cellWidth;
    }
    if (j == m_gridWidth - 1) {
        return ((this->*computeH)(i, j) - (this->*computeH)(i, j - 1)) / m_cellWidth;
    }

    //return (getTerrainH(i, j) - getTerrainH(i, j - 1)) / m_cellWidth;
    if (g_typeGradient_t==0) return ((this->*computeH)(i, j + 1) - (this->*computeH)(i, j - 1)) / (2.f * m_cellWidth);

    else if (g_typeGradient_t==1) return ((this->*computeH)(i, j) - (this->*computeH)(i, j - 1)) / m_cellWidth;

    else return ((this->*computeH)(i, j+1) - (this->*computeH)(i, j)) / m_cellWidth;
}

glm::vec2 Mesh::getGradient(int i, int j, bool withWater) const {
    return glm::vec2(getIDerivate(i, j, withWater), getJDerivate(i, j, withWater));
}


void Mesh::thermalErosionA(float thetaLimit, float erosionCoeff, float dt, bool neighbourReceiver, bool descentDirection, bool typeErosion, bool connexity8) {
    float tangentLimit = glm::tan(thetaLimit);
    std::vector<float> newLayersThickness = m_layersThickness; //vecteur mémoire temporaire
    //neigbourReceiver = 0 : tous les voisins plus bas recoivent dh, 1 : voisin dans la direction de descente
    //descentDirection = 0 : gradient direction, 1 : lowest neighbour
    //typeErosion = 1 : L'érosion donne du sable, 0 : l'érosion donne le même layer (et passe sous le sable)


    for (int i = 0; i < m_gridHeight; i++) {
        for (int j = 0; j < m_gridWidth; j++)
        {
            glm::vec2 directionDescent = -getGradient(i, j, false);
            float slope = glm::length(directionDescent);

            //if (slope > 0) {directionDescent = directionDescent / slope;} //normalise la direction de descente

            //condition d'érosion
            if (slope > tangentLimit) {
                //index du layer qui donne
                int layerIndexCurrentCell = getTerrainTopLayerId(i, j);

                //index du layer qui reçoit
                int newLayerIndex;
                //Si la matière s'érode et donne du sable
                if (typeErosion != 0) {
                    newLayerIndex = m_nbOfLayers - 2;
                }
                //La matière érodée donne la même matière
                else {
                    newLayerIndex = layerIndexCurrentCell;
                }

                float dh = -erosionCoeff * (slope -tangentLimit) * dt; //<=0

                //Si dh est plus grand que l'epaisseur du layer
                if (-dh > getLayerThickness(layerIndexCurrentCell, i, j)) {
                    dh = -getLayerThickness(layerIndexCurrentCell, i, j);
                }

                //on erode si on est pas le layer le plus bas ou de l'eau
                if (layerIndexCurrentCell > 0) {
                    

                    float newThicknessCurrentCell = getLayerThickness(layerIndexCurrentCell, i, j) + dh;

                    newLayersThickness[layerIndexCurrentCell * m_gridHeight * m_gridWidth + i * m_gridWidth + j] = (newThicknessCurrentCell > 0) ? newThicknessCurrentCell : 0.f;
                    
                    //On donne la matière à un voisin (dans la direction de descente)
                    if (neighbourReceiver!=0) {
                        //std::cout << "Un seul voisin reçoit la matière" << std::endl;

                        int nextI;
                        int nextJ;
                        glm::uvec2 nextCell;

                        //Solution 1 : On donne le layer au voisin dans la direction de descente : peut poser pb 
                        if (descentDirection ==0) {
                            //std::cout << "Le voisin qui reçoit la matière est celui dans la direction de plus fort gradient" << std::endl;
                            nextCell = glm::vec2(i, j) + directionDescent; 
                            nextI = round(nextCell.x); //pour obtenir la cellule i,j dans laquelle on atterit
                            nextJ = round(nextCell.y);
                        }

                        //Solution 2 : On donne le layer au voisin le plus bas
                        else {
                            //std::cout << "Le voisin qui reçoit la matière est le plus bas" << std::endl;
                            nextCell = getLowestNeighbor(i, j);
                            nextI = nextCell.x; //pour obtenir la cellule i,j dans laquelle on atterit
                            nextJ = nextCell.y;
                        }

                        //Etude du cas (26, 34) de qui il recoit de la matiere (i == 26 && j == 34) || (i == 27 && j == 34) || (i == 26 && j == 33) || (i == 26 && j == 35) || (i == 27 && j == 35)
                        //if ((i == 61 && j == 62) || (i == 36 && j == 37)) {
                        //    if (getTerrainH(nextI, nextJ) > getTerrainH(i, j)) {
                        //        std::cout << "PB ";
                        //    }
                        //    std::cout << "i "<< i << " j " << j << " Dir - grad " << directionDescent.x << " " << directionDescent.y << -dh << " nextI " << nextI << " NextJ " << nextJ << std::endl;
                        //}

                        //gérer les bords, on ne transfère la matière que sur des cellules à l'intérieur
                        if (nextI >= 0 && nextJ >= 0 && nextI < m_gridHeight && nextJ < m_gridWidth) {
                            //float newThicknessNextCell = getLayerThickness(layerIndexCurrentCell, nextI, nextJ) - dh;
                            newLayersThickness[newLayerIndex * m_gridHeight * m_gridWidth + nextI * m_gridWidth + nextJ] -= dh;
                        }
  
                    }

                    //On distribue la matière dans tous les voisins, plus un voisin est bas plus reçoit de matière
                    //ici, le voisin i de hauteur hi reçoit -dh*(hi-h_j)/(somme de h_i - h_j) (>=0)
                    //ainsi, plus la différence de niveau est grande, plus le voisin reçoit de la matière.
                    else {
                        //std::cout << "Tous les voisins inférieurs reçoivent la matière" << std::endl;

                        std::vector<glm::uvec2> vectorOfLowNeighbors = getAllLowNeighbors(i, j, connexity8);
                        int nextI;
                        int nextJ;
                        float sumOfDifferences = 0;

                        //calcul de somme des h_i-h_j
                        for (unsigned int i = 0; i < vectorOfLowNeighbors.size(); i++)
                        {

                            glm::uvec2 nextCell = vectorOfLowNeighbors.at(i);
                            nextI = nextCell.x; //pour obtenir la cellule i,j dans laquelle on atterit
                            nextJ = nextCell.y;
                            sumOfDifferences += getTerrainH(i,j)-getTerrainH(nextI, nextJ);
                        }


                        //distribution de dh chez les voisins inférieurs
                        for (unsigned int i = 0; i < vectorOfLowNeighbors.size(); i++)
                        {
                            glm::uvec2 nextCell = vectorOfLowNeighbors.at(i);
                            nextI = nextCell.x; //pour obtenir la cellule i,j dans laquelle on atterit
                            nextJ = nextCell.y;
                            
                            if (sumOfDifferences > 0) {
                                //nextI et nextJ respectent déjà les conditions au bord (cf la fonction getAllLowNeighbors)
                                newLayersThickness[newLayerIndex * m_gridHeight * m_gridWidth + nextI * m_gridWidth + nextJ] -= dh * (getTerrainH(i,j)-getTerrainH(nextI, nextJ)) / sumOfDifferences;
                            }
                            //sumOfDifferences vaut 0 si le pixel central n'a pas de voisins donc il reçoit juste toute la matière qu'il a perdu
                            else {
                                newLayersThickness[newLayerIndex * m_gridHeight * m_gridWidth + nextI * m_gridWidth + nextJ] -= dh;
                            }
                            
                        }
                    } 
                }

            }

        }
    }

    m_layersThickness = newLayersThickness;

    //mesh->init();

    
}

void Mesh::thermalErosionB(float thetaLimit, float erosionCoeff, float dt, bool connexity8 = true) {
    int topLayerIndexCurrentCell = 0.f;
    std::vector<float> dHOut; //8 connexite (i-1, j-1) ; (i-1, j) ; (i-1, j+1) ; (i, j-1) ; (i, j+1) ; (i+1, j-1) ; (i+1, j) ; (i+1, j+1) ;
    float dHOutTot = 0.f;
    std::vector<glm::ivec2> neighborTranslations; //translation forme (dI, dJ)

    //En fonction de la connexite change la liste des deplacements pour trouver les voisins
    if (connexity8) {
        dHOut.resize(8);
        neighborTranslations.resize(8);
        neighborTranslations[0] = glm::vec2(-1, -1);
        neighborTranslations[1] = glm::vec2(-1 , 0);
        neighborTranslations[2] = glm::vec2(-1, 1);
        neighborTranslations[3] = glm::vec2(0, -1);
        neighborTranslations[4] = glm::vec2(0, 1);
        neighborTranslations[5] = glm::vec2(1, -1);
        neighborTranslations[6] = glm::vec2(1, 0);
        neighborTranslations[7] = glm::vec2(1, 1);
    }
    else {
        dHOut.resize(4);
        neighborTranslations.resize(4);
        neighborTranslations[0] = glm::vec2(-1, 0);
        neighborTranslations[1] = glm::vec2(0, -1);
        neighborTranslations[2] = glm::vec2(0, 1);
        neighborTranslations[3] = glm::vec2(1, 0);
    }
    
    float K = 1.f; //Coeff de normalisation si le dHOut superieure à H
    int nextCellI = 0;
    int nextCellJ = 0;
    float distanceToNextCell = 0.f;
    float elevationLimit = glm::tan(thetaLimit);
    std::vector<float> newLayersThickness = m_layersThickness;

    //Il faut faire des blouces sur i et j des int et pas unsigned int car i - 1 peut-être négatif quand i = 0 !!
    for (int i = 0; i < m_gridHeight; i++) {
        for (int j = 0; j < m_gridWidth; j++)
        {
            topLayerIndexCurrentCell = getTerrainTopLayerId(i, j);

            if (topLayerIndexCurrentCell > 0) { //Erosion il y a que si le layer afleurent n'est pas de la bedrock
                //On calcule les dH
                dHOutTot = 0.f;
                for (unsigned int k = 0; k < neighborTranslations.size(); k ++) {
                    nextCellI = i + neighborTranslations[k].x;
                    nextCellJ = j + neighborTranslations[k].y;
                    distanceToNextCell = sqrt(pow(m_cellHeight * neighborTranslations[k].x, 2) + pow(m_cellWidth * neighborTranslations[k].y, 2));
                    dHOut[k] = std::max(0.f, erosionCoeff * ((getTerrainH(i, j) - getTerrainH(nextCellI, nextCellJ)) / distanceToNextCell - elevationLimit) * dt);
                    dHOutTot += dHOut[k];
                }

                //Calcule du coeff de normalisation K, pour eviter de perdre plus de matiere que la colonne courante du layer affleurant
                K = std::min(1.f, getLayerThickness(topLayerIndexCurrentCell, i, j) / dHOutTot);

                //Maj les dHOut
                dHOutTot = 0.f;
                for (unsigned int k = 0; k < neighborTranslations.size(); k++) {
                    dHOut[k] *= K;
                    dHOutTot += dHOut[k];
                }

                //Maj des epaisseurs dans la liste intermediaire
                //De la cellule courante on enlève de la matière
                newLayersThickness[getIndex(topLayerIndexCurrentCell, i, j)] -= dHOutTot;
                //Des cellule voisine il faut tester si on ne sort pas de la grille, on ajoute de la matière
                for (unsigned int k = 0; k < neighborTranslations.size(); k++) {
                    nextCellI = i + neighborTranslations[k].x;
                    nextCellJ = j + neighborTranslations[k].y;
                    if (nextCellI >= 0 && nextCellJ >= 0 && nextCellI < m_gridHeight && nextCellJ < m_gridWidth) {
                        newLayersThickness[getIndex(topLayerIndexCurrentCell, nextCellI, nextCellJ)] += dHOut[k];
                    }
                }
            }
        }
    }

    m_layersThickness = newLayersThickness;
    //mesh->init();
}

void Mesh::applyNThermalErosion(unsigned int N, float thetaLimit, float erosionCoeff, float dt, bool neighbourReceiver, bool descentDirection, bool typeErosion, bool connexity8, bool strategyB) {
    if (strategyB) {
        for (unsigned int i = 0; i < N; i += 1) {
            thermalErosionB(thetaLimit, erosionCoeff, dt, connexity8);
        }
    }
    else {
        for (unsigned int i = 0; i < N; i += 1) {
            thermalErosionA(thetaLimit, erosionCoeff, dt, neighbourReceiver, descentDirection, typeErosion, connexity8);
        }
    } 
}

/*float Mesh::getLayersFlux(unsigned k, unsigned int i, unsigned int j) {
    //Clamp (i, j) si en dehors de la grille
    if (i < 0) i = 0;
    if (i >= m_gridHeight) i = m_gridHeight - 1;
    if (j < 0) j = 0;
    if (j >= m_gridWidth) j = m_gridWidth - 1;

    switch (k) {
        case 0: m_layersFlux[i * m_gridWidth + j].x;
        case 1: m_layersFlux[i * m_gridWidth + j].y;
        case 2: m_layersFlux[i * m_gridWidth + j].z;
        case 3: m_layersFlux[i * m_gridWidth + j].w;

        default:
            return 0.f;
    }
}*/

void Mesh::setLayersFlux(float value, unsigned int k, unsigned int i, unsigned int j) {
    if (value < 0.f) {
        value = 0.f;
    }

    switch (k) {
    case 0: 
        m_layersFlux[i * m_gridWidth + j].x = value;
        break;
    case 1: 
        m_layersFlux[i * m_gridWidth + j].y = value;
        break;
    case 2: 
        m_layersFlux[i * m_gridWidth + j].z = value;
        break;
    case 3: 
        m_layersFlux[i * m_gridWidth + j].w = value;
        break;
    }
}

glm::vec4 Mesh::getLayersFlux(int i, int j) {
    //Il ne faut pas clamp les flux mais les mettre à 0 si en dehors de la grille, car sinon les cellules sur le bords recevoivent de l'eau de leur voisin en dehors de la grille
    if (i < 0 || i >= m_gridHeight || j < 0 || j >= m_gridWidth) return glm::vec4(0.f);
    
    return m_layersFlux[i * m_gridWidth + j];
}

void Mesh::setLayersFlux(glm::vec4 values, unsigned int i, unsigned int j) {
    m_layersFlux[i * m_gridWidth + j] = values;
}

float Mesh::waterArriving(unsigned int i, unsigned int j) const {
    if (i < 0 || i >= m_gridHeight || j < 0 || j >= m_gridWidth) return 0.f;

    if (i == g_water_coords[0] && j == g_water_coords[1] && g_waterarriving) return g_water_flow;

    return 0.f;
}

glm::vec2 Mesh::getLayersVelocity(int i, int j) {
    if (i < 0 || i >= m_gridHeight || j < 0 || j >= m_gridWidth) return glm::vec2(0.f);

    return m_layersVelocity[i * m_gridWidth + j];
}

void Mesh::setLayersVelocity(glm::vec2 values, unsigned int i, unsigned int j) {
    m_layersVelocity[i * m_gridWidth + j] = values;
}

float Mesh::getLayersSediment(int i, int j) {
    if (i < 0 || i >= m_gridHeight || j < 0 || j >= m_gridWidth) return 0.f;

    return m_layersSediment[i * m_gridWidth + j];
}

void Mesh::setLayersSediment(float values, unsigned int i, unsigned int j) {
    m_layersSediment[i * m_gridWidth + j] = values;
}

void Mesh::hydraulicErosion(unsigned int N, float dt, float k[3]) {
    //terrain height : b
    //water height : d
    //suspended sediment amount : s
    //outflow flux f
    //velocity v

    float g = 9.81f;
    glm::vec4 flux;
    glm::vec2 vel;
    float b, d, dh, K, dV, dWx, dWy, u, v, d2, d_, C, st, sin_alpha, grad;
    //float kc = 3.f; //Sediment capacity
    //float ks = 0.3f; //Dissolving
    //float kd = 0.3f; //Deposition
    float kc = k[0]; //Sediment capacity
    float ks = k[1]; //Dissolving
    float kd = k[2]; //Deposition
    float A = m_cellWidth * m_cellHeight;
    float dtAg = dt * A * g;


    for (unsigned int k=0; k < N; k++) {
        //Il faut faire des boucles sur i et j des int et pas unsigned int car i - 1 peut-être négatif quand i = 0 !!
        for (int i = 0; i < m_gridHeight; i++) {
            for (int j = 0; j < m_gridWidth; j++) {
                //Water increment
                setLayerThickness(getLayerThickness(m_nbOfLayers - 1, i, j) + dt * waterArriving(i, j), m_nbOfLayers - 1, i, j);

                //Outflow Flux
                //Top
                dh = getH(i, j) - getH(i + 1, j);
                setLayersFlux(std::max(0.f, getLayersFlux(i, j).x + dtAg * dh / m_cellWidth), 0, i, j);
                //Left
                dh = getH(i, j) - getH(i, j - 1);
                setLayersFlux(std::max(0.f, getLayersFlux(i, j).y + dtAg * dh / m_cellHeight), 1, i, j);
                //Right
                dh = getH(i, j) - getH(i, j + 1);
                setLayersFlux(std::max(0.f, getLayersFlux(i, j).z + dtAg * dh / m_cellWidth), 2, i, j);
                //Bottom
                dh = getH(i, j) - getH(i - 1, j);
                setLayersFlux(std::max(0.f, getLayersFlux(i, j).w + dtAg * dh / m_cellHeight), 3, i, j);

                flux = getLayersFlux(i, j);
                K = std::min(1.f, getLayerThickness(m_nbOfLayers - 1, i, j) * m_cellWidth * m_cellHeight / (flux.x + flux.y + flux.z + flux.w) / dt);
                setLayersFlux(K * flux, i, j);
            }
        }

        //Erosion
        for (int i = 0; i < m_gridHeight; i++) {
            for (int j = 0; j < m_gridWidth; j++) {
                flux = getLayersFlux(i, j);
                dV = dt * (getLayersFlux(i - 1, j).x + getLayersFlux(i, j + 1).y + getLayersFlux(i + 1, j).w + getLayersFlux(i, j - 1).z - flux.x - flux.y - flux.z - flux.w);
                d2 = getLayerThickness(m_nbOfLayers - 1, i, j) + dV / A;
                d_ = getLayerThickness(m_nbOfLayers - 1, i, j) + d2;
                setLayerThickness(d2, m_nbOfLayers - 1, i, j);
                //std::cout << d_ << std::endl;

                dWx = getLayersFlux(i, j - 1).y - getLayersFlux(i, j).z + getLayersFlux(i, j).y - getLayersFlux(i, j + 1).z;
                dWy = getLayersFlux(i + 1, j).x - getLayersFlux(i, j).w + getLayersFlux(i, j).x - getLayersFlux(i - 1, j).w;

                //std::cout <<d_<< std::endl;
                d_ = std::max(0.0001f, d_);

                if (d_ <= 0) {
                    u = 0.f;
                    v = 0.f;
                } else {
                    u = dWx / (d_ * m_cellWidth);
                    v = dWy / (d_ * m_cellHeight);
                }



                //std::cout << u << " " << dWx <<  " "<< d_ * m_cellHeight << std::endl;

                //setLayersVelocity(glm::vec2(u, v), i, j);
                //C = kc * glm::length(glm::vec2(u, v));
                grad = glm::length(getGradient(i, j, false));
                sin_alpha = grad / sqrt(grad * grad + 1);
                C = kc * (sin_alpha + 0.01f) * std::max(0.01f, glm::length(glm::vec2(u, v))) * getLayerThickness(m_nbOfLayers - 1, i, j) / (g_h_max - g_h_min);
                //std::cout << getLayersVelocity(i, j).x << std::endl;
                //std::cout << u << " " << v <<  " " << C << std::endl;

                st = getLayersSediment(i, j);

                if (C > st) {
                    setLayerThickness(getLayerThickness(1, i, j) - ks * (C - st), 1, i, j);
                    setLayersSediment(st + ks * (C - st), i, j);
                    //std::cout << u << " " << v << " " << C << std::endl;
                    //std::cout << st << " " << getLayersSediment(i, j) << std::endl;
                    //std::cout << "true";
                } else {
                    setLayerThickness(getLayerThickness(1, i, j) + kd * (st - C), 1, i, j);
                    setLayersSediment(st - kd * (st - C), i, j);
                    //std::cout << C << " " <<st<< " "<<kd * (st - C) << std::endl;
                    //std::cout << "false";
                    //std::cout << st << " " << getLayersSediment(i, j) << " " << i << " " << j << std::endl;

                }
                //std::cout << C << " " << st << std::endl;
            }
        }

        for (int i = 0; i < m_gridHeight; i++) {
            for (int j = 0; j < m_gridWidth; j++) {
                vel = getLayersVelocity(i, j);
                setLayersSediment(getLayersSediment((int) (i + vel.x * dt), (int) (j - vel.y * dt)), i, j);
            }
        }
    }
    //init();
}

void Mesh::applyFault(const int& mode, const int& n_iter, const float& _dy) {

    float d = sqrt(m_gridWidth * m_gridWidth + m_gridHeight * m_gridHeight);
    float dy = _dy;
    int src_index;
    float n_iterf = (float) n_iter;
    float dy0 = 0.2 * dy;
    float dyn = 0.01 * dy;
    float w = 10.f;
    float v, a, b, c, dist;

    if (g_fault_mode == 10) {
        dy = -dy;
    }

    for (int k = 0; k < n_iterf; k++) {

        if (g_fault_circlemode < 1) {

            v = rand();
            a = cos(v);
            b = sin(v);
            c = ((float)rand() / RAND_MAX) * d - d / 2.0f;

            if (mode == 1) {
                if (k < n_iterf) {
                    dy = dy0 + ((float)k / n_iterf) * (dyn / dy0);
                }
                else {
                    dy = dyn;
                }
            }

            for (int i = 0; i < m_gridHeight; i++) {
                for (int j = 0; j < m_gridWidth; j++) {

                    src_index = j + m_gridWidth * i;
                    dist = a * i + b * j - c;

                    if (mode == 2) {
                        dy = atan(dist) * 0.064f;

                    }
                    else if (mode == 3) {
                        dy = exp(-(pow(dist, 2)) / 10.f) * 0.1f;

                    }
                    else if (mode == 4) {
                        dy = -exp(-(pow(dist, 2)) / 10.f) * 0.1f;

                    }
                    else if (mode == 5) {
                        dy = (0.1f - exp(-(pow(dist, 2)) / 10.f)) * 0.1f;

                    }
                    else if (mode == 6) {
                        if (dist < w) {
                            dy = cos((dist * 3.14f / w - 3.14)) * 0.1f;

                        }
                        else {
                            dy = 0;
                        }
                    }
                    else if (mode == 6) {
                        if (dist < w) {
                            dy = sin((dist * 3.14f / w - 3.14)) * 0.1f;

                        }
                        else {
                            dy = 0;
                        }
                    }

                    if (mode < 2) {
                        if (a * i + b * j > c) {

                            mesh->setLayerThickness(mesh->getLayerThickness(1, i, j) + dy, 1, i, j);
                        }
                        else {
                            mesh->setLayerThickness(mesh->getLayerThickness(1, i, j) - dy, 1, i, j);
                        }

                    }
                    else {
                        mesh->setLayerThickness(mesh->getLayerThickness(1, i, j) + dy, 1, i, j);

                    }
                }
            }
        }
        else {

            dy = 0.3f;

            int randX = rand() % (m_gridWidth + 1);
            int randZ = rand() % (m_gridHeight + 1);
            int randCircSize = rand() % ((m_gridWidth + m_gridHeight) / 10); 
            for (int i = 0; i < m_gridHeight; i++) {
                for (int j = 0; j < m_gridWidth; j++) {
                    float pd = sqrt((randX - j) * (randX - j) + (randZ - i) * (randZ - i)) * 2.0f / randCircSize; 
                    if (fabs(pd) <= 1.0f) {
                        float diff = (dy / 2.0f + cos(pd * 3.14f) * dy / 2.0f);
                        mesh->setLayerThickness(mesh->getLayerThickness(1, i, j) + diff, 1, i, j);
                    }
                }
            }
        }
    }
    //mesh->init();
}

void Mesh::createGridPerlinNoise() {

    randomGradients.clear();

    int newGridWidth = m_gridWidth + 1;
    int newGridHeight = m_gridHeight + 1;

    for (unsigned int i = 0; i < newGridHeight; i++)
    {
        for (unsigned int j = 0; j < newGridWidth; j++)
        {
            float thetaRandom = 2*3.1415*((float)rand() / RAND_MAX); //entre 0 et 2pi
            glm::vec2 randomGradient(cos(thetaRandom), sin(thetaRandom));//norme 1
            randomGradients.push_back(randomGradient);
        }
    }
}

float Mesh::perlinNoise(float x, float y, int resolutionX, int resolutionY) {

    int newGridWidth = m_gridWidth + 1;

    //le début de la grille de notre mesh commence en -5,-5 et non en 0,0
    //donc on va décaler les indices pour que la grille des gradients 
    //commence en même temps 
    glm::vec2 topLeftGrid(-5, -5);
    float deltaWidth = -topLeftGrid.x;
    float deltaHeight = -topLeftGrid.y;

    //pour se caler sur notre grille
    x += deltaWidth;
    y += deltaHeight;

    //Pour avoir un bruit plus lisse, on rapproche les sommets lors du calculs
    x /= resolutionX;
    y /= resolutionY;

    glm::vec2 vertex(x, y); //notre sommet du mesh

    //on regarde les 4 coins de la grille autour du sommet
    int leftJ = (int)vertex.x;
    int rightJ = leftJ + 1;
    int topI = (int)vertex.y;
    int bottomI = topI + 1;

    glm::vec2 topRight(rightJ, topI);
    glm::vec2 bottomRight(rightJ, bottomI);
    glm::vec2 topLeft(leftJ, topI);
    glm::vec2 bottomLeft(leftJ, bottomI);

    //Coefficients d'interpolation
    float dx = x - leftJ;
    float dy = y - topI;

    //std::cout << "x = "<<x<< " ,leftJ = " << leftJ << " ,dx = " << dx << " ,y = "<<y<< " ,topI = " << topI << " ,dy = " << dy << std::endl;
    
    //calcul des produits scalaires
    //d'abord on calcule l'offset avec chaque coin de la cellule
    glm::vec2 offsetTopRight(topRight - vertex);
    //ensuite on prend le gradient du sommet en haut à droite du sommet courant
    glm::vec2 gradientTopRight(randomGradients.at(rightJ + topI * newGridWidth)); 
    //puis on prend le produit scalaire
    float dotTopRight = glm::dot(offsetTopRight, gradientTopRight);

    glm::vec2 offsetBottomRight(bottomRight - vertex);
    glm::vec2 gradientBottomRight(randomGradients.at(rightJ + bottomI * newGridWidth)); 
    float dotBottomRight = glm::dot(offsetBottomRight, gradientBottomRight);

    glm::vec2 offsetTopLeft(topLeft-vertex);
    glm::vec2 gradientTopLeft(randomGradients.at(leftJ + topI * newGridWidth));
    float dotTopLeft = glm::dot(offsetTopLeft, gradientTopLeft);
    
    glm::vec2 offsetBottomLeft(bottomLeft - vertex);
    glm::vec2 gradientBottomLeft(randomGradients.at(leftJ + bottomI * newGridWidth));
    float dotBottomLeft = glm::dot(offsetBottomLeft, gradientBottomLeft);

    //interpolation bilinéaire sur la grille
    float noise = dx * dy * dotBottomRight + dx * (1 - dy) * dotTopRight + dy * (1 - dx) * dotBottomLeft + (1 - dy) * (1 - dx) * dotTopLeft;

    //std::cout << "noise = " << noise << std::endl;
    return noise;
}


float Mesh::fractionalBrownianMotion(float x, float y, int numberOfOctaves) {
    float lacunarity = 2;
    float attenuation = 0.5;
    float noise = 0;

    float amplitude = 1;
    float frequency = 1;

    for (unsigned int k = 0; k < numberOfOctaves; k++)
    {
        noise += amplitude * perlinNoise(x,y,1,1); //on le fait sur le bruit de perlin
        frequency *= lacunarity;
        amplitude *= attenuation;
    }

    return noise;
}

void Mesh::applyNoise(const int& n_iter, int layerID, float maxNoise, int typeNoise, int numberOfOctaves) {

    for (int n = 0; n < n_iter; n++)
    {
        //on créé la carte de bruit de perlin à chaque itération si le bruit est perlin ou fractional brownian motion 
        if (typeNoise != 0) createGridPerlinNoise();

        int indexVertex= 0;
        std::vector<float> newThicknessTab;

        for (int i = 0; i < m_gridHeight; i++)
        {
            for (int j = 0; j < m_gridWidth; j++)
            {
                float noise = 0;

                if (typeNoise == 0) noise = 2 * ((float)rand() / RAND_MAX) - 1; //random : bruit entre -1 et 1;

                float x = m_vertexPositions[indexVertex];
                indexVertex += 2;
                float y = m_vertexPositions[indexVertex]; //y correspond au z de m_vertexPositions
                indexVertex += 1;

                if (typeNoise == 1) noise = perlinNoise(x,y, g_resolutionX, g_resolutionY); //perlin : bruit entre -1 et 1;

                if (typeNoise == 2) noise = fractionalBrownianMotion(x, y, numberOfOctaves); //fractional brownian motion 

                //map noise
                noise *= maxNoise; 
                
                float newThickness = mesh->getLayerThickness(layerID, i, j) + noise;

                //on sauvegarde pour mapper entre l'épaisseur min et max après
                newThicknessTab.push_back(newThickness);
            }
        }

        //on map l'épaisseur entre le min et le max
        float minThickness = *std::min_element(newThicknessTab.begin(), newThicknessTab.end());
        float maxThickness = *std::max_element(newThicknessTab.begin(), newThicknessTab.end());

        for (int i = 0; i < m_gridHeight; i++)
        {
            for (int j = 0; j < m_gridWidth; j++)
            {
                float newThickness = (g_max_height_perlin - g_min_height_perlin)*(newThicknessTab.at(j + i * m_gridWidth)-minThickness)/(maxThickness-minThickness)+ g_min_height_perlin;
                mesh->setLayerThickness(newThickness, layerID, i, j);
            }
        }
    }
    //mesh->init();

}


GLuint loadTextureFromFileToGPU(const std::string &filename) {
  int width, height, numComponents;
  // Loading the image in CPU memory using stbd_image
  unsigned char *data = stbi_load(
    filename.c_str(),
    &width, &height,
    &numComponents, // 1 for a 8 bit greyscale image, 3 for 24bits RGB image, 4 for 32bits RGBA image
    0);
  GLuint texID;
  glGenTextures(1, &texID); // generate an OpenGL texture container
  glBindTexture(GL_TEXTURE_2D, texID); // activate the texture
  // The following lines setup the texture filtering option and repeat mode; check www.opengl.org for details.
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  // fills the GPU texture with the data stored in the CPU image
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

  // Freeing the now useless CPU memory
  stbi_image_free(data);
  glBindTexture(GL_TEXTURE_2D, 0); // unbind the texture

  return texID;
}

void printHelp()
{
    std::cout <<
        "> Help:" << std::endl <<
        "    Mouse commands:" << std::endl <<
        "    * Left button: rotate camera" << std::endl <<
        "    * Middle button: zoom" << std::endl <<
        "    * Right button: pan camera" << std::endl <<
        "    Keyboard commands:" << std::endl <<
        "    * H: print this help" << std::endl <<
        "    * E: Thermal Erosion" << std::endl <<
        "    * ESC: quit the program" << std::endl;
}

// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window.
void windowSizeCallback(GLFWwindow* window, int width, int height) {
  g_camera.setAspectRatio(static_cast<float>(width)/static_cast<float>(height));
  glViewport(0, 0, (GLint)width, (GLint)height); // Dimension of the rendering region in the window
}

// Executed each time a key is entered.
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
  if(action == GLFW_PRESS && key == GLFW_KEY_W) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  } else if(action == GLFW_PRESS && key == GLFW_KEY_F) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  } else if(action == GLFW_PRESS && (key == GLFW_KEY_ESCAPE || key == GLFW_KEY_Q)) {
    glfwSetWindowShouldClose(window, true); // Closes the application if the escape key is pressed
  } else if (action == GLFW_PRESS && key == GLFW_KEY_H) {
     printHelp();
  } else if (action == GLFW_PRESS && key == GLFW_KEY_E) {
      mesh->applyNThermalErosion(100, 0.52, 0.3, 0.001, true, true, true, true, true); //0.52 = 30 degrés  dt = 0.001
      std::cout << "Thermal Erosion" << std::endl;
  }

}


// Called each time the mouse cursor moves
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos)
{
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    const float normalizer = static_cast<float>((width + height) / 2);
    const float dx = static_cast<float>((g_baseX - xpos) / normalizer);
    const float dy = static_cast<float>((ypos - g_baseY) / normalizer);
    if (g_rotatingP) {
        const glm::vec3 dRot(-dy * PI, dx * PI, 0.0);
        g_camera.setRotation(g_baseRot + dRot);
    }
    else if (g_panningP) {
        g_camera.setPosition(g_baseTrans + g_meshScale * glm::vec3(dx, dy, 0.0));
    }
    else if (g_zoomingP) {
        g_camera.setPosition(g_baseTrans + g_meshScale * glm::vec3(0.0, 0.0, dy));
    }
}

// Called each time a mouse button is pressed
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        if (!g_rotatingP) {
            g_rotatingP = true;
            glfwGetCursorPos(window, &g_baseX, &g_baseY);
            g_baseRot = g_camera.getRotation();
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        g_rotatingP = false;
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        if (!g_panningP) {
            g_panningP = true;
            glfwGetCursorPos(window, &g_baseX, &g_baseY);
            g_baseTrans = g_camera.getPosition();
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
        g_panningP = false;
    }
    else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS) {
        if (!g_zoomingP) {
            g_zoomingP = true;
            glfwGetCursorPos(window, &g_baseX, &g_baseY);
            g_baseTrans = g_camera.getPosition();
        }
    }
    else if (button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_RELEASE) {
        g_zoomingP = false;
    }
}

void errorCallback(int error, const char *desc) {
  std::cout <<  "Error " << error << ": " << desc << std::endl;
}

void initGLFW() {
  glfwSetErrorCallback(errorCallback);

  // Initialize GLFW, the library responsible for window management
  if(!glfwInit()) {
    std::cerr << "ERROR: Failed to init GLFW" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Before creating the window, set some option flags
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create the window
  g_window = glfwCreateWindow(
      g_windowWidth, g_windowHeight,
      "Projet IGR 205",
      nullptr, nullptr);

  if(!g_window) {
    std::cerr << "ERROR: Failed to open window" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  // Load the OpenGL context in the GLFW window using GLAD OpenGL wrangler
  glfwMakeContextCurrent(g_window);
  glfwSetWindowSizeCallback(g_window, windowSizeCallback);
  glfwSetKeyCallback(g_window, keyCallback);
  glfwSetCursorPosCallback(g_window, cursorPosCallback);
  glfwSetMouseButtonCallback(g_window, mouseButtonCallback);
  glfwSwapInterval(1); // Enable vsync
}

void initOpenGL() {
  // Load extensions for modern OpenGL
  if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cerr << "ERROR: Failed to initialize OpenGL context" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }
  
  //glCullFace(GL_BACK); // specifies the faces to cull (here the ones pointing away from the camera)
  //glEnable(GL_CULL_FACE); // enables face culling (based on the orientation defined by the cw/ccw enumeration).
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthFunc(GL_LESS);   // Specify the depth test for the z-buffer
  glEnable(GL_DEPTH_TEST);      // Enable the z-buffer test in the rasterization
  glClearColor(0.7f, 0.7f, 0.7f, 1.0f); // specify the background color, used any time the framebuffer is cleared
}

// Loads the content of an ASCII file in a standard C++ string
std::string file2String(const std::string &filename) {
  std::ifstream t(filename.c_str());
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

// Loads and compile a shader, before attaching it to a program
void loadShader(GLuint program, GLenum type, const std::string &shaderFilename) {
  GLuint shader = glCreateShader(type); // Create the shader, e.g., a vertex shader to be applied to every single vertex of a mesh
  std::string shaderSourceString = file2String(shaderFilename); // Loads the shader source from a file to a C++ string
  const GLchar *shaderSource = (const GLchar *)shaderSourceString.c_str(); // Interface the C++ string through a C pointer
  glShaderSource(shader, 1, &shaderSource, NULL); // load the vertex shader code
  glCompileShader(shader);
  glAttachShader(program, shader);
  glDeleteShader(shader);
}

void initGPUprogram() {
  g_program = glCreateProgram(); // Create a GPU program, i.e., two central shaders of the graphics pipeline
  loadShader(g_program, GL_VERTEX_SHADER, "../vertexShader.glsl");
  loadShader(g_program, GL_FRAGMENT_SHADER, "../fragmentShader.glsl");
  glLinkProgram(g_program); // The main GPU program is ready to be handle streams of polygons

  glUseProgram(g_program);
}


void initCamera() {
  int width, height;
  glfwGetWindowSize(g_window, &width, &height);
  g_camera.setAspectRatio(static_cast<float>(width)/static_cast<float>(height));

  g_camera.setPosition(glm::vec3(0.0, 0.0, 10.0));
  g_camera.setNear(0.1);
  g_camera.setFar(20.0);
}

void render();
void clear();

void renderImGui() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::Begin("Control Widget");
    ImGui::SetWindowSize(ImVec2(0, 0));

    ImGui::Separator();

    ImGui::Text("Erosion:");

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::Button("Start thermal erosion")) {
        g_nbOfIterations_t = g_iter_t;
        g_nbOfIterations_t_max = g_iter_t;
    }

    ImGui::ProgressBar(1.f - g_nbOfIterations_t / (float) g_nbOfIterations_t_max, ImVec2(0, 20));

    if (ImGui::TreeNode("Parameters thermal")) {

        ImGui::Spacing();
        ImGui::Text("Erosion's strategy :");
        ImGui::Spacing();

        
        if (ImGui::RadioButton("Strategy A", &g_strategyErosion_t, 0)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Strategy B", &g_strategyErosion_t, 1)) {
        }

        ImGui::Spacing();
        ImGui::Text("Type of erosion :");
        ImGui::Spacing();

        
        if (ImGui::RadioButton("Layer breaks into sand", &g_typeErosion_t, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Layer falls down to same layer", &g_typeErosion_t, 0)) {
        }

        ImGui::Spacing();
        ImGui::Text("Connexity :");
        ImGui::Spacing();

        
        if (ImGui::RadioButton("8", &g_connexity_t, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("4", &g_connexity_t, 0)) {
        }

        ImGui::Spacing();
        ImGui::Text("Descent's direction (strategy A) :");
        ImGui::Spacing();

        
        if (ImGui::RadioButton("Lowest neighbour", &g_descentDirection_t, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Heighest gradient", &g_descentDirection_t, 0)) {
        }

        ImGui::Spacing();
        ImGui::Text("Gradient's type (strategy A) :");
        ImGui::Spacing();

        
        if (ImGui::RadioButton("i+1 | i", &g_typeGradient_t, 2)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("i | i-1", &g_typeGradient_t, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("i+1 | i-1", &g_typeGradient_t, 0)) {
        }

        ImGui::Spacing();
        ImGui::Text("Matter's receiver (strategy A) :");
        ImGui::Spacing();

        
        if (ImGui::RadioButton("With descent direction", &g_neighbourReceiver_t, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("All neighbours", &g_neighbourReceiver_t, 0)) {
        }

        ImGui::SliderFloat("Theta", &g_thetaLimit_t, 0, PI / 2);
        ImGui::SliderFloat("Erosion coefficient", &g_erosionCoeff_t, 0, 1);
        ImGui::SliderFloat("Dt thermal", &g_dt_t, 0.000001f, 0.1f, "%f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderInt("Number of iterations thermal", &g_iter_t, 1, 1000);
        ImGui::TreePop();
    }

    ImGui::Spacing();

    if (ImGui::Button("Start hydraulic erosion")) {
        g_nbOfIterations_h = g_iter_h;
        g_nbOfIterations_h_max = g_iter_h;
    }

    ImGui::ProgressBar(1.f - g_nbOfIterations_h / (float)g_nbOfIterations_h_max, ImVec2(0, 20));

    if (ImGui::Button("Water")) {
        if (g_waterarriving) g_waterarriving = false;
        else g_waterarriving = true;
    }

    if (ImGui::TreeNode("Parameters hydraulic")) {

        ImGui::SliderFloat("Dt hydraulic", &g_dt_h, 0.000001f, 0.1f, "%f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderInt("Number of iterations hydraulic", &g_iter_h, 1, 1000);
        ImGui::SliderFloat3("kc, ks, kd", g_k_h, 0.0f, 1.0f, "%f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderInt2("Water coords", g_water_coords, 0, 100);
        ImGui::SliderFloat("Water flow", &g_water_flow, -100.f, 100.f);
        ImGui::TreePop();
    } 

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    if (ImGui::Button("Fault Algorithm")) {

        if (g_fault_mode == 10) {
            g_fault_nbOfIterations = 1;
            g_fault_nbOfIterations_max = 1;
        }
        else {
            g_fault_nbOfIterations = g_fault_niter;
            g_fault_nbOfIterations_max = g_fault_niter;
        }
    }

    ImGui::ProgressBar(1.f - g_fault_nbOfIterations / (float) g_fault_nbOfIterations_max, ImVec2(0, 20));

    if (ImGui::TreeNode("Parameters fault algorithm")) {

        if (ImGui::RadioButton("sin", &g_fault_mode, 7)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("cos", &g_fault_mode, 6)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("1-exp", &g_fault_mode, 5)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("-exp", &g_fault_mode, 4)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("exp", &g_fault_mode, 3)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("atan", &g_fault_mode, 2)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("dy variant", &g_fault_mode, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("dy constant", &g_fault_mode, 0)) {
        }

        if (ImGui::RadioButton("Line", &g_fault_circlemode, 0)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Circle", &g_fault_circlemode, 1)) {
        }
        ImGui::SliderInt("Number of iterations", &g_fault_niter, 1, 1000);
        ImGui::SliderFloat("Decrement", &g_fault_dy, 0, 2.f);

        ImGui::TreePop();
    }

    if (ImGui::Button("Noise Algorithm")) {

        if (g_fault_mode == 10) {
            g_perl_nbOfIterations = 1;
            g_perl_nbOfIterations_max = 1;
        }
        else {
            g_perl_nbOfIterations = g_perl_niter;
            g_perl_nbOfIterations_max = g_perl_niter;
        }
    }

    ImGui::ProgressBar(1.f - g_perl_nbOfIterations / (float) g_perl_nbOfIterations_max, ImVec2(0, 20));

    if (ImGui::TreeNode("Parameters noise")) {

        ImGui::Spacing();
        ImGui::Text("Type of noise");
        ImGui::Spacing();


        if (ImGui::RadioButton("Random", &g_noise_type, 0)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Perlin", &g_noise_type, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Fractional Brownian Motion", &g_noise_type, 2)) {
        }
        ImGui::Spacing();
        ImGui::Text("Layer to apply noise");
        ImGui::Spacing();


        if (ImGui::RadioButton("Layer 0", &g_perl_layer, 0)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Layer 1", &g_perl_layer, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Layer 2", &g_perl_layer, 2)) {
        }

        ImGui::SliderInt("Number of iterations", &g_perl_niter, 1, 1000);
        ImGui::SliderFloat("Max noise", &g_max_noise, 0, abs(g_h_max-g_h_min)*3); //abs(g_h_max-g_h_min)*3 altitude max des 3 layers
        ImGui::SliderInt("Max height", &g_max_height_perlin, 1, 6);
        ImGui::SliderInt("Min height", &g_min_height_perlin, 0, g_max_height_perlin);
        ImGui::SliderInt("Resolution en X", &g_resolutionX, 1, 10);
        ImGui::SliderInt("Resolution en Y", &g_resolutionY, 1, 10);
        ImGui::SliderInt("Number of octaves (for FBM)", &g_number_octaves, 1, 100);
        

        ImGui::TreePop();
    }

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    if (ImGui::Button("Noise terrain generator")) {
    }

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::Button("Restart")) {
        g_nbOfIterations_t = 0;
        g_nbOfIterations_h = 0;
        g_fault_nbOfIterations = 0;
        g_perl_nbOfIterations = 0;
        //mesh = new Mesh({ "../data/simpleB.png", "../data/simpleS.png" }, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f), glm::vec3(148.f / 255.f, 124.f / 255.f, 48.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(g_h_min, g_h_max)); //cpu
        //mesh->init();
        mesh = new Mesh(mesh->getLayersFileNames(), mesh->getLayersColors(), glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(g_h_min, g_h_max)); //cpu
        //mesh->init();
    }

    if (ImGui::TreeNode("Layers")) {

        static char heightMapName[255] = "../data/simpleB.png";
        ImGui::InputText("input text", heightMapName, IM_ARRAYSIZE(heightMapName));
        static float colorAdd[] = { 255.f / 255.f, 255.f / 255.f, 255.f / 255.f, 255.f / 255.f };
        ImGui::ColorEdit4("Color", colorAdd);

        ImGui::TextDisabled("(?)");

        if (ImGui::IsItemHovered())
        {
            ImGui::BeginTooltip();
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
            ImGui::TextUnformatted("dddd");
            ImGui::PopTextWrapPos();
            ImGui::EndTooltip();
        }

        if (ImGui::RadioButton("Layer 0", &g_layer, 0)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Layer 1", &g_layer, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Layer 2", &g_layer, 2)) {
        }

        if (ImGui::Button("Add")) {

            int width, height, channels;

            unsigned char* heightMap = stbi_load(
                ((std::string) heightMapName).c_str(),
                &width, &height,
                &channels, 
                0);

            if (heightMap == NULL) {
                printf("Error in loading the height map image.\n");

            } else {
                std::vector<std::string> fileNames = mesh->getLayersFileNames();
                fileNames[g_layer] = heightMapName;
                std::vector<glm::vec4> colors = mesh->getLayersColors();
                colors[g_layer] = glm::vec4(colorAdd[0], colorAdd[1], colorAdd[2], colorAdd[3]);
                mesh = new Mesh(fileNames, colors, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(g_h_min, g_h_max));
                //mesh->init();
            }
        }

        ImGui::DragFloatRange2("Range height", &g_h_min, &g_h_max, 0.01f, 0.0f, 10.0f, "Min: %.01f %%", "Max: %.01f %%", ImGuiSliderFlags_AlwaysClamp);

        if (ImGui::Button("+")) {
            mesh->applyFault(10, 1, g_dy);
        }
        ImGui::SameLine();
        if (ImGui::Button("-")) {
            mesh->applyFault(10, 1, -g_dy);
        }

        ImGui::TreePop();
    }

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();
    ImGui::Spacing();


    if (ImGui::TreeNode("Colors")) {

        static float color1[] = { 237.f / 255.f, 224.f / 255.f, 81.f / 255.f };
        static float color0[] = { 120.f / 255.f, 135.f / 255.f, 124.f / 255.f };
        ImGui::ColorEdit3("Color layer 1", color1);
        ImGui::ColorEdit3("Color layer 0", color0);
        if (ImGui::Button("Update color")) {
            glfwSetWindowTitle(g_window, "Un projet incroyable");
            mesh->setLayersColors(1, color1);
            mesh->setLayersColors(0, color0);
            //mesh->init(true);
        }
        ImGui::TreePop();
    }

    ImGui::Separator();

    if (ImGui::TreeNode("Display")) {
        static int e = 0;
        if (ImGui::RadioButton("Edged", &e, 1)) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Full", &e, 0)) {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }
        ImGui::TreePop();
    }

    ImGui::Separator();

    /* (ImGui::TreeNode("Range Widgets"))
    {
        static float begin = 10, end = 90;
        //static int begin_i = 100, end_i = 1000;
        //ImGui::DragFloatRange2("range float", &g_h_min, &g_h_max, 0.01f, 0.0f, 10.0f, "Min: %.01f %%", "Max: %.01f %%", ImGuiSliderFlags_AlwaysClamp);
        //ImGui::DragIntRange2("range int", &begin_i, &end_i, 5, 0, 1000, "Min: %d units", "Max: %d units");
        //ImGui::DragIntRange2("range int (no bounds)", &begin_i, &end_i, 5, 0, 0, "Min: %d units", "Max: %d units");
        ImGui::TreePop();
    }*/

    ImGui::Separator();

    ImGui::Spacing();
    ImGui::Spacing();

    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

    ImGui::End();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void initImGui() {
    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    // Setup Platform/Renderer bindings
    ImGui_ImplGlfw_InitForOpenGL(g_window, true);
    ImGui_ImplOpenGL3_Init("#version 150");
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
}

void init() {
    initGLFW();
    initOpenGL();
    //mesh = new Mesh({ "../data/heightmap5.png" }, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f)}, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 1.f)); //cpu
    //mesh = new Mesh({ "../data/simpleB.png", "../data/simpleS.png", "../data/simpleB.png" }, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f), glm::vec3(148.f / 255.f, 124.f / 255.f, 48.f / 255.f), glm::vec3(0.f / 255.f, 0.f / 255.f, 255.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 5.f)); //cpu
    mesh = new Mesh({ "../data/simpleB.png", "../data/sand-with-water.png","../data/sand-with-water.png" }, { glm::vec4(120.f / 255.f, 135.f / 255.f, 124.f / 255.f, 1.f), glm::vec4(148.f / 255.f, 124.f / 255.f, 48.f / 255.f, 1.f), glm::vec4(39.f / 255.f, 112.f / 255.f, 125.f / 255.f, 100.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 2.f)); //cpu
    //mesh = new Mesh({ "../data/simpleB.png", "../data/simpleStr.png","../data/simpleB.png" }, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f), glm::vec3(148.f / 255.f, 124.f / 255.f, 48.f / 255.f), glm::vec3(0.f / 255.f, 0.f / 255.f, 255.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 5.f)); //cpu
    //mesh = new Mesh({ "../data/simpleB.png", "../data/simpleW.png" }, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f), glm::vec3(20.f / 255.f, 107.f / 255.f, 150.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 1.f)); //cpu
    initGPUprogram();
    //g_sunID = loadTextureFromFileToGPU("../data/heightmap3.jpg");
    //mesh->init(); //gpu
    initCamera();
    initImGui();
}

void clear() {
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  glDeleteProgram(g_program);

  glfwDestroyWindow(g_window);
  glfwTerminate();
}

// The main rendering call
void render() {

    const glm::mat4 viewMatrix = g_camera.computeViewMatrix();
    const glm::mat4 projMatrix = g_camera.computeProjectionMatrix();

    glUniformMatrix4fv(glGetUniformLocation(g_program, "viewMat"), 1, GL_FALSE, glm::value_ptr(viewMatrix)); // compute the view matrix of the camera and pass it to the GPU program
    glUniformMatrix4fv(glGetUniformLocation(g_program, "projMat"), 1, GL_FALSE, glm::value_ptr(projMatrix)); // compute the projection matrix of the camera and pass it to the GPU program
    
    mesh->render();
}

int main(int argc, char ** argv) {
    init(); // Your initialization code (user interface, OpenGL states, scene with geometry, material, lights, etc)

    while(!glfwWindowShouldClose(g_window)) {

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Erase the color and z buffers.

        if (g_nbOfIterations_t > 0) {
            mesh->applyNThermalErosion(1, g_thetaLimit_t, g_erosionCoeff_t, g_dt_t, g_neighbourReceiver_t, g_descentDirection_t, g_typeErosion_t, g_connexity_t, g_strategyErosion_t);
            g_nbOfIterations_t -= 1;
        }

        if (g_fault_nbOfIterations > 0) {
            mesh->applyFault(g_fault_mode, 1, g_fault_dy);            
            g_fault_nbOfIterations -= 1;
        }

        if (g_nbOfIterations_h > 0) {
            mesh->hydraulicErosion(1, g_dt_h, g_k_h); //J'ai mis à 1 pour avoir une iteration par frame c'est plus smooth
            g_nbOfIterations_h -= 1;
        }

        if (g_perl_nbOfIterations > 0) {
            mesh->applyNoise(1, g_perl_layer, g_max_noise, g_noise_type, g_number_octaves);
            g_perl_nbOfIterations -= 1;
        }
        
        mesh->init(false);
        render();

        mesh->init(true);
        render();


        renderImGui();

        glfwSwapBuffers(g_window);
        glfwPollEvents();
    }

    // Cleanup
    clear();
    return EXIT_SUCCESS;
}

