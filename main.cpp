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

// Window parameters
GLFWwindow *g_window = nullptr;
int g_windowWidth = 1024;
int g_windowHeight = 768;

//Simulation parameters
float g_thetaLimit_t = 0.3;
float g_erosionCoeff_t = 0.3;
float g_dt_t = 0.0500;
int g_iter_t = 5;
int g_strategyErosion_t = 0; //stratégie A
int g_typeErosion_t = 1; //Layer breaks into sand
int g_connexity_t = 0; // 4 connexité
int g_descentDirection_t = 0; //Heighest gradient
int g_typeGradient_t = 0; //Heighest gradient
int g_neighbourReceiver_t = 0; //all neighbors
unsigned int g_nbOfIterations_t = 1;
unsigned int g_nbOfIterations_t_max = 1;

//Fault parameters
int g_fault_mode = 2;
int g_fault_niter = 1;
int g_fault_circlemode = 0;
unsigned int g_fault_nbOfIterations = 1;
unsigned int g_fault_nbOfIterations_max = 1;

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
    Mesh(const std::vector<std::string>& filenames, const std::vector<glm::vec3> layersColor, const glm::vec4& corners, const glm::vec2& e);
    Mesh(const int nbOfLayers, const int width, const int height, const int mode, const std::vector<glm::vec3> layersColor, const glm::vec4& corners, const glm::vec2& e);

    void init(); //Génère la surface à partir des différentes épaisseur de niveau, ces normales, puis envoie l'info au GPU
    void render();
    unsigned int getIndex(unsigned int k, unsigned int i, unsigned int j);
    float getH(unsigned int i, unsigned int j) const; //Donne la hauteur du terrain eau comprise
    float getTerrainH(unsigned int i, unsigned int j) const; //Pour avoir la hauteur issue des différentes epaisseurs des layers au point (i,j) de la grille
    float getLayerThickness(unsigned int k, unsigned int i, unsigned int j) const;
    float getLayerH(unsigned int k, unsigned int i, unsigned int j) const;
    void setLayerThickness(float value, unsigned int k, unsigned int i, unsigned int j);
    unsigned int getTopLayerId(unsigned int i, unsigned int j) const; //Id de 0 à m_nbOfLayers - 1 correspond à indice dans m_layersColor
    unsigned int getTerrainTopLayerId(unsigned int i, unsigned int j) const; //Id de 0 à m_nbOfLayers - 2 correspond on exclu l'eau
    glm::uvec2 getLowestNeighbor(unsigned int i, unsigned int j) const;
    float getIDerivate(unsigned int i, unsigned int j) const; //Dervive par rapport à i c'est à dire quand passe de ligne i à ligne i + 1 (Z)
    float getJDerivate(unsigned int i, unsigned int j) const; //Derive par rapport à j (X)
    glm::vec2 getGradient(unsigned int i, unsigned int j) const;
    void setLayersColors(int layer, float color[]);
    std::vector<glm::uvec2> getAllLowNeighbors(unsigned int i, unsigned int j, bool connexity8) const;
    void thermalErosionA(float thetaLimit,float erosionCoeff,float dt, bool neighbourReceiver, bool descentDirection, bool typeErosion, bool connexity8);
    void thermalErosionB(float thetaLimit, float erosionCoeff, float dt, bool connexity8);
    void applyNThermalErosion(unsigned int N, float thetaLimit, float erosionCoeff, float dt, bool neighbourReceiver, bool descentDirection, bool typeErosion, bool connexity8, bool strategyB);
    void applyFault(const int& mode, const int& n_iter);


private:
    unsigned int m_gridWidth = 0; //Nb de colonnes de la grille de discretisation (image)
    unsigned int m_gridHeight = 0; //Nb de lignes de la grille de discretisation (image)
    unsigned int m_nbOfLayers = 0;
    glm::vec2 m_gridTopLeftCorner;
    glm::vec2 m_gridBottomRightCorner;
    float m_cellWidth = 0;
    float m_cellHeight = 0;
    std::vector<float> m_layersThickness;
    std::vector<glm::vec3> m_layersColor;

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
    unsigned char* Mesh::createHeightMapFault(const int& width, const int& height, const int& mode, const int& n_iter);
    unsigned char* loadHeightMapFromFile(const std::string& filename, int& width, int& height, int& channels);
};

//Mesh CPU
Mesh* mesh;

void Mesh::init() {
    //Update cpu
    float ax = m_gridTopLeftCorner.x; //a coin en haut gauche
    float az = m_gridTopLeftCorner.y;
    float bx = m_gridBottomRightCorner.x; //b coin en bas droite
    float bz = m_gridBottomRightCorner.y;
    glm::vec2 grad(0.f, 0.f);
    glm::vec3 normal(0.f, 0.f, 0.f);
    glm::vec3 color(0.f, 0.f, 0.f);
    unsigned int ind = 0;
    for (int i = 0; i < m_gridHeight; i++) {
        for (int j = 0; j < m_gridWidth; j++) {//On calcule les vecteurs normaux, en utilisant le gradient de la fonction d'élévation
            grad = getGradient(i, j);
            normal = glm::normalize(glm::vec3(-grad.y, 1, -grad.x)); //On fait attention bien mettre dans bon ordre i.e. derive par rapport à x, correspond à derive par rapport à j...
            color = m_layersColor[getTopLayerId(i, j)];
            //glm::vec2 gradN = (grad.x >0 && grad.y > 0) ? glm::normalize(grad) : grad;
            //color.x = 1 / 2.f * (gradN.x + 1);
            //color.y = 1 / 2.f * (gradN.y + 1);
            //color.z = 1 / 2.f;
            //color = glm::vec3(0.5, 0.5, 0.5);
            //color = m_layersColor[getTopLayerId(i, j)];

            m_vertexPositions[ind] = (ax + (bx - ax) * j / (m_gridWidth - 1)); //x
            m_vertexNormals[ind] = normal.x;
            m_vertexColors[ind] = color.x;
            ind += 1;

            m_vertexPositions[ind] = getH(i, j); //y
            m_vertexNormals[ind] = normal.y;
            m_vertexColors[ind] = color.y;
            ind += 1;

            m_vertexPositions[ind] = (az + (bz - az) * i / (m_gridWidth - 1)); //z
            m_vertexNormals[ind] = normal.z;
            m_vertexColors[ind] = color.z;
            ind += 1;
        }
    }

    //Init gpu
    // Create a single handle that joins together attributes (vertex positions,
    // normals) and connectivity (triangles indices)
    glCreateVertexArrays(1, &m_vao);
    glBindVertexArray(m_vao);

    // Generate a GPU buffer to store the positions of the vertices
    size_t vertexBufferSize = sizeof(float) * m_vertexPositions.size(); // Gather the size of the buffer from the CPU-side vector
    glCreateBuffers(1, &m_posVbo);
    glNamedBufferStorage(m_posVbo, vertexBufferSize, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    glNamedBufferSubData(m_posVbo, 0, vertexBufferSize, m_vertexPositions.data()); // Fill the data storage from a CPU array
    glBindBuffer(GL_ARRAY_BUFFER, m_posVbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

    //Buffer for normal
    size_t vertexBufferSize2 = sizeof(float) * m_vertexNormals.size(); // Gather the size of the buffer from the CPU-side vector
    glCreateBuffers(1, &m_normalVbo);
    glNamedBufferStorage(m_normalVbo, vertexBufferSize2, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    glNamedBufferSubData(m_normalVbo, 0, vertexBufferSize2, m_vertexNormals.data()); // Fill the data storage from a CPU array
    glBindBuffer(GL_ARRAY_BUFFER, m_normalVbo);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

    //Buffer for color
    size_t vertexBufferSize3 = sizeof(float) * m_vertexColors.size(); // Gather the size of the buffer from the CPU-side vector
    glCreateBuffers(1, &m_colVbo);
    glNamedBufferStorage(m_colVbo, vertexBufferSize3, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    glNamedBufferSubData(m_colVbo, 0, vertexBufferSize3, m_vertexColors.data()); // Fill the data storage from a CPU array
    glBindBuffer(GL_ARRAY_BUFFER, m_colVbo);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), 0);

    //Buffer for texture
    //size_t vertexBufferSize3 = sizeof(float) * m_vertexTexCoords.size(); // Gather the size of the buffer from the CPU-side vector
    //glCreateBuffers(1, &m_texCoordsVbo);
    //glNamedBufferStorage(m_texCoordsVbo, vertexBufferSize3, NULL, GL_DYNAMIC_STORAGE_BIT); // Create a data storage on the GPU
    //glNamedBufferSubData(m_texCoordsVbo, 0, vertexBufferSize3, m_vertexTexCoords.data()); // Fill the data storage from a CPU array
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

    glBindVertexArray(0); // deactivate the VAO for now, will be activated at rendering time

    //unsigned int fbo;
    //glGenFramebuffers(1, &fbo);
    //glBindFramebuffer(GL_FRAMEBUFFER, fbo);

}

void Mesh::setLayersColors(int layer, float color[]) {
    std::cout << color[0] << " " << color[1] << std::endl;
    m_layersColor[layer] = glm::vec3(color[0], color[1], color[2]);

}


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

unsigned char* Mesh::createHeightMapFault(const int& width, const int& height, const int& mode, const int& n_iter) {

    size_t map_size = width * height;
    unsigned char * heightMap = (unsigned char*) malloc(map_size);

    for (unsigned char * p = heightMap, *pg = heightMap; p != heightMap + map_size; p += 1, pg += 1) {
        *pg = (unsigned char) 115;
    }

    float d = sqrt(width * width + height * height);
    float dy;
    int src_index;
    float n_iterf = (float) n_iter;
    float dy0 = 0.1f;
    float dyn = 1.f;
    float v, a, b, c;
    float dist;
    float w = 10.f;


    for (int k = 0; k < n_iterf; k++) {

        v = rand();
        a = cos(v);
        b = sin(v);

        c = (cos(rand()) / 2.f + 0.5f) * d - d / 2;
        //std::cout << v << " c " << c << " d " << d << " " << a << " " << b << std::endl;

        switch (mode) {
            case 0:
                dy = 0.f;
                break;

            case 1:
                if (k < n_iterf) {
                    dy = dy0 + (k / n_iterf) * (dyn / dy0);
                }
                else {
                    dy = dyn;
                }
                break;

            default:
                dy = 1.f;
        }

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {

                src_index = j + width * i;
                dist = (float) (a * i + b * j - c);

                if (mode == 2) {
                    dy = atan(dist) * 10.f;

                } else if (mode == 3) {
                    dy = exp(- (pow(dist, 2)) / 10.f) * 10.f;

                } else if (mode == 4) {
                    dy = - exp(- (dist) / 10.f) * 10.f;

                } else if (mode == 5) {
                    if (dist < w) {
                        //dy = cos((dist * 3.14f / w - 3.14f)) * 10.f;
                        dy = cos((dist * 3.14f / w - 3.14)) * 10.f;

                    } else {
                        dy = 0;
                    }
                }

                if (a * i + b * j > c) {

                    mesh->setLayerThickness(mesh->getLayerThickness(1, i, j) + dy, 1, i, j);

                    //if ((float) heightMap[src_index] < 256.f - dy) {
                        //heightMap[src_index] += (unsigned char) dy;
                    //}
                }

                else {
                    //if ((float) heightMap[src_index] > dy) {
                      //  heightMap[src_index] -= (unsigned char) dy;
                    //}
                    mesh->setLayerThickness(mesh->getLayerThickness(1, i, j) - dy, 1, i, j);
                }
                    //std::cout << dy << std::endl;
            }
            //std::cout << i << " " << 0 << " " << (float) heightMap[src_index] << " " << (a * i > c) << std::endl;
        }
    }

    //std::cout << " " << dy << " " << dist;

    return heightMap;
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


Mesh::Mesh(const std::vector<std::string>& filenames, const std::vector<glm::vec3> layersColor, const glm::vec4& corners, const glm::vec2& e) {
    int width = 0, height = 0, channels = 0;
    m_nbOfLayers = filenames.size();
    m_layersColor = layersColor;
    m_gridTopLeftCorner = glm::vec2(corners.x, corners.y);
    m_gridBottomRightCorner = glm::vec2(corners.z, corners.w);
    float emin = e.x;
    float emax = e.y;

    for (unsigned int k = 0; k < m_nbOfLayers; k = k + 1) {
        unsigned char* gray_img = loadHeightMapFromFile(filenames[k], width, height, channels);

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
            m_vertexTexCoords.push_back(float(j) / (width - 1)); //x
            m_vertexTexCoords.push_back(float(i) / (height - 1)); //y
        }
    }
}

Mesh::Mesh(const int nbOfLayers, const int width, const int height, const int mode, const std::vector<glm::vec3> layersColor, const glm::vec4& corners, const glm::vec2& e) {
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
}

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

float Mesh::getLayerH(unsigned int k, unsigned int i, unsigned int j) const {
    float h = 0;
    for (unsigned int l = 0; l <= k; l = l + 1) {
        h += getLayerThickness(l, i, j);
    }
    return h;
}

glm::uvec2 Mesh::getLowestNeighbor(unsigned int i, unsigned int j) const {
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

std::vector<glm::uvec2> Mesh::getAllLowNeighbors(unsigned int i, unsigned int j, bool connexity8) const {
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

float Mesh::getLayerThickness(unsigned int k, unsigned int i, unsigned int j) const {
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

float Mesh::getH(unsigned int i, unsigned int j) const {
    //Donne la hauteur en un point avec l'eau comprise (pour hydraulic simulation)
    float h = 0;
    for (unsigned int k = 0; k < m_nbOfLayers; k = k + 1) {
        h += getLayerThickness(k, i, j);
    }
    return h;
}

float Mesh::getTerrainH(unsigned int i, unsigned int j) const {
    //Donne la hauteur en un point du terrain, dont exclu l'eau qui est par convention dans l'index m_nbOfLayers - 1 (pour thermal erosion)
    float h = 0;
    for (unsigned int k = 0; k < m_nbOfLayers - 1; k = k + 1) {
        h += getLayerThickness(k, i, j);
    }
    return h;
}

float Mesh::getIDerivate(unsigned int i, unsigned int j) const {
    //On fait attention au bord
    if (i == 0) {
        return (getTerrainH(i + 1, j) - getTerrainH(i, j)) / m_cellHeight;
    }
    if (i == m_gridHeight - 1) {
        return (getTerrainH(i, j) - getTerrainH(i - 1, j)) / m_cellHeight;
    }
    
    //return (getTerrainH(i, j) - getTerrainH(i - 1, j)) / m_cellHeight;
    return (getTerrainH(i + 1, j) - getTerrainH(i - 1, j)) / (2.f * m_cellHeight);
}

float Mesh::getJDerivate(unsigned int i, unsigned int j) const {
    //On fait attention au bord
    if (j == 0) {
        return (getTerrainH(i, j + 1) - getTerrainH(i, j)) / m_cellWidth;
    }
    if (j == m_gridWidth - 1) {
        return (getTerrainH(i, j) - getTerrainH(i, j - 1)) / m_cellWidth;
    }

    //return (getTerrainH(i, j) - getTerrainH(i, j - 1)) / m_cellWidth;
    return (getTerrainH(i, j + 1) - getTerrainH(i, j - 1)) / (2.f * m_cellWidth);
}

glm::vec2 Mesh::getGradient(unsigned int i, unsigned int j) const {
    return glm::vec2(getIDerivate(i, j), getJDerivate(i, j));
}


void Mesh::thermalErosionA(float thetaLimit, float erosionCoeff, float dt, bool neighbourReceiver, bool descentDirection, bool typeErosion, bool connexity8) {
    float tangentLimit = glm::tan(thetaLimit);
    std::vector<float> newLayersThickness = m_layersThickness; //vecteur mémoire temporaire
    //neigbourReceiver = 0 : tous les voisins plus bas recoivent dh, 1 : voisin dans la direction de descente
    //descentDirection = 0 : gradient direction, 1 : lowest neighbour
    //typeErosion = 1 : L'érosion donne du sable, 0 : l'érosion donne le même layer (et passe sous le sable)


    for (unsigned int i = 0; i < m_gridHeight; i++) {
        for (unsigned int j = 0; j < m_gridWidth; j++)
        {
            glm::vec2 directionDescent = -getGradient(i, j);
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

                float currentThickness = newLayersThickness[layerIndexCurrentCell * m_gridHeight * m_gridWidth + i * m_gridWidth + j];

                //Si dh est plus grand que l'epaisseur du layer
                if (-dh > currentThickness) {
                    dh = -currentThickness;
                }

                //on erode si on est pas le layer le plus bas ou de l'eau
                if (layerIndexCurrentCell > 0) {
                    
                    
                    float newThicknessCurrentCell = currentThickness + dh;

                    newLayersThickness[layerIndexCurrentCell * m_gridHeight * m_gridWidth + i * m_gridWidth + j] = (newThicknessCurrentCell > 0) ? newThicknessCurrentCell : 0.f;
                    
                    //On donne la matière à un voisin (dans la direction de descente)
                    if (neighbourReceiver!=0) {
                        //std::cout << "Un seul voisin reçoit la matière" << std::endl;

                        unsigned int nextI;
                        unsigned int nextJ;
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
                        unsigned int nextI;
                        unsigned int nextJ;
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

    mesh->init();

    
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

    for (unsigned int i = 0; i < m_gridHeight; i++) {
        for (unsigned int j = 0; j < m_gridWidth; j++)
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
    mesh->init();
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

void Mesh::applyFault(const int& mode, const int& n_iter) {

    /*size_t map_size = width * height;
    unsigned char* heightMap = (unsigned char*)malloc(map_size);

    for (unsigned char* p = heightMap, *pg = heightMap; p != heightMap + map_size; p += 1, pg += 1) {
        *pg = (unsigned char)115;
    }*/

    float d = sqrt(m_gridWidth * m_gridWidth + m_gridHeight * m_gridHeight);
    float dy;
    int src_index;
    float n_iterf = (float) n_iter;
    float dy0 = 0.2f;
    float dyn = 0.01f;
    float v, a, b, c;
    float dist;
    float w = 10.f;


    for (int k = 0; k < n_iterf; k++) {

        if (g_fault_circlemode < 1) {

            v = rand();
            a = cos(v);
            b = sin(v);

            //c = (cos(rand()) / 2.f + 0.5f) * d - d / 2.f;
            c = ((float)rand() / RAND_MAX) * d - d / 2.0f;

            switch (mode) {
            case 0:
                dy = 0.05f;
                break;

            case 1:
                if (k < n_iterf) {
                    dy = dy0 + ((float)k / n_iterf) * (dyn / dy0);
                }
                else {
                    dy = dyn;
                }
                break;

            default:
                dy = 0.1f;
            }

            for (int i = 0; i < m_gridHeight; i++) {
                for (int j = 0; j < m_gridWidth; j++) {

                    src_index = j + m_gridWidth * i;
                    //dist = (float) (a * (float) i + b * (float) j - c) / d;
                    dist = a * i + b * j - c;

                    if (mode == 2) {
                        dy = atan(dist) * 0.064f;
                        //dy = dist / (float) 10.f;

                        if (false) {
                            std::cout << dy << " " << dist << std::endl;
                        }

                    }
                    else if (mode == 3) {
                        dy = exp(-(pow(dist, 2)) / 10.f) * 0.1f;

                    }
                    else if (mode == 4) {
                        dy = -exp(-(pow(dist, 2)) / 10.f) * 0.1f;

                    }
                    else if (mode == 5) {
                        dy = (1.f - exp(-(dist) / 10.f)) * 0.1f;

                    }
                    else if (mode == 6) {
                        if (dist < w) {
                            dy = cos((dist * 3.14f / w - 3.14)) * 0.1f;

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

                    //std::cout << a * i + b * j - c << std::endl;


                    //std::cout << dy << std::endl;
                }
                //std::cout << i << " " << 0 << " " << (float) heightMap[src_index] << " " << (a * i > c) << std::endl;

            }
        }
        else {

            dy = 0.3f;

            int randX = rand() % (m_gridWidth + 1);
            int randZ = rand() % (m_gridHeight + 1);
            int randCircSize = rand() % ((m_gridWidth + m_gridHeight) / 10); // circle diameter
            for (int i = 0; i < m_gridHeight; i++) {
                for (int j = 0; j < m_gridWidth; j++) {
                    float pd = sqrt((randX - j) * (randX - j) + (randZ - i) * (randZ - i)) * 2.0f / randCircSize; // pd = distanceFromCircle*2/size
                    if (fabs(pd) <= 1.0f) { // if the vertex is within the circle, displace it upwards
                        float diff = (dy / 2.0f + cos(pd * 3.14f) * dy / 2.0f);
                        mesh->setLayerThickness(mesh->getLayerThickness(1, i, j) + diff, 1, i, j);

                        //terrain[z][x] += (terrain[z][x] + diff > maxHeight) ? 0 : diff; // constrain to maxHeight
                    }
                }
            }
        }
        //std::cout << c << " " << a << std::endl;
    }

    //std::cout << " " << mesh->getLayerThickness(1, 50, 50) << " " << dist;
    mesh->init();
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

    if (ImGui::TreeNode("Parameters")) {

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
        ImGui::SliderFloat("Dt", &g_dt_t, 0.000001f, 0.1f, "%f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderInt("Number of iterations", &g_iter_t, 1, 1000);
        ImGui::TreePop();
    }

    ImGui::Spacing();
    /*
    static float thetaLimit_h = 0.3;
    static float erosionCoeff_h = 0.3;
    static float dt_h = 0.0500;
    static int iter_h = 5;

    if (ImGui::Button("Start hydraulic erosion")) {
        mesh->applyNThermalErosion(iter_h, thetaLimit_h, erosionCoeff_h, dt_h, true);
    }

    if (ImGui::TreeNode("Parameters")) {

        //ImGui::SliderFloat("Theta", &thetaLimit_h, 0, PI / 2);
        //ImGui::SliderFloat("Erosion coefficient", &erosionCoeff_h, 0, 1);
        //ImGui::SliderFloat("Dt", &dt_h, 0.000001f, 0.1f, "%f", ImGuiSliderFlags_Logarithmic);
        //ImGui::SliderInt("Number of iterations", &iter_h, 1, 1000);
        ImGui::TreePop();
        
    }*/

    

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    if (ImGui::Button("Fault Algorithm")) {
        //mesh = new Mesh(2, 100, 100, g_fault_mode, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f), glm::vec3(237.f / 255.f, 224.f / 255.f, 81.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 1.f));
        //mesh->init();

        //mesh->applyFault(g_fault_mode, g_fault_niter);
        //mesh->init();
        g_fault_nbOfIterations = g_fault_niter;
        g_fault_nbOfIterations_max = g_fault_niter;
    }

    if (ImGui::TreeNode("Parameters fault algorithm")) {

        if (ImGui::RadioButton("Mode5", &g_fault_mode, 5)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Mode4", &g_fault_mode, 4)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Mode3", &g_fault_mode, 3)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Mode2", &g_fault_mode, 2)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Mode1", &g_fault_mode, 1)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Mode0", &g_fault_mode, 0)) {
        }

        if (ImGui::RadioButton("Line", &g_fault_circlemode, 0)) {
        }
        ImGui::SameLine();

        if (ImGui::RadioButton("Circle", &g_fault_circlemode, 1)) {
        }
        ImGui::ProgressBar(1.f - g_fault_nbOfIterations / (float) g_fault_nbOfIterations_max, ImVec2(0, 20));
        ImGui::SliderInt("Number of iterations", &g_fault_niter, 1, 1000);

        ImGui::TreePop();
    }

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::Button("Restart")) {
        g_nbOfIterations_t = 1;
        g_fault_nbOfIterations = 1;
        mesh = new Mesh({ "../data/simpleB.png", "../data/simpleS.png", "../data/simpleB.png" }, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f), glm::vec3(148.f / 255.f, 124.f / 255.f, 48.f / 255.f), glm::vec3(0.f / 255.f, 0.f / 255.f, 255.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 5.f)); //cpu
        mesh->init();
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
            mesh->init();
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

    /*if (ImGui::TreeNode("Range Widgets"))
    {
        static float begin = 10, end = 90;
        static int begin_i = 100, end_i = 1000;
        ImGui::DragFloatRange2("range float", &begin, &end, 0.25f, 0.0f, 100.0f, "Min: %.1f %%", "Max: %.1f %%", ImGuiSliderFlags_AlwaysClamp);
        ImGui::DragIntRange2("range int", &begin_i, &end_i, 5, 0, 1000, "Min: %d units", "Max: %d units");
        ImGui::DragIntRange2("range int (no bounds)", &begin_i, &end_i, 5, 0, 0, "Min: %d units", "Max: %d units");
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
    mesh = new Mesh({ "../data/simpleB.png", "../data/sand-with-water.png","../data/water-around-sand.png" }, { glm::vec3(120.f / 255.f, 135.f / 255.f, 124.f / 255.f), glm::vec3(148.f / 255.f, 124.f / 255.f, 48.f / 255.f), glm::vec3(0.f / 255.f, 0.f / 255.f, 255.f / 255.f) }, glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 2.f)); //cpu    initGPUprogram();
    //g_sunID = loadTextureFromFileToGPU("../data/heightmap3.jpg");
    mesh->init(); //gpu
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
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Erase the color and z buffers.

    const glm::mat4 viewMatrix = g_camera.computeViewMatrix();
    const glm::mat4 projMatrix = g_camera.computeProjectionMatrix();

    glUniformMatrix4fv(glGetUniformLocation(g_program, "viewMat"), 1, GL_FALSE, glm::value_ptr(viewMatrix)); // compute the view matrix of the camera and pass it to the GPU program
    glUniformMatrix4fv(glGetUniformLocation(g_program, "projMat"), 1, GL_FALSE, glm::value_ptr(projMatrix)); // compute the projection matrix of the camera and pass it to the GPU program
    
    mesh->render();
}

int main(int argc, char ** argv) {
    init(); // Your initialization code (user interface, OpenGL states, scene with geometry, material, lights, etc)

    while(!glfwWindowShouldClose(g_window)) {
        if (g_nbOfIterations_t > 0) {
            mesh->applyNThermalErosion(1, g_thetaLimit_t, g_erosionCoeff_t, g_dt_t, g_neighbourReceiver_t, g_descentDirection_t, g_typeErosion_t, g_connexity_t, g_strategyErosion_t);
            g_nbOfIterations_t -= 1;
        }

        if (g_fault_nbOfIterations > 0) {
            mesh->applyFault(g_fault_mode, 1);            
            g_fault_nbOfIterations -= 1;
        }

        render();
        renderImGui();

        glfwSwapBuffers(g_window);
        glfwPollEvents();
    }

    // Cleanup
    clear();
    return EXIT_SUCCESS;
}

