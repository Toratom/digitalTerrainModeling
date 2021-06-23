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


//extern "C" {
//    _declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
//}

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
#include <string>

#define _USE_MATH_DEFINES
#include <cmath>
#include <memory>

#include <Windows.h> 

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define PI 3.141592

#define PATCH_HEIGHT 32
#define PATCH_WIDTH 32
#define NB_OF_LAYERS 3
unsigned int g_nbGroupsX = 0;
unsigned int g_nbGroupsY = 0;

// Window parameters
GLFWwindow *g_window = nullptr;
int g_windowWidth = 1024;
int g_windowHeight = 768;

GLFWwindow* g_window2 = nullptr;

//Simulation parameters
float g_dt = 0.005; //0.0005


//GPU objects - Program
GLuint g_program = 0; // A GPU program contains at least a vertex shader and a fragment shader
GLuint g_computeProgram = 0; //GPU program for compute shader
GLuint g_computeProgramHydraulicA = 0;
GLuint g_computeProgramHydraulicB = 0;
GLuint g_computeForRendering = 0;
//GPU objects - Buffers
//GLuint g_colVbo = 0;
GLuint g_gridPosVbo = 0; //Un vbo qui donne le z (i) et x (j) des points de la grille
GLuint g_gridNormalsVbo = 0; //Normale au sol calculer par un compute shader
GLuint g_gridHeightsVbo = 0; //Un vbo qui donne la hauteur du terrain en (i,j), car le vertex shader ne peut pas la calculer à partir des layers
GLuint g_gridColorsVbo = 0;
GLuint g_gridLayersHeightVboR = 0; //b0(i,j), b1(i,j), ..., bNbOfLayers(i,j) ... Donne l'epaisseur des différents layers AU MAX 4 LAYERS (sans compter l'eau)
GLuint g_gridLayersHeightVboW = 0; //R aura le role de t et W de t + dt (puis swap) : en pratique on lit dans celui pointé par R et on écrit dans celui pointé par W
GLuint g_waterFlowsVboR = 0; //B (i - 1, j) , L (i, j - 1), R (i, j + 1), T (i + 1, j)
GLuint g_waterFlowsVboW = 0;
GLuint g_waterVelocity = 0; //En R et W : v (vitesse direction i), u (vitesse direction j)
GLuint g_waterSedimentR = 0;
GLuint g_waterSedimentW = 0;
//... Cf article pour autre buffer necessaire
GLuint g_ibo = 0;
//Variable pour l'affichage
bool g_displayWater = true;
bool g_displayVel = true;
bool g_displaySed = true;

//Debug
float* points;

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
    Mesh(const std::vector<std::string>& filenames, const std::vector<glm::vec4> layersColor, 
        const std::vector<float> layersErosionCoeffs, const std::vector<float> layersThetasLimit,
        const glm::vec4& corners, const glm::vec2& e);
    const std::vector<unsigned int>& getTriangleIndices() const;
    float getLayerThickness(unsigned int k, unsigned int i, unsigned int j) const;
    float getH(unsigned int i, unsigned int j) const;
    glm::vec2 getGridPos(unsigned int i, unsigned int j) const;
    unsigned int getGridWidth() const;
    unsigned int getGridHeight() const;
    unsigned int getNbOfLayers() const;
    float getCellWidth() const;
    float getCellHeight() const;
    std::vector<float> getLayersColor() const; //Attention renvoie un tableau de float, il faut le lire par groupe de 3 pour avoir les couleurs
    std::vector<float> getLayersErosionCoeffs() const;
    std::vector<float> getLayersThetaLimit() const;
    void render();

private:
    unsigned int m_gridWidth = 0; //Nb de colonnes de la grille de discretisation (image)
    unsigned int m_gridHeight = 0; //Nb de lignes de la grille de discretisation (image)
    unsigned int m_nbOfLayers = 0; //AU PLUS 4
    float m_cellWidth = 0;
    float m_cellHeight = 0;
    glm::vec2 m_gridTopLeftCorner;
    glm::vec2 m_gridBottomRightCorner;
    std::vector<float> m_layersThickness;
    std::vector<glm::vec4> m_layersColor;
    std::vector<float> m_layersErosionCoeffs;
    std::vector<float> m_layersThetaLimits;
    std::vector<glm::vec2> m_gridPositions;
    std::vector<unsigned int> m_triangleIndices;
    unsigned int m_vao;

    unsigned char* Mesh::createHeightMapFault(const int& width, const int& height);
    unsigned char* loadHeightMapFromFile(const std::string& filename, int& width, int& height, int& channels);
    void setLayerThickness(float value, unsigned int k, unsigned int i, unsigned int j);
};

//Mesh CPU
Mesh* mesh;

std::vector<float> Mesh::getLayersColor() const {
    std::vector<float> tmp;
    tmp.resize(4 * NB_OF_LAYERS);
    unsigned int k2 = 0;
    glm::vec4 currentColor;
    for (unsigned int k1 = 0; k1 < NB_OF_LAYERS; k1++) {
        currentColor = m_layersColor[k1];
        k2 = 4 * k1;
        tmp[k2] = currentColor.x;
        tmp[k2 + 1] = currentColor.y;
        tmp[k2 + 2] = currentColor.z;
        tmp[k2 + 3] = currentColor.w;
    }
    return tmp;
}

std::vector<float> Mesh::getLayersErosionCoeffs() const {
    return m_layersErosionCoeffs;
}
std::vector<float> Mesh::getLayersThetaLimit() const {
    return m_layersThetaLimits;
}

const std::vector<unsigned int>& Mesh::getTriangleIndices() const {
    return m_triangleIndices;
}

float Mesh::getCellHeight() const {
    return m_cellHeight;
}

float Mesh::getCellWidth() const {
    return m_cellWidth;
}

unsigned int Mesh::getGridHeight() const {
    return m_gridHeight;
}

unsigned int Mesh::getGridWidth() const {
    return m_gridWidth;
}

unsigned int Mesh::getNbOfLayers() const {
    return m_nbOfLayers;
}

float Mesh::getLayerThickness(unsigned int k, unsigned int i, unsigned int j) const {
    if (i < 0 || i >= m_gridHeight || j < 0 || j >= m_gridWidth) {
        std::cout << "WARNING : You are loooking outside the grid" << std::endl;
        return 0;
    }
    return m_layersThickness[k * m_gridHeight * m_gridWidth + i * m_gridWidth + j];
}

float Mesh::getH(unsigned int i, unsigned int j) const {
    float h = 0;
    for (unsigned int k = 0; k < m_nbOfLayers; k = k + 1) {
        h += getLayerThickness(k, i, j);
    }
    return h;
}

void Mesh::setLayerThickness(float value, unsigned int k, unsigned int i, unsigned int j) {
    if (value < 0) {
        value = 0.f;
    }

    m_layersThickness[k * m_gridHeight * m_gridWidth + i * m_gridWidth + j] = value;
}

glm::vec2 Mesh::getGridPos(unsigned int i, unsigned int j) const {
    //Return posI, posJ donc Z puis X
    return m_gridPositions[i * m_gridWidth + j];
}

void Mesh::render() {
    glm::mat4 model = glm::mat4(1.f); //matrice modele

    //On donne la matrice du modèle au vertex shader et les couleurs au fragment shader
    glUniformMatrix4fv(glGetUniformLocation(g_program, "modelMat"), 1, GL_FALSE, glm::value_ptr(model)); 

    glBindVertexArray(m_vao);
    
    // bind the buffers
    glBindBuffer(GL_ARRAY_BUFFER, g_gridPosVbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (void*) 0);

    glBindBuffer(GL_ARRAY_BUFFER, g_gridHeightsVbo);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 1 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, g_gridNormalsVbo);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0); //Attention ici normal en vec4 cf explication dans verxterShader

    glBindBuffer(GL_ARRAY_BUFFER, g_gridColorsVbo);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ibo);


    //Draw
    glDrawElements(GL_TRIANGLES, m_triangleIndices.size(), GL_UNSIGNED_INT, 0); // Call for rendering: stream the current GPU geometry through the current GPU program

    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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

Mesh::Mesh(const std::vector<std::string>& filenames, const std::vector<glm::vec4> layersColor,
    const std::vector<float> layersErosionCoeffs, const std::vector<float> layersThetasLimit,
    const glm::vec4& corners, const glm::vec2& e) {
    int width = 0, height = 0, channels = 0;
    m_nbOfLayers = filenames.size();
    if (m_nbOfLayers != NB_OF_LAYERS) {
        std::cout << "ERROR : wrong nb of layer should correpond to NB_OF_LAYERS = " << NB_OF_LAYERS  << std::endl;
        exit(1);
    }
    m_layersColor = layersColor;
    m_layersErosionCoeffs = layersErosionCoeffs;
    m_layersThetaLimits = layersThetasLimit;
    m_gridTopLeftCorner = glm::vec2(corners.x, corners.y);
    m_gridBottomRightCorner = glm::vec2(corners.z, corners.w);
    float ax = m_gridTopLeftCorner.x; //a coin en haut gauche
    float az = m_gridTopLeftCorner.y;
    float bx = m_gridBottomRightCorner.x; //b coin en bas droite
    float bz = m_gridBottomRightCorner.y;
    float emin = e.x;
    float emax = e.y;

    glGenVertexArrays(1, &m_vao);

    for (unsigned int k = 0; k < m_nbOfLayers; k = k + 1) {
        unsigned char* gray_img = loadHeightMapFromFile(filenames[k], width, height, channels);

        //Met a jour la taille de la grille avec la valeur de la première map
        if (k == 0) {
            m_gridWidth = width;
            m_gridHeight = height;

            m_cellWidth = abs(m_gridBottomRightCorner.x - m_gridTopLeftCorner.x) / (width - 1.f);
            m_cellHeight = abs(m_gridBottomRightCorner.y - m_gridTopLeftCorner.y) / (height - 1.f);

            std::cout << m_cellHeight << " " << m_cellWidth << std::endl;

            m_layersThickness.resize(m_nbOfLayers * m_gridWidth * m_gridHeight);
            m_gridPositions.resize(m_gridWidth * m_gridHeight);

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

            unsigned int ind = 0;
            for (int i = 0; i < m_gridHeight; i++) {
                for (int j = 0; j < m_gridWidth; j++) {
                    m_gridPositions[ind] = glm::vec2((az + (bz - az) * i / (m_gridWidth - 1)), (ax + (bx - ax) * j / (m_gridHeight - 1))); //i puis j donc z puis x
                    ind += 1;
                }
            }
        }

        //Verfier que bonne dim
        if (width != m_gridWidth || height != m_gridHeight) {
            std::cout << "ERROR : wrong size" << std::endl;
            exit(1);
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
    //for (int i = 0; i < height; i++) {
    //    for (int j = 0; j < width; j++) {
    //        m_vertexTexCoords.push_back(float(j) / (width - 1)); //x
    //        m_vertexTexCoords.push_back(float(i) / (height - 1)); //y
    //    }
    //}
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

void openglCallbackFunction(GLenum source,
    GLenum type,
    GLuint id,
    GLenum severity,
    GLsizei length,
    const GLchar* message,
    const void* userParam) {

    if (id == 131185) {
        return;
    }

    std::cout << "---------------------opengl-callback-start------------" << std::endl;
    std::cout << "message: " << message << std::endl;
    std::cout << "type: ";
    switch (type) {
    case GL_DEBUG_TYPE_ERROR:
        std::cout << "ERROR";
        break;
    case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
        std::cout << "DEPRECATED_BEHAVIOR";
        break;
    case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
        std::cout << "UNDEFINED_BEHAVIOR";
        break;
    case GL_DEBUG_TYPE_PORTABILITY:
        std::cout << "PORTABILITY";
        break;
    case GL_DEBUG_TYPE_PERFORMANCE:
        std::cout << "PERFORMANCE";
        break;
    case GL_DEBUG_TYPE_OTHER:
        std::cout << "OTHER";
        break;
    }
    std::cout << std::endl;

    std::cout << "id: " << id << std::endl;
    std::cout << "severity: ";
    switch (severity) {
    case GL_DEBUG_SEVERITY_LOW:
        std::cout << "LOW";
        break;
    case GL_DEBUG_SEVERITY_MEDIUM:
        std::cout << "MEDIUM";
        break;
    case GL_DEBUG_SEVERITY_HIGH:
        std::cout << "HIGH";
        break;
    }
    std::cout << std::endl;
    if (type != GL_DEBUG_TYPE_OTHER && id != 131218) {
        std::cout << "break point..." << std::endl;
        exit(1);
    }
    std::cout << "---------------------opengl-callback-end--------------" << std::endl;
}

void initOpenGL() {
    // Load extensions for modern OpenGL
    if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cerr << "ERROR: Failed to initialize OpenGL context" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
    }
  
    //Pour debug
    glEnable(GL_DEBUG_OUTPUT);
    glDebugMessageCallback(openglCallbackFunction, nullptr);
    GLuint unusedIds = 0;
    glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, &unusedIds, true);

    //glCullFace(GL_BACK); // specifies the faces to cull (here the ones pointing away from the camera)
    //glEnable(GL_CULL_FACE); // enables face culling (based on the orientation defined by the cw/ccw enumeration).
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Prend le alpha de l'eau et melange avec le fond en (1 - alpha)
    glDepthFunc(GL_LESS);   // Specify the depth test for the z-buffer
    glEnable(GL_DEPTH_TEST);      // Enable the z-buffer test in the rasterization
    glClearColor(0.7f, 0.7f, 0.7f, 1.0f); // specify the background color, used any time the framebuffer is cleared

    std::cout << "---INIT OPENGL DONE---" << std::endl;
}

// Loads the content of an ASCII file in a standard C++ string
std::string file2String(const std::string &filename) {
  std::ifstream t(filename.c_str());
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

//Test if there are errors in the bind shader
void CheckGlErrors(std::string caller) {
    unsigned int glerr = glGetError();
    if (glerr == GL_NO_ERROR)
        return;

    std::cout << "GL Error discovered from caller " << caller << std::endl;
    switch (glerr) {
    case GL_INVALID_ENUM:
        std::cout << "Invalid enum" << std::endl;
        break;
    case GL_INVALID_VALUE:
        std::cout << "Invalid value" << std::endl;
        break;
    case GL_INVALID_OPERATION:
        std::cout << "Invalid Operation" << std::endl;
        break;
    case GL_STACK_OVERFLOW:
        std::cout << "Stack overflow" << std::endl;
        break;
    case GL_STACK_UNDERFLOW:
        std::cout << "Stack underflow" << std::endl;
        break;
    case GL_OUT_OF_MEMORY:
        std::cout << "Out of memory" << std::endl;
        break;
    default:
        std::cout << "Unknown OpenGL error: " << glerr << std::endl;
    }
}

// Loads and compile a shader, before attaching it to a program
void loadShader(GLuint program, GLenum type, const std::string &shaderSourceString, const std::string &shaderName) {
    int status = 0;
    int logLength = 0;
    GLuint shader = glCreateShader(type); // Create the shader, e.g., a vertex shader to be applied to every single vertex of a mesh
    const GLchar *shaderSource = (const GLchar *)shaderSourceString.c_str(); // Interface the C++ string through a C pointer
    glShaderSource(shader, 1, &shaderSource, NULL); // load the vertex shader code
    glCompileShader(shader);
    CheckGlErrors(shaderName + " 1");

    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        std::cout << "Shader compilation failed : " << shaderName << std::endl;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);
        GLchar* log = new GLchar[logLength];
        glGetShaderInfoLog(shader, logLength, NULL, log);
        std::cout << "Log (len) " << logLength << " : " << log << std::endl;
        delete[] log;
        exit(1);
    }

    CheckGlErrors(shaderName + " 2");
    glAttachShader(program, shader);
    glDeleteShader(shader);
}

void initGPUprograms() {
    int status = 0;
    int logLength = 0;
    std::string  shaderSourceString;

    //Programm for Vertex and Fragment Shaders
    g_program = glCreateProgram(); // Create a GPU program, i.e., two central shaders of the graphics pipeline
    shaderSourceString = file2String("../vertexShader.glsl"); // Loads the shader source from a file to a C++ string
    loadShader(g_program, GL_VERTEX_SHADER, shaderSourceString, "Vertex Shader");
    shaderSourceString = file2String("../fragmentShader.glsl");
    loadShader(g_program, GL_FRAGMENT_SHADER, shaderSourceString, "Fragment Shader");
    glLinkProgram(g_program); // The main GPU program is ready to be handle streams of polygons
    CheckGlErrors("Shader Program 1");
    glGetProgramiv(g_program, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        std::cout << "Link failed Shader Program " << std::endl;
        glGetProgramiv(g_program, GL_INFO_LOG_LENGTH, &logLength);
        GLchar* log = new GLchar[logLength];
        glGetProgramInfoLog(g_program, logLength, NULL, log);
        std::cout << "Log (len) " << logLength << " : " << log << std::endl;
        delete[] log;
        exit(1);
    }
    CheckGlErrors("Shader Program 2");


    //Programm for compute shader
    g_computeProgram = glCreateProgram();
    shaderSourceString = file2String("../computeShader/thermalErosion.glsl");
    //Ajout de l'entête
    std::string computeShaderHeader = std::string("#version 430 \n #define NB_OF_LAYERS ") + std::to_string(NB_OF_LAYERS)
        + std::string("\n #define PATCH_HEIGHT ") + std::to_string(PATCH_HEIGHT)
        + std::string("\n #define PATCH_WIDTH ") + std::to_string(PATCH_WIDTH) + std::string("\n");
    shaderSourceString = computeShaderHeader + shaderSourceString;
    loadShader(g_computeProgram, GL_COMPUTE_SHADER, shaderSourceString, "Thermal Erosion");
    glLinkProgram(g_computeProgram);
    CheckGlErrors("Thermal Erosion 1");
    glGetProgramiv(g_computeProgram, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        std::cout << "Link failed Thermal Erosion" << std::endl;
        glGetProgramiv(g_computeProgram, GL_INFO_LOG_LENGTH, &logLength);
        GLchar* log = new GLchar[logLength];
        glGetProgramInfoLog(g_computeProgram, logLength, NULL, log);
        std::cout << "Log (len) " << logLength << " : " << log << std::endl;
        delete[] log;
        exit(1);
    }
    CheckGlErrors("Thermal Erosion 2");

    g_computeProgramHydraulicA = glCreateProgram();
    shaderSourceString = file2String("../computeShader/hydraulicErosionA.glsl");
    //Ajout de l'entête
    shaderSourceString = computeShaderHeader + shaderSourceString;
    loadShader(g_computeProgramHydraulicA, GL_COMPUTE_SHADER, shaderSourceString, "Hydraulic Erosion A");
    glLinkProgram(g_computeProgramHydraulicA);
    CheckGlErrors("Hydraulic Erosion A 1");
    glGetProgramiv(g_computeProgramHydraulicA, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        std::cout << "Link failed Hydraulic Erosion A" << std::endl;
        glGetProgramiv(g_computeProgramHydraulicA, GL_INFO_LOG_LENGTH, &logLength);
        GLchar* log = new GLchar[logLength];
        glGetProgramInfoLog(g_computeProgramHydraulicA, logLength, NULL, log);
        std::cout << "Log (len) " << logLength << " : " << log << std::endl;
        delete[] log;
        exit(1);
    }
    CheckGlErrors("Hydraulic Erosion A 2");

    g_computeProgramHydraulicB = glCreateProgram();
    shaderSourceString = file2String("../computeShader/hydraulicErosionBForward.glsl"); //hydraulicErosionBForward.glsl
    //Ajout de l'entête
    shaderSourceString = computeShaderHeader + shaderSourceString;
    loadShader(g_computeProgramHydraulicB, GL_COMPUTE_SHADER, shaderSourceString, "Hydraulic Erosion B");
    glLinkProgram(g_computeProgramHydraulicB);
    CheckGlErrors("Hydraulic Erosion B 1");
    glGetProgramiv(g_computeProgramHydraulicB, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        std::cout << "Link failed Hydraulic Erosion B" << std::endl;
        glGetProgramiv(g_computeProgramHydraulicB, GL_INFO_LOG_LENGTH, &logLength);
        GLchar* log = new GLchar[logLength];
        glGetProgramInfoLog(g_computeProgramHydraulicB, logLength, NULL, log);
        std::cout << "Log (len) " << logLength << " : " << log << std::endl;
        delete[] log;
        exit(1);
    }
    CheckGlErrors("Hydraulic Erosion B 2");

    //Program for computing normals, height, color of each point of the grid
    g_computeForRendering = glCreateProgram();
    shaderSourceString = file2String("../computeShader/computeForRendering.glsl");
    //Ajout de l'entête
    shaderSourceString = computeShaderHeader + shaderSourceString;
    loadShader(g_computeForRendering, GL_COMPUTE_SHADER, shaderSourceString, "Compute for rendering");
    glLinkProgram(g_computeForRendering);
    CheckGlErrors("Compute for rendering Shader 1");
    glGetProgramiv(g_computeForRendering, GL_LINK_STATUS, &status);
    if (status == GL_FALSE)
    {
        std::cout << "Link failed Compute for rendering terrain Shader" << std::endl;
        glGetProgramiv(g_computeForRendering, GL_INFO_LOG_LENGTH, &logLength);
        GLchar* log = new GLchar[logLength];
        glGetProgramInfoLog(g_computeForRendering, logLength, NULL, log);
        std::cout << "Log (len) " << logLength << " : " << log << std::endl;
        delete[] log;
        exit(1);
    }
    CheckGlErrors("Compute for rendering Shader 2");


    //Calcul des dims de l'espace d'invocation
    g_nbGroupsX = (mesh->getGridHeight() + PATCH_HEIGHT - 1) / PATCH_HEIGHT;
    g_nbGroupsY = (mesh->getGridWidth() + PATCH_WIDTH - 1) / PATCH_WIDTH;

    std::cout << "---INIT GPU PROGRAMS DONE---" << std::endl;
}

void initBuffersAndUniforms() {
    //Uniforms
    glUseProgram(g_computeProgram); //Il faut "bind" le program afin de pouvoir set des variables uniforme avec glUniform...
    GLint loc = glGetUniformLocation(g_computeProgram, "gridHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridHeight CP" << std::endl;
    glUniform1i(loc, mesh->getGridHeight());
    loc = glGetUniformLocation(g_computeProgram, "gridWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridWidth CP" << std::endl;
    glUniform1i(loc, mesh->getGridWidth());
    loc = glGetUniformLocation(g_computeProgram, "cellHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellHeight CP" << std::endl;
    glUniform1f(loc, mesh->getCellHeight());
    loc = glGetUniformLocation(g_computeProgram, "cellWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellWidth CP" << std::endl;
    glUniform1f(loc, mesh->getCellWidth());
    loc = glGetUniformLocation(g_computeProgram, "erosionCoeffs");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM erosionCoeffs CP" << std::endl;
    glUniform1fv(loc, NB_OF_LAYERS, mesh->getLayersErosionCoeffs().data());
    loc = glGetUniformLocation(g_computeProgram, "thetasLimit");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM thetasLimit CP" << std::endl;
    glUniform1fv(loc, NB_OF_LAYERS, mesh->getLayersThetaLimit().data());
    loc = glGetUniformLocation(g_computeProgram, "dt");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM dt CP" << std::endl;
    glUniform1f(loc, g_dt);
    glUseProgram(0);

    glUseProgram(g_computeProgramHydraulicA);
    loc = glGetUniformLocation(g_computeProgramHydraulicA, "gridHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridHeight CP HA" << std::endl;
    glUniform1i(loc, mesh->getGridHeight());
    loc = glGetUniformLocation(g_computeProgramHydraulicA, "gridWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridWidth CP HA" << std::endl;
    glUniform1i(loc, mesh->getGridWidth());
    loc = glGetUniformLocation(g_computeProgramHydraulicA, "cellHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellHeight CP HA" << std::endl;
    glUniform1f(loc, mesh->getCellHeight());
    loc = glGetUniformLocation(g_computeProgramHydraulicA, "cellWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellWidth CP HA" << std::endl;
    glUniform1f(loc, mesh->getCellWidth());
    loc = glGetUniformLocation(g_computeProgramHydraulicA, "dt");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM dt CP HA" << std::endl;
    glUniform1f(loc, g_dt);
    glUseProgram(0);

    
    glUseProgram(g_computeProgramHydraulicB);
    loc = glGetUniformLocation(g_computeProgramHydraulicB, "gridHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridHeight CP HB" << std::endl;
    glUniform1i(loc, mesh->getGridHeight());
    loc = glGetUniformLocation(g_computeProgramHydraulicB, "gridWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridWidth CP HB" << std::endl;
    glUniform1i(loc, mesh->getGridWidth());
    loc = glGetUniformLocation(g_computeProgramHydraulicB, "cellHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellHeight CP HB" << std::endl;
    glUniform1f(loc, mesh->getCellHeight());
    loc = glGetUniformLocation(g_computeProgramHydraulicB, "cellWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellWidth CP HB" << std::endl;
    glUniform1f(loc, mesh->getCellWidth());
    loc = glGetUniformLocation(g_computeProgramHydraulicB, "dt");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM dt CP HB" << std::endl;
    glUniform1f(loc, g_dt);
    glUseProgram(0);

    glUseProgram(g_computeForRendering); //Il faut "bind" le program afin de pouvoir set des variables uniforme avec glUniform...
    loc = glGetUniformLocation(g_computeForRendering, "gridHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridHeight R" << std::endl;
    glUniform1i(loc, mesh->getGridHeight());
    loc = glGetUniformLocation(g_computeForRendering, "gridWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM gridWidth R" << std::endl;
    glUniform1i(loc, mesh->getGridWidth());
    loc = glGetUniformLocation(g_computeForRendering, "cellHeight");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellHeight R" << std::endl;
    glUniform1f(loc, mesh->getCellHeight());
    loc = glGetUniformLocation(g_computeForRendering, "cellWidth");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM cellWidth R" << std::endl;
    glUniform1f(loc, mesh->getCellWidth());
    loc = glGetUniformLocation(g_computeForRendering, "layersColor");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM layersColor R" << std::endl;
    glUniform4fv(loc, NB_OF_LAYERS, mesh->getLayersColor().data());
    loc = glGetUniformLocation(g_computeForRendering, "displayVel");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM displayVel R" << std::endl;
    glUniform1i(loc, g_displayVel);
    glUseProgram(0);

    std::cout << "---INIT UNIFORMS DONE---" << std::endl;


    //Buffers :
    std::vector<float> gridPosTmp;
    for (unsigned int i = 0; i < mesh->getGridHeight(); i += 1) {
        for (unsigned int j = 0; j < mesh->getGridWidth(); j += 1) {
            gridPosTmp.push_back(mesh->getGridPos(i, j).x);
            gridPosTmp.push_back(mesh->getGridPos(i, j).y);
        }
    }
    size_t vertexBufferSize = sizeof(float) * gridPosTmp.size();
    glCreateBuffers(1, &g_gridPosVbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_gridPosVbo);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, gridPosTmp.data(), 0); //Le flag plutot 0 que GL_DYNAMIC_STORAGE_BIT ?

    std::vector<float> terrainPropsTmp;
    for (unsigned int i = 0; i < mesh->getGridHeight(); i += 1) {
        for (unsigned int j = 0; j < mesh->getGridWidth(); j += 1) {
            for (unsigned int k = 0; k < mesh->getNbOfLayers(); k += 1) {
                terrainPropsTmp.push_back(mesh->getLayerThickness(k, i, j));
            }
        }
    }
    vertexBufferSize = sizeof(float) * terrainPropsTmp.size();
    glCreateBuffers(1, &g_gridLayersHeightVboR);
    glBindBuffer(GL_ARRAY_BUFFER, g_gridLayersHeightVboR);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, terrainPropsTmp.data(), 0);

    for (unsigned int i = 0; i < terrainPropsTmp.size(); i += 1) {
        terrainPropsTmp[i] = 0.f;
    }
    vertexBufferSize = sizeof(float) * terrainPropsTmp.size();
    glCreateBuffers(1, &g_gridLayersHeightVboW);
    glBindBuffer(GL_ARRAY_BUFFER, g_gridLayersHeightVboW);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, terrainPropsTmp.data(), 0);

    std::vector<float> flowTmp(4 * mesh->getGridHeight() * mesh->getGridWidth(), 0);
    vertexBufferSize = sizeof(float) * flowTmp.size();
    glCreateBuffers(1, &g_waterFlowsVboR);
    glBindBuffer(GL_ARRAY_BUFFER, g_waterFlowsVboR);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, flowTmp.data(), 0);

    glCreateBuffers(1, &g_waterFlowsVboW);
    glBindBuffer(GL_ARRAY_BUFFER, g_waterFlowsVboW);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, flowTmp.data(), 0);

    vertexBufferSize = 2 * sizeof(float) * mesh->getGridHeight() * mesh->getGridWidth();
    glCreateBuffers(1, &g_waterVelocity);
    glBindBuffer(GL_ARRAY_BUFFER, g_waterVelocity);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, NULL, 0);

    vertexBufferSize = sizeof(float) * mesh->getGridHeight() * mesh->getGridWidth();
    glCreateBuffers(1, &g_waterSedimentR);
    glBindBuffer(GL_ARRAY_BUFFER, g_waterSedimentR);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, NULL, 0);

    glCreateBuffers(1, &g_waterSedimentW);
    glBindBuffer(GL_ARRAY_BUFFER, g_waterSedimentW);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, NULL, 0);

    vertexBufferSize = 4 * sizeof(float) * mesh->getGridHeight() * mesh->getGridWidth();
    glCreateBuffers(1, &g_gridNormalsVbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_gridNormalsVbo);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, NULL, 0);

    vertexBufferSize = sizeof(float) * mesh->getGridHeight() * mesh->getGridWidth();
    glCreateBuffers(1, &g_gridHeightsVbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_gridHeightsVbo);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, NULL, 0);

    vertexBufferSize = 4 * sizeof(float) * mesh->getGridHeight() * mesh->getGridWidth();
    glCreateBuffers(1, &g_gridColorsVbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_gridColorsVbo);
    glBufferStorage(GL_ARRAY_BUFFER, vertexBufferSize, NULL, 0);

    vertexBufferSize = sizeof(unsigned int) * mesh->getTriangleIndices().size();
    glCreateBuffers(1, &g_ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ibo);
    glBufferStorage(GL_ELEMENT_ARRAY_BUFFER, vertexBufferSize, mesh->getTriangleIndices().data(), 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    std::cout << "---INIT BUFFERS DONE---" << std::endl;
}

void setDisplayVelocity() {
    glUseProgram(g_computeForRendering);
    GLuint loc = glGetUniformLocation(g_computeForRendering, "displayVel");
    if (loc == -1) std::cout << "ERROR WITH UNIFORM displayVel R" << std::endl;
    glUniform1i(loc, g_displayVel);

    glUseProgram(0);
}

void initCamera() {
    int width, height;
    glfwGetWindowSize(g_window, &width, &height);
    g_camera.setAspectRatio(static_cast<float>(width)/static_cast<float>(height));

    g_camera.setPosition(glm::vec3(0.0, 0.0, 10.0));
    g_camera.setNear(0.1);
    g_camera.setFar(20.0);

    std::cout << "---INIT CAMERA DONE---" << std::endl;
}

void render();
void clear();

void renderImGui() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::Begin("Control Widget");
    ImGui::SetWindowSize(ImVec2(0, 0));
    static float layer1_col = 0.0;
    //ImGui::SliderFloat("Layer 1 color", &layer1_col, 0, 255.f);
    static float layer2_col = 0.0;
    //ImGui::SliderFloat("Layer 2 color", &layer2_col, 0, 255.f);
    static float translation[] = {0.0, 0.0};
    //ImGui::SliderFloat2("position", translation, -1.0, 1.0);

    ImGui::Separator();

    ImGui::Text("Erosion:");

    ImGui::Spacing();
    ImGui::Spacing();

    static float thetaLimit = 0.3;
    static float erosionCoeff = 0.3;
    static float dt = 0.0001;

    ImGui::SliderFloat("Theta", &thetaLimit, 0, PI/2);
    ImGui::SliderFloat("Erosion coefficient", &erosionCoeff, 0, 1);
    ImGui::SliderFloat("Dt", &dt, 0, 1); //Faudrait plutot 0.000001 à 0.1 genre juste les puissance de 10 si possible

    //if (ImGui::Button("Start thermal erosion")) {
    //    mesh->thermalErosion(thetaLimit, erosionCoeff, dt);
    //}

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

    if (ImGui::TreeNode("Range Widgets"))
    {
        static float begin = 10, end = 90;
        static int begin_i = 100, end_i = 1000;
        ImGui::DragFloatRange2("range float", &begin, &end, 0.25f, 0.0f, 100.0f, "Min: %.1f %%", "Max: %.1f %%", ImGuiSliderFlags_AlwaysClamp);
        ImGui::DragIntRange2("range int", &begin_i, &end_i, 5, 0, 1000, "Min: %d units", "Max: %d units");
        ImGui::DragIntRange2("range int (no bounds)", &begin_i, &end_i, 5, 0, 0, "Min: %d units", "Max: %d units");
        ImGui::TreePop();
    }

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
    ImGui_ImplOpenGL3_Init("#version 130");
    // Setup Dear ImGui style
    ImGui::StyleColorsDark();

    std::cout << "---INIT IMGUI DONE---" << std::endl;
}

void init() {
    initGLFW();
    initOpenGL();

    //mesh = new Mesh({ "../data/bedrock.png", "../data/sand.png", "../data/bedrock2.png" }, 
    //    { glm::vec4(120.f / 255.f, 135.f / 255.f, 124.f / 255.f, 1.f), glm::vec4(148.f / 255.f, 124.f / 255.f, 48.f / 255.f, 1.f), glm::vec4(20.f / 255.f, 107.f / 255.f, 150.f / 255.f, 0.5f) },
    //    { 0.f , 0.3f, 0.f},
    //    { 0.f, 0.3f, 0.f},
    //    glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 2.f)); //cpu
    mesh = new Mesh({ "../data/simpleB.png", "../data/sand-with-water.png", "../data/water-around-sand.png" },
        { glm::vec4(120.f / 255.f, 135.f / 255.f, 124.f / 255.f, 1.f), glm::vec4(148.f / 255.f, 124.f / 255.f, 48.f / 255.f, 1.f), glm::vec4(39.f / 255.f, 112.f / 255.f, 125.f / 255.f, 0.2) },
        { 0.f , 0.05f, 0.f },
        { 0.f, 0.3f, 0.f },
        glm::vec4(-5.f, -5.f, 5.f, 5.f), glm::vec2(0.f, 2.f)); //cpu
                                                               
                                                               
    //g_sunID = loadTextureFromFileToGPU("../data/heightmap3.jpg");
    initGPUprograms(); //Init aussi les dimension de l'espace d'invocation
    initBuffersAndUniforms(); //Aprs gpu programs car fait aussi uniform
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
    const glm::vec3 camPosition = g_camera.getPosition();

    glUniformMatrix4fv(glGetUniformLocation(g_program, "viewMat"), 1, GL_FALSE, glm::value_ptr(viewMatrix)); // compute the view matrix of the camera and pass it to the GPU program
    glUniformMatrix4fv(glGetUniformLocation(g_program, "projMat"), 1, GL_FALSE, glm::value_ptr(projMatrix)); // compute the projection matrix of the camera and pass it to the GPU program
    glUniform3f(glGetUniformLocation(g_program, "camPos"), camPosition[0], camPosition[1], camPosition[2]);

    mesh->render();
}

int main(int argc, char ** argv) {
    init(); // Your initialization code (user interface, OpenGL states, scene with geometry, material, lights, etc)
    
    GLint loc = 0;
    GLuint swapBuffT = 0;
    GLuint swapBuffF = 0;
    GLuint swapBuffS = 0;
    unsigned int nbOfItT = 500000;
    unsigned int nbOfItH = 500000;
    while(!glfwWindowShouldClose(g_window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Erase the color and z buffers.

        if (nbOfItT > 0) {
            //Phase de calculs :
            glUseProgram(g_computeProgram);
            //Bind les buffer du compute shader :
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_gridLayersHeightVboR);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_gridLayersHeightVboW);
            //Appelle au compute shader
            glDispatchCompute(g_nbGroupsX, g_nbGroupsY, 1); //Dimension 2D de l'espace d'invocation x correspond à i et y à j
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

            //Swap les buffer R et W
            swapBuffT = g_gridLayersHeightVboR;
            g_gridLayersHeightVboR = g_gridLayersHeightVboW;
            g_gridLayersHeightVboW = swapBuffT;

            nbOfItT -= 1;
        }

        if (nbOfItH > 0) {
            //Etape A
            //Phase de calculs :
            glUseProgram(g_computeProgramHydraulicA);
            //Bind les buffer du compute shader :
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_gridLayersHeightVboR);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_gridLayersHeightVboW);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, g_waterFlowsVboR);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, g_waterFlowsVboW);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, g_waterVelocity);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, g_waterSedimentR);
            //Appelle au compute shader
            glDispatchCompute(g_nbGroupsX, g_nbGroupsY, 1); //Dimension 2D de l'espace d'invocation x correspond à i et y à j
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

            //Swap les buffer R et W
            swapBuffT = g_gridLayersHeightVboR;
            g_gridLayersHeightVboR = g_gridLayersHeightVboW;
            g_gridLayersHeightVboW = swapBuffT;

            swapBuffF = g_waterFlowsVboR;
            g_waterFlowsVboR = g_waterFlowsVboW;
            g_waterFlowsVboW = swapBuffF;

            //Etape B:
            glUseProgram(g_computeProgramHydraulicB);
            //Bind les buffer du compute shader :
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_waterSedimentR);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_waterSedimentW);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, g_waterVelocity);
            
            //Appelle au compute shader
            glDispatchCompute(g_nbGroupsX, g_nbGroupsY, 1); //Dimension 2D de l'espace d'invocation x correspond à i et y à j
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

            //Swap les buffer R et W
            swapBuffS = g_waterSedimentR;
            g_waterSedimentR = g_waterSedimentW;
            g_waterSedimentW = swapBuffS;

            nbOfItH -= 1;
        }


        //Phase de rendering du terrain:
        //Update avant render (normales et hauteur)
        glUseProgram(g_computeForRendering);
        loc = glGetUniformLocation(g_computeForRendering, "renderWater");
        glUniform1i(loc, false);
        //Bind les buffer du compute shader :
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_gridLayersHeightVboR);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_gridNormalsVbo);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, g_gridHeightsVbo);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, g_gridColorsVbo);
        //Appelle au compute shader
        glDispatchCompute(g_nbGroupsX, g_nbGroupsY, 1); //Dimension 2D de l'espace d'invocation x correspond à i et y à j
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        glUseProgram(g_program);
        render();

        ////Phase de rendering de l'eau:
        ////Update avant render (normales et hauteur)
        glUseProgram(g_computeForRendering);
        loc = glGetUniformLocation(g_computeForRendering, "renderWater");
        glUniform1i(loc, true);
        //Bind les buffer du compute shader :
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, g_gridLayersHeightVboR);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, g_gridNormalsVbo);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, g_gridHeightsVbo);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, g_gridColorsVbo);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, g_waterVelocity);
        //Appelle au compute shader
        glDispatchCompute(g_nbGroupsX, g_nbGroupsY, 1); //Dimension 2D de l'espace d'invocation x correspond à i et y à j
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        glUseProgram(g_program);
        render();

        //Render Imgui
        renderImGui();


        glfwSwapBuffers(g_window);
        glfwPollEvents();

        //Sleep(500);
    }

    // Cleanup
    clear();
    return EXIT_SUCCESS;
}

