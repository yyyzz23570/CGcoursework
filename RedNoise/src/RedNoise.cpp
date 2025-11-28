#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>

#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>
#include <sstream>
#include <glm/glm.hpp>
#include <limits>
#include <SDL.h>

#define WIDTH 640
#define HEIGHT 480

std::vector<ModelTriangle> modelTriangles;
std::vector<float> depthBuffer(WIDTH * HEIGHT, -std::numeric_limits<float>::infinity());

enum RenderMode {
    WIREFRAME = 0,
    RASTERIZED = 1,
    RAYTRACED_SHADOW = 2
};

RenderMode renderMode = RASTERIZED;

glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
float focalLength = 3.5f;

float sceneRotationY = 0.0f;
glm::vec3 sceneCenter(0.0f, 0.0f, 0.0f);
glm::vec3 lightPosition(0.0f, 0.9f, 0.0f);
glm::mat3 rotateY(float angle) {
    float c = cos(angle);
    float s = sin(angle);
    return glm::mat3(
        c,  0.0f, s,
        0.0f, 1.0f, 0.0f,
        -s, 0.0f, c
    );
}

glm::vec3 rotateVertexAroundPoint(const glm::vec3& vertex, const glm::vec3& center, float angleY) {
    glm::vec3 relativePos = vertex - center;
    glm::mat3 rotationMatrix = rotateY(angleY);
    glm::vec3 rotatedRelative = rotationMatrix * relativePos;
    return center + rotatedRelative;
}

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
    std::vector<float> result;
    if (numberOfValues <= 1) {
        result.push_back(from);
        return result;
    }

    float step = (to - from) / (numberOfValues - 1);
    for (int i = 0; i < numberOfValues; i++) {
        result.push_back(from + i * step);
    }
    return result;
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;

    int numberOfSteps = static_cast<int>(std::max(std::abs(xDiff), std::abs(yDiff))) + 1;
    if (numberOfSteps <= 0) {
        numberOfSteps = 1;
    }

    std::vector<float> xValues = interpolateSingleFloats(from.x, to.x, numberOfSteps);
    std::vector<float> yValues = interpolateSingleFloats(from.y, to.y, numberOfSteps);

    uint32_t pixelColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;

    for (int i = 0; i < numberOfSteps; i++) {
        int x = round(xValues[i]);
        int y = round(yValues[i]);
        if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
            window.setPixelColour(x, y, pixelColour);
        }
    }
}

std::map<std::string, Colour> loadMTLFile(const std::string &filename) {
    std::map<std::string, Colour> materials;
    std::ifstream file;
    std::string line;
    std::string currentMaterial;

    std::vector<std::string> paths = {
        filename,
        "../" + filename,
        "../../" + filename,
        "./" + filename
    };

    bool opened = false;
    for (const auto &path : paths) {
        file.open(path);
        if (file.is_open()) {
            opened = true;
            break;
        }
    }

    if (!opened) {
        std::cerr << "Failed to open MTL: " << filename << std::endl;
        return materials;
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "newmtl") {
            iss >> currentMaterial;
        } else if (token == "Kd" && !currentMaterial.empty()) {
            float r, g, b;
            iss >> r >> g >> b;
            materials[currentMaterial] = Colour(
                currentMaterial,
                static_cast<int>(r * 255),
                static_cast<int>(g * 255),
                static_cast<int>(b * 255)
            );
        }
    }

    file.close();
    return materials;
}

std::vector<ModelTriangle> loadOBJFile(const std::string &filename, const std::string &mtlPath, float scaleFactor = 1.0f, glm::vec3 translate = glm::vec3(0.0f, 0.0f, 0.0f)) {
    std::vector<ModelTriangle> triangles;
    std::vector<glm::vec3> vertices;
    std::map<std::string, Colour> materials = loadMTLFile(mtlPath);
    Colour currentColour(255, 255, 255);

    std::ifstream file;
    std::string line;

    std::vector<std::string> paths = {
        filename,
        "../" + filename,
        "../../" + filename,
        "./" + filename
    };

    bool opened = false;
    for (const auto &path : paths) {
        file.open(path);
        if (file.is_open()) {
            opened = true;
            break;
        }
    }

    if (!opened) {
        std::cerr << "Failed to open OBJ: " << filename << std::endl;
        return triangles;
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "v") {
            float x, y, z;
            iss >> x >> y >> z;
            glm::vec3 vertex = glm::vec3(x * scaleFactor, y * scaleFactor, z * scaleFactor) + translate;
            vertices.push_back(vertex);
        } else if (token == "usemtl") {
            std::string materialName;
            iss >> materialName;
            if (materials.find(materialName) != materials.end()) {
                currentColour = materials[materialName];
            }
        } else if (token == "f") {
            std::vector<int> vertexIndices;

            std::string vertexStr;
            std::string remainingLine = (line.length() > 1) ? line.substr(1) : "";
            std::istringstream faceIss(remainingLine);

            while (faceIss >> vertexStr) {
                if (vertexStr.empty()) continue;

                size_t slashPos = vertexStr.find('/');
                if (slashPos != std::string::npos) {
                    vertexStr = vertexStr.substr(0, slashPos);
                }

                if (!vertexStr.empty()) {
                    try {
                        int idx = std::stoi(vertexStr) - 1;
                        if (idx >= 0 && idx < static_cast<int>(vertices.size())) {
                            vertexIndices.push_back(idx);
                        }
                    } catch (const std::exception &e) {
                        continue;
                    }
                }
            }

            if (vertexIndices.size() >= 3) {
                for (size_t i = 1; i + 1 < vertexIndices.size(); i++) {
                    ModelTriangle triangle(
                        vertices[vertexIndices[0]],
                        vertices[vertexIndices[i]],
                        vertices[vertexIndices[i + 1]],
                        currentColour
                    );
                    triangles.push_back(triangle);
                }
            }
        }
    }

    file.close();
    return triangles;
}

CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 cameraPosition, float focalLength, glm::vec3 vertexPosition) {
    float x = vertexPosition.x - cameraPosition.x;
    float y = vertexPosition.y - cameraPosition.y;
    float z = cameraPosition.z - vertexPosition.z;

    if (z < 0.1f) z = 0.1f;

    float u = focalLength * (x / z) * 160.0f + WIDTH * 0.5f;
    float v = HEIGHT - (focalLength * (y / z) * 160.0f + HEIGHT * 0.5f);
    float depth = 1.0f / z;

    return CanvasPoint(u, v, depth);
}

void drawStrokedTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
    drawLine(window, triangle[0], triangle[1], colour);
    drawLine(window, triangle[1], triangle[2], colour);
    drawLine(window, triangle[2], triangle[0], colour);
}

glm::vec3 computeBarycentric(const CanvasPoint &a, const CanvasPoint &b, const CanvasPoint &c, float px, float py) {
    float denom = ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
    if (std::abs(denom) < 1e-6f) {
        return glm::vec3(-1.0f, -1.0f, -1.0f);
    }

    float w0 = ((b.y - c.y) * (px - c.x) + (c.x - b.x) * (py - c.y)) / denom;
    float w1 = ((c.y - a.y) * (px - c.x) + (a.x - c.x) * (py - c.y)) / denom;
    float w2 = 1.0f - w0 - w1;

    return glm::vec3(w0, w1, w2);
}

void drawFilledTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
    CanvasPoint p0 = triangle[0];
    CanvasPoint p1 = triangle[1];
    CanvasPoint p2 = triangle[2];

    uint32_t fillColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;

    float minXf = std::min({ p0.x, p1.x, p2.x });
    float maxXf = std::max({ p0.x, p1.x, p2.x });
    float minYf = std::min({ p0.y, p1.y, p2.y });
    float maxYf = std::max({ p0.y, p1.y, p2.y });

    int minX = std::max(0, static_cast<int>(std::floor(minXf)));
    int maxX = std::min(static_cast<int>(window.width) - 1, static_cast<int>(std::ceil(maxXf)));
    int minY = std::max(0, static_cast<int>(std::floor(minYf)));
    int maxY = std::min(static_cast<int>(window.height) - 1, static_cast<int>(std::ceil(maxYf)));

    if (minX > maxX || minY > maxY) return;

    for (int y = minY; y <= maxY; ++y) {
        for (int x = minX; x <= maxX; ++x) {
            float px = x + 0.5f;
            float py = y + 0.5f;

            glm::vec3 bary = computeBarycentric(p0, p1, p2, px, py);
            const float eps = -1e-4f;

            if (bary.x >= eps && bary.y >= eps && bary.z >= eps) {
                float depth = bary.x * p0.depth + bary.y * p1.depth + bary.z * p2.depth;

                int pixelIndex = y * window.width + x;
                if (depth > depthBuffer[pixelIndex]) {
                    depthBuffer[pixelIndex] = depth;
                    window.setPixelColour(x, y, fillColour);
                }
            }
        }
    }
}

void drawWireframe(DrawingWindow &window) {
    Colour white(255, 255, 255);

    for (const auto &modelTriangle : modelTriangles) {
        glm::vec3 v0 = rotateVertexAroundPoint(modelTriangle.vertices[0], sceneCenter, sceneRotationY);
        glm::vec3 v1 = rotateVertexAroundPoint(modelTriangle.vertices[1], sceneCenter, sceneRotationY);
        glm::vec3 v2 = rotateVertexAroundPoint(modelTriangle.vertices[2], sceneCenter, sceneRotationY);
        
        CanvasPoint p0 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, v0);
        CanvasPoint p1 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, v1);
        CanvasPoint p2 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, v2);

        CanvasTriangle canvasTriangle(p0, p1, p2);
        drawStrokedTriangle(window, canvasTriangle, white);
    }
}

void drawRasterized(DrawingWindow &window) {
    for (const auto &modelTriangle : modelTriangles) {
        glm::vec3 v0 = rotateVertexAroundPoint(modelTriangle.vertices[0], sceneCenter, sceneRotationY);
        glm::vec3 v1 = rotateVertexAroundPoint(modelTriangle.vertices[1], sceneCenter, sceneRotationY);
        glm::vec3 v2 = rotateVertexAroundPoint(modelTriangle.vertices[2], sceneCenter, sceneRotationY);
        
        CanvasPoint p0 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, v0);
        CanvasPoint p1 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, v1);
        CanvasPoint p2 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, v2);

        CanvasTriangle canvasTriangle(p0, p1, p2);
        drawFilledTriangle(window, canvasTriangle, modelTriangle.colour);
    }
}

void clearDepthBuffer() {
    std::fill(depthBuffer.begin(), depthBuffer.end(), -std::numeric_limits<float>::infinity());
}

glm::vec3 generateRayDirection(int pixelX, int pixelY) {
    float scale = 160.0f;

    float x = (pixelX + 0.5f - WIDTH * 0.5f) / scale;
    float y = -(pixelY + 0.5f - HEIGHT * 0.5f) / scale;

    glm::vec3 dir(x, y, -focalLength);
    return glm::normalize(dir);
}

RayTriangleIntersection getClosestValidIntersection(const glm::vec3 &origin,
                                                    const glm::vec3 &rayDirection,
                                                    int excludeIndex = -1) {
    RayTriangleIntersection closest;
    closest.distanceFromCamera = std::numeric_limits<float>::infinity();
    closest.triangleIndex = -1;

    glm::vec3 d = glm::normalize(rayDirection);

    for (size_t i = 0; i < modelTriangles.size(); i++) {
        if (excludeIndex >= 0 && static_cast<int>(i) == excludeIndex) {
            continue;
        }
        const ModelTriangle &tri = modelTriangles[i];

        glm::vec3 v0 = rotateVertexAroundPoint(tri.vertices[0], sceneCenter, sceneRotationY);
        glm::vec3 v1 = rotateVertexAroundPoint(tri.vertices[1], sceneCenter, sceneRotationY);
        glm::vec3 v2 = rotateVertexAroundPoint(tri.vertices[2], sceneCenter, sceneRotationY);

        glm::vec3 e0 = v1 - v0;
        glm::vec3 e1 = v2 - v0;
        glm::vec3 SPVector = origin - v0;

        glm::mat3 DEMatrix(-d, e0, e1);
        float det = glm::determinant(DEMatrix);
        if (std::abs(det) < 0.0001f) continue;

        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
        float t = possibleSolution.x;
        float u = possibleSolution.y;
        float v = possibleSolution.z;

        if (t > 0.001f &&
            u >= 0.0f && u <= 1.0f &&
            v >= 0.0f && v <= 1.0f &&
            (u + v) <= 1.0f) {

            if (t < closest.distanceFromCamera) {
                glm::vec3 intersectionPoint = origin + t * d;

                ModelTriangle rotatedTri = tri;
                rotatedTri.vertices[0] = v0;
                rotatedTri.vertices[1] = v1;
                rotatedTri.vertices[2] = v2;
                glm::vec3 n = glm::normalize(glm::cross(e1, e0));
                rotatedTri.normal = n;

                closest = RayTriangleIntersection(intersectionPoint, t, rotatedTri, i);
            }
        }
    }

    return closest;
}

void drawRayTracedHardShadow(DrawingWindow &window) {
    float ambient = 0.3f;

    for (int y = 0; y < window.height; ++y) {
        for (int x = 0; x < window.width; ++x) {
            glm::vec3 rayDir = generateRayDirection(x, y);
            RayTriangleIntersection hit = getClosestValidIntersection(cameraPosition, rayDir);

            if (hit.triangleIndex == -1) {
                continue;
            }

            glm::vec3 hitPoint = hit.intersectionPoint;
            glm::vec3 toLight = lightPosition - hitPoint;

            float distanceToLight = glm::length(toLight);
            if (distanceToLight < 0.001f) distanceToLight = 0.001f;

            glm::vec3 lightDir = glm::normalize(toLight);

            glm::vec3 shadowOrigin = hitPoint;
            glm::vec3 n = hit.intersectedTriangle.normal;
            if (glm::length(n) > 0.0f) {
                n = glm::normalize(n);
                shadowOrigin += n * 0.0001f;
            } else {
                shadowOrigin += lightDir * 0.0001f;
            }

            RayTriangleIntersection shadowHit =
                getClosestValidIntersection(shadowOrigin, lightDir, static_cast<int>(hit.triangleIndex));

            bool inShadow = false;
            if (shadowHit.triangleIndex != -1 &&
                shadowHit.distanceFromCamera < distanceToLight) {
                inShadow = true;
            }

            Colour base = hit.intersectedTriangle.colour;
            float brightness = inShadow ? ambient : 1.0f;

            int r = static_cast<int>(base.red   * brightness);
            int g = static_cast<int>(base.green * brightness);
            int b = static_cast<int>(base.blue  * brightness);
            r = glm::clamp(r, 0, 255);
            g = glm::clamp(g, 0, 255);
            b = glm::clamp(b, 0, 255);

            uint32_t colourValue =
                (255 << 24) + (r << 16) + (g << 8) + b;
            window.setPixelColour(x, y, colourValue);
        }
    }
}

void draw(DrawingWindow &window) {
    window.clearPixels();
    clearDepthBuffer();
    if (renderMode == WIREFRAME) {
        drawWireframe(window);
    } else if (renderMode == RASTERIZED) {
        drawRasterized(window);
    } else if (renderMode == RAYTRACED_SHADOW) {
        drawRayTracedHardShadow(window);
    }
}

bool handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_1) {
            renderMode = WIREFRAME;
            return true;
        }
        else if (event.key.keysym.sym == SDLK_2) {
            renderMode = RASTERIZED;
            return true;
        }
        else if (event.key.keysym.sym == SDLK_3) {
            renderMode = RAYTRACED_SHADOW;
            return true;
        }
    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    return false;
}

bool updateCamera(float deltaTime) {
    const Uint8 *keyState = SDL_GetKeyboardState(NULL);
    const float cameraSpeed = 2.0f * deltaTime;
    const float rotationSpeed = 2.0f * deltaTime;
    bool cameraMoved = false;
    bool objectRotated = false;

    if (keyState[SDL_SCANCODE_LEFT]) {
        cameraPosition.x -= cameraSpeed;
        cameraMoved = true;
    }
    if (keyState[SDL_SCANCODE_RIGHT]) {
        cameraPosition.x += cameraSpeed;
        cameraMoved = true;
    }
    if (keyState[SDL_SCANCODE_UP]) {
        cameraPosition.y += cameraSpeed;
        cameraMoved = true;
    }
    if (keyState[SDL_SCANCODE_DOWN]) {
        cameraPosition.y -= cameraSpeed;
        cameraMoved = true;
    }
    if (keyState[SDL_SCANCODE_W]) {
        cameraPosition.z -= cameraSpeed;
        cameraMoved = true;
    }
    if (keyState[SDL_SCANCODE_S]) {
        cameraPosition.z += cameraSpeed;
        cameraMoved = true;
    }
    if (keyState[SDL_SCANCODE_A]) {
        sceneRotationY += rotationSpeed;
        objectRotated = true;
    }
    if (keyState[SDL_SCANCODE_D]) {
        sceneRotationY -= rotationSpeed;
        objectRotated = true;
    }

    return cameraMoved || objectRotated;
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

    float boxScale = 0.35f;
    std::vector<ModelTriangle> boxTriangles = loadOBJFile("cornell-box.obj", "cornell-box.mtl", boxScale);
    modelTriangles.insert(modelTriangles.end(), boxTriangles.begin(), boxTriangles.end());

    float sphereScale = 0.25f;
    float floorY = -2.74f * boxScale;
    float sphereCenterY = 1.5f * sphereScale;
    float sphereRadius = 1.0f * sphereScale;
    float translateY = floorY - (sphereCenterY - sphereRadius);
    glm::vec3 sphereTranslate(-0.7f, translateY, 0.0f);
    std::vector<ModelTriangle> sphereTriangles = loadOBJFile("sphere.obj", "", sphereScale, sphereTranslate);

    Colour pink(255, 192, 203);
    for (auto &triangle : sphereTriangles) {
        triangle.colour = pink;
    }
    modelTriangles.insert(modelTriangles.end(), sphereTriangles.begin(), sphereTriangles.end());

    Uint32 lastFrameTime = SDL_GetTicks();
    const float targetFPS = 60.0f;
    const float frameTime = 1000.0f / targetFPS;
    const float maxDeltaTime = 1.0f / 30.0f;

    while (true) {
        Uint32 currentTime = SDL_GetTicks();
        float deltaTime = (currentTime - lastFrameTime) / 1000.0f;
        deltaTime = std::min(deltaTime, maxDeltaTime);

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                window.exitCleanly();
            }
            if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE) {
                window.exitCleanly();
            }
            handleEvent(event, window);
        }

        updateCamera(deltaTime);
        draw(window);
        window.renderFrame();

        Uint32 elapsed = SDL_GetTicks() - currentTime;
        if (elapsed < frameTime) {
            SDL_Delay(static_cast<Uint32>(frameTime - elapsed));
        }

        lastFrameTime = currentTime;
    }
}
