#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <ModelTriangle.h>

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

// 全局数据
std::vector<ModelTriangle> modelTriangles;
// 深度缓冲：初始化为 -∞，表示“还没被任何东西画过”
std::vector<float> depthBuffer(WIDTH * HEIGHT, -std::numeric_limits<float>::infinity());

enum RenderMode {
    WIREFRAME = 0,
    RASTERIZED = 1
};

RenderMode renderMode = RASTERIZED;

glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
float focalLength = 2.0f;

// ----------------- 工具函数 -----------------

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
    std::ifstream file(filename);
    std::string line;
    std::string currentMaterial;

    if (!file.is_open()) {
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

    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
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

// 把 3D 顶点投影到屏幕坐标
CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 cameraPosition, float focalLength, glm::vec3 vertexPosition) {
    float x = vertexPosition.x - cameraPosition.x;
    float y = vertexPosition.y - cameraPosition.y;
    float z = cameraPosition.z - vertexPosition.z;

    // 简单 near-plane，防止除 0 和 depth 过大
    if (z < 0.1f) z = 0.1f;

    float u = focalLength * (x / z) * 160.0f + WIDTH * 0.5f;
    float v = HEIGHT - (focalLength * (y / z) * 160.0f + HEIGHT * 0.5f);

    float depth = 1.0f / z;  // z 越近，depth 越大

    return CanvasPoint(u, v, depth);
}

void drawStrokedTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
    drawLine(window, triangle[0], triangle[1], colour);
    drawLine(window, triangle[1], triangle[2], colour);
    drawLine(window, triangle[2], triangle[0], colour);
}

// 计算重心坐标的工具函数
glm::vec3 computeBarycentric(const CanvasPoint &a, const CanvasPoint &b, const CanvasPoint &c, float px, float py) {
    float denom = ((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
    if (std::abs(denom) < 1e-6f) {
        // 退化三角形
        return glm::vec3(-1.0f, -1.0f, -1.0f);
    }

    float w0 = ((b.y - c.y) * (px - c.x) + (c.x - b.x) * (py - c.y)) / denom;
    float w1 = ((c.y - a.y) * (px - c.x) + (a.x - c.x) * (py - c.y)) / denom;
    float w2 = 1.0f - w0 - w1;

    return glm::vec3(w0, w1, w2);
}

// 使用重心坐标 + Z-buffer 填充三角形（修复毛刺）
void drawFilledTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
    CanvasPoint p0 = triangle[0];
    CanvasPoint p1 = triangle[1];
    CanvasPoint p2 = triangle[2];

    uint32_t fillColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;

    // 计算包围盒
    float minXf = std::min({ p0.x, p1.x, p2.x });
    float maxXf = std::max({ p0.x, p1.x, p2.x });
    float minYf = std::min({ p0.y, p1.y, p2.y });
    float maxYf = std::max({ p0.y, p1.y, p2.y });

    int minX = std::max(0, static_cast<int>(std::floor(minXf)));
    int maxX = std::min(static_cast<int>(window.width) - 1, static_cast<int>(std::ceil(maxXf)));
    int minY = std::max(0, static_cast<int>(std::floor(minYf)));
    int maxY = std::min(static_cast<int>(window.height) - 1, static_cast<int>(std::ceil(maxYf)));

    // 防止越界
    if (minX > maxX || minY > maxY) return;

    for (int y = minY; y <= maxY; ++y) {
        for (int x = minX; x <= maxX; ++x) {
            // 用像素中心 (x+0.5, y+0.5) 计算
            float px = x + 0.5f;
            float py = y + 0.5f;

            glm::vec3 bary = computeBarycentric(p0, p1, p2, px, py);
            const float eps = -1e-4f; // 允许一点点负数误差

            if (bary.x >= eps && bary.y >= eps && bary.z >= eps) {
                // 插值深度
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
        CanvasPoint p0 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, modelTriangle.vertices[0]);
        CanvasPoint p1 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, modelTriangle.vertices[1]);
        CanvasPoint p2 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, modelTriangle.vertices[2]);

        CanvasTriangle canvasTriangle(p0, p1, p2);

        drawStrokedTriangle(window, canvasTriangle, white);
    }
}

void drawRasterized(DrawingWindow &window) {
    for (const auto &modelTriangle : modelTriangles) {
        CanvasPoint p0 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, modelTriangle.vertices[0]);
        CanvasPoint p1 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, modelTriangle.vertices[1]);
        CanvasPoint p2 = projectVertexOntoCanvasPoint(cameraPosition, focalLength, modelTriangle.vertices[2]);

        CanvasTriangle canvasTriangle(p0, p1, p2);

        drawFilledTriangle(window, canvasTriangle, modelTriangle.colour);
    }
}

void clearDepthBuffer() {
    std::fill(depthBuffer.begin(), depthBuffer.end(), -std::numeric_limits<float>::infinity());
}

void draw(DrawingWindow &window) {
    window.clearPixels();
    clearDepthBuffer();
    if (renderMode == WIREFRAME) {
        drawWireframe(window);
    } else {
        drawRasterized(window);
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
    }
    else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
    return false;
}

bool updateCamera(float deltaTime) {
    const Uint8 *keyState = SDL_GetKeyboardState(NULL);
    const float cameraSpeed = 2.0f * deltaTime; // 每秒2个单位
    bool cameraMoved = false;

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

    return cameraMoved;
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
    const float maxDeltaTime = 1.0f / 30.0f; // 限制最大deltaTime为30FPS

    while (true) {
        Uint32 currentTime = SDL_GetTicks();
        float deltaTime = (currentTime - lastFrameTime) / 1000.0f;
        deltaTime = std::min(deltaTime, maxDeltaTime); // 限制最大deltaTime

        // 处理窗口事件
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

        // 更新相机位置（基于键盘状态）
        updateCamera(deltaTime);

        // 每帧都渲染
        draw(window);
        window.renderFrame();

        // 帧率控制
        Uint32 elapsed = SDL_GetTicks() - currentTime;
        if (elapsed < frameTime) {
            SDL_Delay(static_cast<Uint32>(frameTime - elapsed));
        }

        lastFrameTime = currentTime;
    }
}
