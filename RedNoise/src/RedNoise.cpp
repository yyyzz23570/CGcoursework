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

#define WIDTH 640
#define HEIGHT 480

std::vector<ModelTriangle> modelTriangles;
std::vector<float> depthBuffer(WIDTH * HEIGHT, std::numeric_limits<float>::max());

enum RenderMode {
	WIREFRAME = 0,
	RASTERIZED = 1
};

RenderMode renderMode = RASTERIZED;

glm::vec3 cameraPosition(0.0f, 0.0f, 4.0f);
float focalLength = 2.0f;

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

std::vector<ModelTriangle> loadOBJFile(const std::string &filename, const std::string &mtlPath, float scaleFactor = 1.0f) {
	std::vector<ModelTriangle> triangles;
	std::vector<glm::vec3> vertices;
	std::map<std::string, Colour> materials = loadMTLFile(mtlPath);
	Colour currentColour(255, 255, 255);
	
	std::ifstream file(filename);
	std::string line;
	
	if (!file.is_open()) {
		return triangles;
	}
	
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string token;
		iss >> token;
		
		if (token == "v") {
			float x, y, z;
			iss >> x >> y >> z;
			vertices.push_back(glm::vec3(x * scaleFactor, y * scaleFactor, z * scaleFactor));
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
	glm::vec3 vertexInCameraSpace = vertexPosition - cameraPosition;
	
	float x = vertexInCameraSpace.x;
	float y = vertexInCameraSpace.y;
	float z = -vertexInCameraSpace.z;  // 反转 z 轴，使相机朝向 -z 方向
	
	if (z <= 0) {
		z = 0.0001f;
	}
	
	float u = focalLength * (x / z) + WIDTH * 0.5f;
	float v = focalLength * (y / z) + HEIGHT * 0.5f;
	
	float imagePlaneScale = 160.0f;
	
	float offsetX = (u - WIDTH * 0.5f) * imagePlaneScale;
	float offsetY = (v - HEIGHT * 0.5f) * imagePlaneScale;
	
	u = offsetX + WIDTH * 0.5f;
	v = offsetY + HEIGHT * 0.5f;
	
	v = HEIGHT - v;
	
	return CanvasPoint(u, v, z);
}

void drawStrokedTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
	drawLine(window, triangle[0], triangle[1], colour);
	drawLine(window, triangle[1], triangle[2], colour);
	drawLine(window, triangle[2], triangle[0], colour);
}

void drawFilledTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
	CanvasPoint p0 = triangle[0];
	CanvasPoint p1 = triangle[1];
	CanvasPoint p2 = triangle[2];
	
	if (p0.y > p1.y) std::swap(p0, p1);
	if (p0.y > p2.y) std::swap(p0, p2);
	if (p1.y > p2.y) std::swap(p1, p2);
	
	uint32_t fillColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
	
	int startY = (int)round(p0.y);
	int endY = (int)round(p2.y);
	
	for (int y = startY; y <= endY; y++) {
		float t;
		float leftX, rightX;
		float leftDepth, rightDepth;
		
		if (y <= p1.y) {
			if (p1.y != p0.y) {
				t = (y - p0.y) / (p1.y - p0.y);
				leftX = p0.x + t * (p1.x - p0.x);
				leftDepth = p0.depth + t * (p1.depth - p0.depth);
			} else {
				leftX = p0.x;
				leftDepth = p0.depth;
			}
			
			if (p2.y != p0.y) {
				t = (y - p0.y) / (p2.y - p0.y);
				rightX = p0.x + t * (p2.x - p0.x);
				rightDepth = p0.depth + t * (p2.depth - p0.depth);
			} else {
				rightX = p0.x;
				rightDepth = p0.depth;
			}
		} else {
			if (p2.y != p1.y) {
				t = (y - p1.y) / (p2.y - p1.y);
				leftX = p1.x + t * (p2.x - p1.x);
				leftDepth = p1.depth + t * (p2.depth - p1.depth);
			} else {
				leftX = p1.x;
				leftDepth = p1.depth;
			}
			
			if (p2.y != p0.y) {
				t = (y - p0.y) / (p2.y - p0.y);
				rightX = p0.x + t * (p2.x - p0.x);
				rightDepth = p0.depth + t * (p2.depth - p0.depth);
			} else {
				rightX = p0.x;
				rightDepth = p0.depth;
			}
		}
		
		if (leftX > rightX) {
			std::swap(leftX, rightX);
			std::swap(leftDepth, rightDepth);
		}
		
		int startX = (int)round(leftX);
		int endX = (int)round(rightX);
		
		for (int x = startX; x <= endX; x++) {
			if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
				size_t pixelIndex = y * window.width + x;
				
				float depth;
				if (endX != startX) {
					float tX = (x - leftX) / (rightX - leftX);
					depth = leftDepth + tX * (rightDepth - leftDepth);
				} else {
					depth = leftDepth;
				}
				
				if (depth < depthBuffer[pixelIndex] && depth > 0) {
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
	std::fill(depthBuffer.begin(), depthBuffer.end(), std::numeric_limits<float>::max());
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

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_1) {
			renderMode = WIREFRAME;
		}
		else if (event.key.keysym.sym == SDLK_2) {
			renderMode = RASTERIZED;
		}
		else if (event.key.keysym.sym == SDLK_LEFT) {
			cameraPosition.x -= 0.1f;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			cameraPosition.x += 0.1f;
		}
		else if (event.key.keysym.sym == SDLK_UP) {
			cameraPosition.y += 0.1f;
		}
		else if (event.key.keysym.sym == SDLK_DOWN) {
			cameraPosition.y -= 0.1f;
		}
		else if (event.key.keysym.sym == SDLK_w) {
			cameraPosition.z -= 0.1f;
		}
		else if (event.key.keysym.sym == SDLK_s) {
			cameraPosition.z += 0.1f;
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	
	float scaleFactor = 0.35f;
	modelTriangles = loadOBJFile("cornell-box.obj", "cornell-box.mtl", scaleFactor);
	
	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		window.renderFrame();
	}
}
