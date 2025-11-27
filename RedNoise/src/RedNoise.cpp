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
std::vector<float> depthBuffer(WIDTH * HEIGHT, 0.0f);

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
	float x = vertexPosition.x - cameraPosition.x;
	float y = vertexPosition.y - cameraPosition.y;
	float z = cameraPosition.z - vertexPosition.z;
	
	if (z <= 0) z = 0.0001f;
	
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
		float leftX, rightX, leftDepth, rightDepth;
		
		if (y <= p1.y) {
			float t1 = (p1.y != p0.y) ? (y - p0.y) / (p1.y - p0.y) : 0;
			leftX = p0.x + t1 * (p1.x - p0.x);
			leftDepth = p0.depth + t1 * (p1.depth - p0.depth);
			
			float t2 = (p2.y != p0.y) ? (y - p0.y) / (p2.y - p0.y) : 0;
			rightX = p0.x + t2 * (p2.x - p0.x);
			rightDepth = p0.depth + t2 * (p2.depth - p0.depth);
		}
		else {
			float t1 = (p2.y != p1.y) ? (y - p1.y) / (p2.y - p1.y) : 0;
			leftX = p1.x + t1 * (p2.x - p1.x);
			leftDepth = p1.depth + t1 * (p2.depth - p1.depth);
			
			float t2 = (p2.y != p0.y) ? (y - p0.y) / (p2.y - p0.y) : 0;
			rightX = p0.x + t2 * (p2.x - p0.x);
			rightDepth = p0.depth + t2 * (p2.depth - p0.depth);
		}
		
		if (leftX > rightX) {
			std::swap(leftX, rightX);
			std::swap(leftDepth, rightDepth);
		}
		
		int startX = (int)round(leftX);
		int endX = (int)round(rightX);
		
		for (int x = startX; x <= endX; x++) {
			if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
				float depth = leftDepth;
				if (endX != startX) {
					float t = (x - leftX) / (rightX - leftX);
					depth = leftDepth + t * (rightDepth - leftDepth);
				}
				
				int pixelIndex = y * window.width + x;
				if (depth > depthBuffer[pixelIndex] && depth > 0) {
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
	for (int i = 0; i < WIDTH * HEIGHT; i++) {
		depthBuffer[i] = 0.0f;
	}
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
