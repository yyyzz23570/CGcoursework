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

#define WIDTH 640
#define HEIGHT 480

std::vector<ModelTriangle> modelTriangles;

enum RenderMode {
	WIREFRAME = 0,
	RASTERIZED = 1
};

RenderMode renderMode = RASTERIZED;

glm::vec3 cameraPosition(0.0f, 0.0f, -3.0f);
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
	
	float numberOfSteps = std::max(std::abs(xDiff), std::abs(yDiff));
	
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
		std::cerr << "Error: Cannot open MTL file: " << filename << std::endl;
		return materials;
	}
	
	std::cout << "Loading materials from: " << filename << std::endl;
	
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
			std::cout << "  Loaded material: " << currentMaterial 
			          << " (RGB: " << static_cast<int>(r * 255) << ", " 
			          << static_cast<int>(g * 255) << ", " 
			          << static_cast<int>(b * 255) << ")" << std::endl;
		}
	}
	
	file.close();
	std::cout << "Total materials loaded: " << materials.size() << std::endl;
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
		std::cerr << "Error: Cannot open OBJ file: " << filename << std::endl;
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
			} else {
				std::cerr << "Warning: Material '" << materialName << "' not found, using default color" << std::endl;
			}
		} else if (token == "f") {
			std::vector<int> vertexIndices;
			
			std::string vertexStr;
			std::string remainingLine = line.substr(1);
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
	std::cout << "Loaded " << triangles.size() << " triangles from " << filename << std::endl;
	return triangles;
}

CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 cameraPosition, float focalLength, glm::vec3 vertexPosition) {
	glm::vec3 vertexInCameraSpace = vertexPosition - cameraPosition;
	
	float x = vertexInCameraSpace.x;
	float y = vertexInCameraSpace.y;
	float z = vertexInCameraSpace.z;
	
	if (z <= 0) {
		z = 0.1f;
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
		
		if (y <= p1.y) {
			if (p1.y != p0.y) {
				t = (y - p0.y) / (p1.y - p0.y);
				leftX = p0.x + t * (p1.x - p0.x);
			} else {
				leftX = p0.x;
			}
			
			if (p2.y != p0.y) {
				t = (y - p0.y) / (p2.y - p0.y);
				rightX = p0.x + t * (p2.x - p0.x);
			} else {
				rightX = p0.x;
			}
		} else {
			if (p2.y != p1.y) {
				t = (y - p1.y) / (p2.y - p1.y);
				leftX = p1.x + t * (p2.x - p1.x);
			} else {
				leftX = p1.x;
			}
			
			if (p2.y != p0.y) {
				t = (y - p0.y) / (p2.y - p0.y);
				rightX = p0.x + t * (p2.x - p0.x);
			} else {
				rightX = p0.x;
			}
		}
		
		if (leftX > rightX) std::swap(leftX, rightX);
		
		int startX = (int)round(leftX);
		int endX = (int)round(rightX);
		
		for (int x = startX; x <= endX; x++) {
			if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
				window.setPixelColour(x, y, fillColour);
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

void draw(DrawingWindow &window) {
	window.clearPixels();
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
			std::cout << "Switched to Wireframe Render" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_2) {
			renderMode = RASTERIZED;
			std::cout << "Switched to Rasterised Render" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_LEFT) {
			cameraPosition.x -= 0.1f;
			std::cout << "Camera moved left. Position: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			cameraPosition.x += 0.1f;
			std::cout << "Camera moved right. Position: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_UP) {
			cameraPosition.y += 0.1f;
			std::cout << "Camera moved up. Position: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_DOWN) {
			cameraPosition.y -= 0.1f;
			std::cout << "Camera moved down. Position: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_w) {
			cameraPosition.z -= 0.1f;
			std::cout << "Camera moved forward. Position: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_s) {
			cameraPosition.z += 0.1f;
			std::cout << "Camera moved backward. Position: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
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
	
	std::cout << "\n=== Verifying loaded triangles ===" << std::endl;
	size_t numToPrint = std::min(static_cast<size_t>(5), modelTriangles.size());
	for (size_t i = 0; i < numToPrint; i++) {
		std::cout << "Triangle " << i << ":" << std::endl;
		std::cout << modelTriangles[i] << std::endl;
		std::cout << "  Colour: " << modelTriangles[i].colour.name 
		          << " (RGB: " << modelTriangles[i].colour.red << ", " 
		          << modelTriangles[i].colour.green << ", " 
		          << modelTriangles[i].colour.blue << ")" << std::endl;
	}
	if (modelTriangles.size() > numToPrint) {
		std::cout << "... and " << (modelTriangles.size() - numToPrint) << " more triangles" << std::endl;
	}
	std::cout << "===================================\n" << std::endl;
	
	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		window.renderFrame();
	}
}
