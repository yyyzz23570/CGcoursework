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
#include <limits>
#include <glm/glm.hpp>

#define WIDTH 640
#define HEIGHT 480

// 全局变量
std::vector<ModelTriangle> modelTriangles;
std::vector<float> depthBuffer;
glm::vec3 cameraPosition(0.0f, 0.0f, -5.0f);
float focalLength = 200.0f;

enum RenderMode {
	WIREFRAME = 0,
	RASTERIZED = 1
};

RenderMode renderMode = RASTERIZED;

// ==================== 工具函数 ====================

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

// ==================== 文件加载函数 ====================

std::map<std::string, Colour> loadMTLFile(const std::string &filename) {
	std::map<std::string, Colour> materials;
	std::ifstream file(filename);
	std::string line;
	std::string currentMaterial;
	
	if (!file.is_open()) {
		std::cerr << "错误: 无法打开MTL文件: " << filename << std::endl;
		return materials;
	}
	
	std::cout << "正在加载材质文件: " << filename << std::endl;
	
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
			std::cout << "  加载材质: " << currentMaterial 
			          << " (RGB: " << static_cast<int>(r * 255) << ", " 
			          << static_cast<int>(g * 255) << ", " 
			          << static_cast<int>(b * 255) << ")" << std::endl;
		}
	}
	
	file.close();
	std::cout << "总共加载材质数: " << materials.size() << std::endl;
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
		std::cerr << "错误: 无法打开OBJ文件: " << filename << std::endl;
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
	std::cout << "从文件加载了 " << triangles.size() << " 个三角形: " << filename << std::endl;
	return triangles;
}

// ==================== 投影函数 ====================

CanvasPoint projectVertexOntoCanvasPoint(glm::vec3 vertexPosition) {
	glm::vec3 vertexInCameraSpace = vertexPosition - cameraPosition;
	
	float x = vertexInCameraSpace.x;
	float y = vertexInCameraSpace.y;
	float z = vertexInCameraSpace.z;
	
	if (z <= 0) {
		z = 0.1f;
	}
	
	// 透视投影：u = f * (x/z) + width/2, v = f * (y/z) + height/2
	float u = focalLength * (x / z) + WIDTH * 0.5f;
	float v = focalLength * (y / z) + HEIGHT * 0.5f;
	
	// 翻转Y轴（屏幕坐标系中Y向下）
	v = HEIGHT - v;
	
	// 计算逆深度（1/z）用于深度缓冲
	float depth = 1.0f / z;
	
	return CanvasPoint(u, v, depth);
}

// ==================== 绘制函数 ====================

void drawStrokedTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
	drawLine(window, triangle[0], triangle[1], colour);
	drawLine(window, triangle[1], triangle[2], colour);
	drawLine(window, triangle[2], triangle[0], colour);
}

void drawFilledTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
	CanvasPoint p0 = triangle[0];
	CanvasPoint p1 = triangle[1];
	CanvasPoint p2 = triangle[2];
	
	glm::vec2 v0(p0.x, p0.y);
	glm::vec2 v1(p1.x, p1.y);
	glm::vec2 v2(p2.x, p2.y);
	
	float minX = std::min({p0.x, p1.x, p2.x});
	float maxX = std::max({p0.x, p1.x, p2.x});
	float minY = std::min({p0.y, p1.y, p2.y});
	float maxY = std::max({p0.y, p1.y, p2.y});
	
	uint32_t fillColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
	
	int startX = (int)std::max(0.0f, std::floor(minX));
	int endX = (int)std::min((float)window.width - 1, std::ceil(maxX));
	int startY = (int)std::max(0.0f, std::floor(minY));
	int endY = (int)std::min((float)window.height - 1, std::ceil(maxY));
	
	for (int y = startY; y <= endY; y++) {
		for (int x = startX; x <= endX; x++) {
			glm::vec2 r((float)x, (float)y);
			glm::vec3 barycentric = convertToBarycentricCoordinates(v0, v1, v2, r);
			
			// 检查点是否在三角形内部
			if (barycentric.x >= 0 && barycentric.y >= 0 && barycentric.z >= 0) {
				// 使用重心坐标进行深度插值
				// barycentric.x 对应 v0/p0
				// barycentric.y 对应 v1/p1
				// barycentric.z 对应 v2/p2
				float invDepth = barycentric.x * p0.depth + barycentric.y * p1.depth + barycentric.z * p2.depth;
				
				int pixelIndex = y * window.width + x;
				if (pixelIndex >= 0 && pixelIndex < static_cast<int>(depthBuffer.size())) {
					// 比较深度值：如果新的深度更近（更大的1/z值），则更新像素
					if (invDepth > depthBuffer[pixelIndex]) {
						depthBuffer[pixelIndex] = invDepth;
						window.setPixelColour(x, y, fillColour);
					}
				}
			}
		}
	}
}

void drawWireframe(DrawingWindow &window) {
	Colour white(255, 255, 255);
	
	for (const auto &modelTriangle : modelTriangles) {
		CanvasPoint p0 = projectVertexOntoCanvasPoint(modelTriangle.vertices[0]);
		CanvasPoint p1 = projectVertexOntoCanvasPoint(modelTriangle.vertices[1]);
		CanvasPoint p2 = projectVertexOntoCanvasPoint(modelTriangle.vertices[2]);
		
		CanvasTriangle canvasTriangle(p0, p1, p2);
		drawStrokedTriangle(window, canvasTriangle, white);
	}
}

void drawRasterized(DrawingWindow &window) {
	// 初始化深度缓冲
	depthBuffer.clear();
	depthBuffer.resize(window.width * window.height, 0.0f);
	
	for (const auto &modelTriangle : modelTriangles) {
		CanvasPoint p0 = projectVertexOntoCanvasPoint(modelTriangle.vertices[0]);
		CanvasPoint p1 = projectVertexOntoCanvasPoint(modelTriangle.vertices[1]);
		CanvasPoint p2 = projectVertexOntoCanvasPoint(modelTriangle.vertices[2]);
		
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

// ==================== 事件处理 ====================

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_1) {
			renderMode = WIREFRAME;
			std::cout << "切换到线框渲染模式" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_2) {
			renderMode = RASTERIZED;
			std::cout << "切换到光栅化渲染模式" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_LEFT) {
			cameraPosition.x -= 0.2f;
			std::cout << "相机左移，位置: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_RIGHT) {
			cameraPosition.x += 0.2f;
			std::cout << "相机右移，位置: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
		else if (event.key.keysym.sym == SDLK_UP) {
			cameraPosition.y += 0.2f;
			std::cout << "相机上移，位置: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
        else if (event.key.keysym.sym == SDLK_DOWN) {
			cameraPosition.y -= 0.2f;
			std::cout << "相机下移，位置: (" << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << ")" << std::endl;
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
		std::cout << "已保存帧为 output.ppm 和 output.bmp" << std::endl;
	}
}

// ==================== 主函数 ====================

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	
	// 加载OBJ文件和MTL文件
	float scaleFactor = 0.35f;
	modelTriangles = loadOBJFile("cornell-box.obj", "cornell-box.mtl", scaleFactor);
	
	std::cout << "\n========== 加载信息 ==========" << std::endl;
	std::cout << "总共加载了 " << modelTriangles.size() << " 个三角形" << std::endl;
	std::cout << "按键1: 线框渲染" << std::endl;
	std::cout << "按键2: 光栅化渲染" << std::endl;
	std::cout << "方向键: 移动相机" << std::endl;
	std::cout << "鼠标点击: 保存帧" << std::endl;
	std::cout << "============================\n" << std::endl;
	
	// 主循环
	while (true) {
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		window.renderFrame();
	}
	
	return 0;
}