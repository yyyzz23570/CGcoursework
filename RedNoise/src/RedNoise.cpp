#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <ModelTriangle.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <map>
#include <sstream>
#include <glm/glm.hpp>

#define WIDTH 640
#define HEIGHT 480

// Global variable to store 3D model triangles
std::vector<ModelTriangle> modelTriangles;

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
	// 计算x和y的差值
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	
	// 选择较大的差值作为步数，确保线条连续
	float numberOfSteps = std::max(std::abs(xDiff), std::abs(yDiff));
	
	// 插值x和y坐标
	std::vector<float> xValues = interpolateSingleFloats(from.x, to.x, numberOfSteps);
	std::vector<float> yValues = interpolateSingleFloats(from.y, to.y, numberOfSteps);
	
	// 将颜色转换为uint32_t格式 (ARGB)
	uint32_t pixelColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
	
	// 绘制线条上的每个点
	for (int i = 0; i < numberOfSteps; i++) {
		int x = round(xValues[i]);
		int y = round(yValues[i]);
		// 边界检查
		if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
			window.setPixelColour(x, y, pixelColour);
		}
	}
}

// Parse MTL file to get material colors
std::map<std::string, Colour> loadMTLFile(const std::string &filename) {
	std::map<std::string, Colour> materials; // Hashmap using material name as key
	std::ifstream file(filename);
	std::string line;
	std::string currentMaterial;
	
	if (!file.is_open()) {
		std::cerr << "Error: Cannot open MTL file: " << filename << std::endl;
		return materials;
	}
	
	std::cout << "Loading materials from: " << filename << std::endl;
	
	while (std::getline(file, line)) {
		// Skip empty lines
		if (line.empty()) continue;
		
		std::istringstream iss(line);
		std::string token;
		iss >> token;
		
		if (token == "newmtl") {
			// New material definition
			iss >> currentMaterial;
		} else if (token == "Kd" && !currentMaterial.empty()) {
			// Diffuse color (Kd)
			float r, g, b;
			iss >> r >> g >> b;
			// Convert from 0-1 range to 0-255 range
			// Use the Colour constructor that sets the name attribute
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

// Parse OBJ file and return vector of ModelTriangles
std::vector<ModelTriangle> loadOBJFile(const std::string &filename, const std::string &mtlPath, float scaleFactor = 1.0f) {
	std::vector<ModelTriangle> triangles;
	std::vector<glm::vec3> vertices;
	std::map<std::string, Colour> materials = loadMTLFile(mtlPath);
	Colour currentColour(255, 255, 255); // Default white
	
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
			// Vertex - apply scaling factor when reading
			float x, y, z;
			iss >> x >> y >> z;
			// Scale the vertex coordinates
			vertices.push_back(glm::vec3(x * scaleFactor, y * scaleFactor, z * scaleFactor));
		} else if (token == "usemtl") {
			// Material reference - look up the color from our palette
			std::string materialName;
			iss >> materialName;
			// Use hashmap lookup for efficient color retrieval
			if (materials.find(materialName) != materials.end()) {
				currentColour = materials[materialName];
				// The ModelTriangle will use this colour when created
			} else {
				std::cerr << "Warning: Material '" << materialName << "' not found, using default color" << std::endl;
			}
		} else if (token == "f") {
			// Face - OBJ indices are 1-based
			std::vector<int> vertexIndices;
			
			std::string vertexStr;
			// Read the rest of the line
			std::string remainingLine = line.substr(1); // Remove 'f' character
			std::istringstream faceIss(remainingLine);
			
			while (faceIss >> vertexStr) {
				// Handle format like "1", "1/", "1/2", "1/2/3", "1//3"
				// We only need the vertex index (first number)
				if (vertexStr.empty()) continue;
				
				size_t slashPos = vertexStr.find('/');
				if (slashPos != std::string::npos) {
					vertexStr = vertexStr.substr(0, slashPos);
				}
				
				if (!vertexStr.empty()) {
					try {
						int idx = std::stoi(vertexStr) - 1; // Convert to 0-based
						if (idx >= 0 && idx < static_cast<int>(vertices.size())) {
							vertexIndices.push_back(idx);
						}
					} catch (const std::exception &e) {
						// Skip invalid indices
						continue;
					}
				}
			}
			
			// Create triangles from face (triangulate if needed)
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

// Simple 3D to 2D projection (perspective projection)
CanvasPoint projectPoint(glm::vec3 point, float focalLength = 2.0f) {
	// Simple perspective projection
	// Move camera backward to see more of the scene
	float cameraZ = -3.5f; // Camera positioned behind the scene
	float offsetZ = point.z - cameraZ;
	
	if (offsetZ <= 0) {
		// Behind camera, don't render (or handle specially)
		offsetZ = 0.1f; // Small epsilon to avoid division by zero
	}
	
	float scale = focalLength / offsetZ;
	float x = point.x * scale;
	float y = point.y * scale;
	
	// Scale to fit the screen
	// After applying 0.35 scale factor, vertices are roughly in -1.0 to 1.0 range
	// Adjust scale factor accordingly
	float scaleFactor = 160.0f;
	
	// Convert to screen coordinates (center at origin, then offset)
	float screenX = x * scaleFactor + WIDTH * 0.5f;
	float screenY = -y * scaleFactor + HEIGHT * 0.5f; // Flip Y axis and center
	
	return CanvasPoint(screenX, screenY, point.z);
}

// Convert ModelTriangle to CanvasTriangle for rendering
CanvasTriangle modelToCanvas(const ModelTriangle &triangle, float focalLength = 2.0f) {
	CanvasPoint p0 = projectPoint(triangle.vertices[0], focalLength);
	CanvasPoint p1 = projectPoint(triangle.vertices[1], focalLength);
	CanvasPoint p2 = projectPoint(triangle.vertices[2], focalLength);
	
	return CanvasTriangle(p0, p1, p2);
}

void drawStrokedTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
	// 依次绘制三条边
	drawLine(window, triangle[0], triangle[1], colour);
	drawLine(window, triangle[1], triangle[2], colour);
	drawLine(window, triangle[2], triangle[0], colour);
}

void drawFilledTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const Colour &colour) {
	// 获取三个顶点
	CanvasPoint p0 = triangle[0];
	CanvasPoint p1 = triangle[1];
	CanvasPoint p2 = triangle[2];
	
	// 按y坐标排序顶点（从上到下）
	if (p0.y > p1.y) std::swap(p0, p1);
	if (p0.y > p2.y) std::swap(p0, p2);
	if (p1.y > p2.y) std::swap(p1, p2);
	
	// 计算填充颜色
	uint32_t fillColour = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
	
	// 扫描线算法：从上到下逐行填充
	int startY = (int)round(p0.y);
	int endY = (int)round(p2.y);
	
	for (int y = startY; y <= endY; y++) {
		float t;
		float leftX, rightX;
		
		if (y <= p1.y) {
			// 上半部分：p0到p1，p0到p2
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
			// 下半部分：p1到p2，p0到p2
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
		
		// 确保leftX <= rightX
		if (leftX > rightX) std::swap(leftX, rightX);
		
		// 绘制水平线段
		int startX = (int)round(leftX);
		int endX = (int)round(rightX);
		
		for (int x = startX; x <= endX; x++) {
			if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
				window.setPixelColour(x, y, fillColour);
			}
		}
	}
	
	// 绘制白色描边
	Colour white(255, 255, 255);
	drawStrokedTriangle(window, triangle, white);
}

void drawTexturedTriangle(DrawingWindow &window, const CanvasTriangle &triangle, const TextureMap &texture) {
	// 获取三个顶点
	CanvasPoint p0 = triangle[0];
	CanvasPoint p1 = triangle[1];
	CanvasPoint p2 = triangle[2];
	
	// 按y坐标排序顶点（从上到下）
	if (p0.y > p1.y) std::swap(p0, p1);
	if (p0.y > p2.y) std::swap(p0, p2);
	if (p1.y > p2.y) std::swap(p1, p2);
	
	// 扫描线算法：从上到下逐行填充
	int startY = (int)round(p0.y);
	int endY = (int)round(p2.y);
	
	for (int y = startY; y <= endY; y++) {
		float t;
		float leftX, rightX;
		float leftTexX, leftTexY, rightTexX, rightTexY;
		
		if (y <= p1.y) {
			// 上半部分：p0到p1，p0到p2
			if (p1.y != p0.y) {
				t = (y - p0.y) / (p1.y - p0.y);
				leftX = p0.x + t * (p1.x - p0.x);
				leftTexX = p0.texturePoint.x + t * (p1.texturePoint.x - p0.texturePoint.x);
				leftTexY = p0.texturePoint.y + t * (p1.texturePoint.y - p0.texturePoint.y);
			} else {
				leftX = p0.x;
				leftTexX = p0.texturePoint.x;
				leftTexY = p0.texturePoint.y;
			}
			
			if (p2.y != p0.y) {
				t = (y - p0.y) / (p2.y - p0.y);
				rightX = p0.x + t * (p2.x - p0.x);
				rightTexX = p0.texturePoint.x + t * (p2.texturePoint.x - p0.texturePoint.x);
				rightTexY = p0.texturePoint.y + t * (p2.texturePoint.y - p0.texturePoint.y);
			} else {
				rightX = p0.x;
				rightTexX = p0.texturePoint.x;
				rightTexY = p0.texturePoint.y;
			}
		} else {
			// 下半部分：p1到p2，p0到p2
			if (p2.y != p1.y) {
				t = (y - p1.y) / (p2.y - p1.y);
				leftX = p1.x + t * (p2.x - p1.x);
				leftTexX = p1.texturePoint.x + t * (p2.texturePoint.x - p1.texturePoint.x);
				leftTexY = p1.texturePoint.y + t * (p2.texturePoint.y - p1.texturePoint.y);
			} else {
				leftX = p1.x;
				leftTexX = p1.texturePoint.x;
				leftTexY = p1.texturePoint.y;
			}
			
			if (p2.y != p0.y) {
				t = (y - p0.y) / (p2.y - p0.y);
				rightX = p0.x + t * (p2.x - p0.x);
				rightTexX = p0.texturePoint.x + t * (p2.texturePoint.x - p0.texturePoint.x);
				rightTexY = p0.texturePoint.y + t * (p2.texturePoint.y - p0.texturePoint.y);
			} else {
				rightX = p0.x;
				rightTexX = p0.texturePoint.x;
				rightTexY = p0.texturePoint.y;
			}
		}
		
		// 确保leftX <= rightX
		if (leftX > rightX) {
			std::swap(leftX, rightX);
			std::swap(leftTexX, rightTexX);
			std::swap(leftTexY, rightTexY);
		}
		
		// 绘制水平线段
		int startX = (int)round(leftX);
		int endX = (int)round(rightX);
		
		for (int x = startX; x <= endX; x++) {
			if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
				// 计算当前像素的纹理坐标
				float texT = (endX != startX) ? (x - startX) / (float)(endX - startX) : 0.0f;
				float texX = leftTexX + texT * (rightTexX - leftTexX);
				float texY = leftTexY + texT * (rightTexY - leftTexY);
				
				// 将纹理坐标转换为像素索引
				int texPixelX = (int)round(texX * (texture.width - 1));
				int texPixelY = (int)round(texY * (texture.height - 1));
				
				// 边界检查
				texPixelX = std::max(0, std::min(texPixelX, (int)texture.width - 1));
				texPixelY = std::max(0, std::min(texPixelY, (int)texture.height - 1));
				
				// 获取纹理像素
				int pixelIndex = texPixelY * texture.width + texPixelX;
				uint32_t textureColour = texture.pixels[pixelIndex];
				
				window.setPixelColour(x, y, textureColour);
			}
		}
	}
	
	// 绘制白色描边
	Colour white(255, 255, 255);
	drawStrokedTriangle(window, triangle, white);
}

void draw(DrawingWindow &window) {
	// Clear the window (black background)
	window.clearPixels();
	
	// Render all 3D triangles
	for (const auto &modelTriangle : modelTriangles) {
		CanvasTriangle canvasTriangle = modelToCanvas(modelTriangle);
		// Draw filled triangle with the triangle's color
		drawFilledTriangle(window, canvasTriangle, modelTriangle.colour);
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window, const TextureMap &texture) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u) {
        	// 随机生成三角形三个顶点
        	CanvasPoint p0(rand() % window.width, rand() % window.height);
        	CanvasPoint p1(rand() % window.width, rand() % window.height);
        	CanvasPoint p2(rand() % window.width, rand() % window.height);
        	CanvasTriangle tri(p0, p1, p2);
        	// 随机颜色
        	Colour col(rand() % 256, rand() % 256, rand() % 256);
        	// 绘制描边三角形
        	drawStrokedTriangle(window, tri, col);
        }
        else if (event.key.keysym.sym == SDLK_f) {
        	// 随机生成三角形三个顶点
        	CanvasPoint p0(rand() % window.width, rand() % window.height);
        	CanvasPoint p1(rand() % window.width, rand() % window.height);
        	CanvasPoint p2(rand() % window.width, rand() % window.height);
        	CanvasTriangle tri(p0, p1, p2);
        	// 随机颜色
        	Colour col(rand() % 256, rand() % 256, rand() % 256);
        	// 绘制填充三角形（包含白色描边）
        	drawFilledTriangle(window, tri, col);
        }
        else if (event.key.keysym.sym == SDLK_t) {
        	// 创建带纹理坐标的三角形顶点
        	CanvasPoint p0(rand() % window.width, rand() % window.height);
        	p0.texturePoint = TexturePoint(0.0, 0.0);  // 左上角
        	
        	CanvasPoint p1(rand() % window.width, rand() % window.height);
        	p1.texturePoint = TexturePoint(1.0, 0.0);  // 右上角
        	
        	CanvasPoint p2(rand() % window.width, rand() % window.height);
        	p2.texturePoint = TexturePoint(0.5, 1.0);  // 底部中心
        	
        	CanvasTriangle tri(p0, p1, p2);
        	// 绘制带纹理的三角形
        	drawTexturedTriangle(window, tri, texture);
        }
        else if (event.key.keysym.sym == SDLK_v) {
        	// 视觉验证测试：使用指定的坐标
        	CanvasPoint p0(160, 10);
        	p0.texturePoint = TexturePoint(195.0/400.0, 5.0/400.0);  // 归一化纹理坐标
        	
        	CanvasPoint p1(300, 230);
        	p1.texturePoint = TexturePoint(395.0/400.0, 380.0/400.0);
        	
        	CanvasPoint p2(10, 150);
        	p2.texturePoint = TexturePoint(65.0/400.0, 330.0/400.0);
        	
        	CanvasTriangle tri(p0, p1, p2);
        	// 绘制带纹理的三角形
        	drawTexturedTriangle(window, tri, texture);
        	std::cout << "视觉验证测试三角形已绘制" << std::endl;
        }
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
    // 随机数种子
    srand(static_cast<unsigned int>(time(nullptr)));
    
    // 加载纹理
    TextureMap texture("texture.ppm");
    std::cout << "纹理加载完成: " << texture.width << "x" << texture.height << std::endl;
	
	// Load 3D model with scaling factor of 0.35 (as recommended for Cornell Box)
	float scaleFactor = 0.35f;
	modelTriangles = loadOBJFile("cornell-box.obj", "cornell-box.mtl", scaleFactor);
	
	// Print out the first few triangles to verify loading
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
	
	// 测试 interpolateSingleFloats 函数
	std::vector<float> result;
	result = interpolateSingleFloats(2.2, 8.5, 7);
	std::cout << "插值测试结果: ";
	for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	std::cout << std::endl;
	
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window, texture);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
