#define GLM_FORCE_SWIZZLE
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdio>
#include <cmath>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtx/string_cast.hpp"
#include "invert.h"
#include "objects.h"
#include "ppm.h"

//functionn for raytracing
void setImage();
glm::vec3 rayTrace(glm::vec4 dirVec, glm::vec4 eye, float depth, int n);
glm::vec3 getPixelColor(Sphere currSphere, float intersection, glm::vec4 dirVec, glm::vec4 eye, int n);
glm::vec3 capRGB(glm::vec3 color);
int getShadowRay(glm::vec4 P, Light light);
float intersects(Sphere sphere, glm::vec4 dirVec, glm::vec4 eyePosition, float ignoreDepth);


//functions for parsing the file
void readFile(char *filename);
void seperateFileLine(char *fileLine, std::vector<std::string> &stringList);
void allocateFileLineData(std::vector<std::string> stringList);
void createSphere(std::vector<std::string> attrList);
void createLight(std::vector<std::string> attrList);
void assignOtherParams(std::vector<std::string> attrList);

//global variables
std::vector<Sphere> spheres;
std::vector<Light> lights;
std::string outputFile;
glm::vec3 backColor;
glm::vec3 ambient;
int screenRes[2];

const int MAX_DEPTH = 30;
const int NUM_REFLECTIONS = 3;
const glm::vec4 EYE_POSITION = glm::vec4(0.0, 0.0, 0.0, 1.0);


int main(int argc, char *argv[]) {
        
	readFile(argv[1]);
	setImage();
	return 0;
}

void setImage() {
	
	int nRows = screenRes[0];
	int nCols = screenRes[1];

	int index = 0; //for pixel color placement
	unsigned char *pixels = new unsigned char [nRows*nCols*3];
	float planeHeight = 2 * screen.top * tan(glm::radians(55.0 / 2.0));
	float planeWidth = 2 * screen.right * tan(glm::radians(55.0 / 2.0));

	// for each pixel get direction vector from camera to pixel
	for (int r = nRows - 1; r > 0; r--) {
		for (int c = 0; c < nCols; c++) {
			
			// calculating coords for dirVec
			float uc = planeWidth*(2*c/(float)nCols - 1);
			float vr = planeHeight*(2*r/(float)nRows - 1);
			float N = screen.near;

			glm::vec4 dirVec = glm::vec4(uc, vr, -N, 0);

			dirVec = glm::normalize(dirVec);
			glm::vec3 color = rayTrace(dirVec, EYE_POSITION, 1.0, NUM_REFLECTIONS+1);
		
				
			int R = (int)(color.x * 255);
			int G = (int)(color.y * 255);
			int B = (int)(color.z * 255);

			pixels[index++] = char(R);
			pixels[index++] = char(G);
			pixels[index++] = char(B);

		}
	}

	char name[outputFile.length() + 1];
	strcpy(name, outputFile.c_str());
	save_imageP3(nCols, nRows, name, pixels);
}

glm::vec3 rayTrace(glm::vec4 dirVec, glm::vec4 eye, float depth, int n) {
	
	std::vector<float> rootValues;
	std::vector<int> sphereNum;
		
	//determining if the ray intersects with sphere(s)
	for(int i = 0; i < spheres.size(); i++) {
		float th = intersects(spheres[i], dirVec, eye, depth);
		
		if (th != 0) {
			rootValues.push_back(th);
			sphereNum.push_back(i);	
		} 	
	}

	// finding the closest intersection among the spheres
	int min = 0;
	for (int i = 0; i < rootValues.size(); i++) {
		if (rootValues[min] > rootValues[i]) {
			min = i;
		}
			
	}


	//if there is an intersection(s) calculate reflection and shadow rays
	if (rootValues.size() > 0) {	
		glm::vec3 pixelColor = getPixelColor(spheres[sphereNum[min]], rootValues[min], dirVec, eye, n);
		return capRGB(pixelColor);
	} else if (n == 4) {	// if n equals three it means it is the first ray. if the ray does not intersect 
				// anything we return the background color
		return backColor;
	}

	return glm::vec3(0.0, 0.0, 0.0);
}


glm::vec3 getPixelColor(Sphere currSphere, float intersection, glm::vec4 dirVec, glm::vec4 eye, int n) {

	glm::vec3 color = glm::vec3(0.0, 0.0, 0.0);

	// get point on original sphere
		
	glm::vec3 cInverse = (currSphere.invA * dirVec).xyz();
	glm::vec3 SInverse = (currSphere.invA * eye).xyz();
	
	// Point
	glm::vec4 P = eye + intersection*dirVec;
	
	//calculate normal of unit sphere then inverse transpose it for the original sphere
	glm::vec4 N = glm::vec4(SInverse + intersection*cInverse, 0);
	N = glm::vec4((glm::transpose(currSphere.invA)*N).xyz(), 0);
	N = glm::normalize(N);
	

	// calculate shadow rays
	glm::vec3 lightImpact = currSphere.ka*ambient*currSphere.getColor();
	for (int i = 0; i < lights.size(); i++) {
		if(getShadowRay(P, lights[i]) == 1) {
	
			//TODO: clean up		
			//helper vectors for specular and diffuse
			glm::vec4 L = glm::normalize(glm::vec4(lights[i].getPosition(), 1) - P);
			glm::vec4 V = glm::normalize(eye - P);
			glm::vec4 R = 2*(glm::dot(N,L))*N - L;
		
			// calculating diffuse color
			float diffuseAngle = glm::dot(N, L);
			if (diffuseAngle < 0) diffuseAngle = 0;	
			glm::vec3 lightDiffuse = currSphere.kd * lights[i].getColor() *
						 diffuseAngle * currSphere.getColor();

			// calculating specular color
			float specularAngle = glm::dot(R, V);
			if (specularAngle < 0) specularAngle = 0;
			specularAngle = pow(specularAngle, currSphere.n);		
			glm::vec3 lightSpecular = currSphere.ks * lights[i].getColor() *
						  specularAngle;

	
			// adding light effects for each light
			lightImpact += lightDiffuse + lightSpecular;	
		}
	} 
	color += lightImpact;
	

	// if we have relected 3 less than 3 times continue to recursively do so	
	if (n != 0) {

		glm::vec4 reflection = -2*(glm::dot(N, dirVec))*N + dirVec;
		glm::vec3 reflectedColor = rayTrace(reflection, P, 0.0001, n-1);

		color += reflectedColor*currSphere.kr;
	}

	return color;
}

int getShadowRay(glm::vec4 P, Light light) {
	
	// get vector from point to light position
	glm::vec4 PLight = glm::vec4(light.getPosition(), 1) - P;

	// testing all spheres to see if they intersect with the light ray
	for (int i = 0; i < spheres.size(); i++) {
		// if the scalar is not equal to zero, there is an object in-between the point and light
		float th = intersects(spheres[i], PLight, P, 0.0001);
		if (th != 0.0) {
			return 0.0;
		} 
	}

	// if there is no intersection, return the lights color
	return 1.0;
}

float intersects(Sphere sphere, glm::vec4 dirVec, glm::vec4 eyePosition, float ignoreDepth) {
	
	//get the inverse of each sphere, apply it to the ray and test intersection on the conanical sphere
	glm::vec3 cInverse = (sphere.invA * dirVec).xyz();
	glm::vec3 SInverse = (sphere.invA * eyePosition).xyz();
	

	// assigning the variables for the quadratic equation
	float A = pow(glm::length(cInverse), 2);
	float B = glm::dot(cInverse, SInverse);
	float C = pow(glm::length(SInverse), 2) - 1;

	float nIntersect = pow(B, 2) - A*C; 
	float t1, t2;

	if (nIntersect == 0) {	
		return -(B/A);

	} else if (nIntersect > 0) {
		// roots of the conanical sphere
		t1 = -(B/A) - sqrt(pow(B,2) - A*C)/A;	
		t2 = -(B/A) + sqrt(pow(B,2) - A*C)/A;	
			
		if (t1 <= ignoreDepth) t1 = -1;
		if (t2 <= ignoreDepth) t2 = -1;
		
		// cannot have two negative roots because it scales behind the 
		// camera or shadow ray, which is undesirable 
		if (t1 > 0 || t2 > 0) {
			// if both roots are positive return the smallest one
			if (t1 > 0 && t2 > 0) {
				if (t1 < t2)  	return t1;
				else  		return t2;
			} else {	// if both roots are not positive, returnt the positive one

				if (t1 > 0)	return t1;
				else 		return t2;
			}	
		}
	}
	
	// return 0 if there are no intersections
	return 0.0;
}

glm::vec3 capRGB(glm::vec3 color) {
	
	if (color.x > 1) 	color.x = 1;
	if (color.y > 1) 	color.y = 1;
	if (color.z > 1) 	color.z = 1;

	return color;

}


void readFile(char *filename) {

	FILE *fptr;


	fptr = fopen(filename, "r");
	char fileLine[100];
	

	std::vector<std::string> stringList;

	while(fgets(fileLine, 100, fptr)) {
		seperateFileLine(fileLine, stringList);
		if (stringList.size() > 0) {
			allocateFileLineData(stringList);
		}
		stringList.clear();
	}

}

void allocateFileLineData(std::vector<std::string> stringList) {
	
	if (stringList[0].compare("SPHERE") == 0) {
		createSphere(stringList);
	} else if (stringList[0].compare("LIGHT") == 0){
		createLight(stringList); 
	} else {
		assignOtherParams(stringList);
	}

}


void seperateFileLine(char *fileLine, std::vector<std::string> &stringList) {
	
	std::string tempString = ""; 
	int inString = 0;
	int i = 0;

	do {
		
		if (fileLine[i] == 9 || fileLine[i] == 32 || fileLine[i] == 10) {
			if (inString == 1) {
				stringList.push_back(tempString);
				tempString = "";
				inString = 0;
			}
		} else { 
			if (inString == 0) {
				inString = 1;
			}
			tempString = tempString + fileLine[i];
		}

		
	} while (fileLine[i++] != 10);
}

void assignOtherParams(std::vector<std::string> attrList) {
	
	if (attrList[0].compare("OUTPUT") == 0) outputFile = attrList[1]; 
	if (attrList[0].compare("LEFT") == 0)  screen.left = std::stoi(attrList[1]);
	if (attrList[0].compare("RIGHT") == 0) screen.right = std::stoi(attrList[1]); 
	if (attrList[0].compare("BOTTOM") == 0) screen.bottom = std::stoi(attrList[1]);
	if (attrList[0].compare("TOP") == 0) screen.top = std::stoi(attrList[1]);
	if (attrList[0].compare("NEAR") == 0) screen.near = std::stoi(attrList[1]);
	if (attrList[0].compare("AMBIENT") == 0) ambient = glm::vec3(std::stof(attrList[1]), std::stof(attrList[2]), std::stof(attrList[3]));   
	if (attrList[0].compare("BACK") == 0) backColor = glm::vec3(std::stof(attrList[1]), std::stof(attrList[2]), std::stof(attrList[3]));
	if (attrList[0].compare("RES") == 0) {
		screenRes[0] = std::stoi(attrList[1]);		
		screenRes[1] = std::stoi(attrList[2]);		
	}	
}


void createSphere(std::vector<std::string> attrList) {

        float ambient, diffuse, specular, Kr, n;
	
	std::string name = attrList[1];
	glm::vec3 position = glm::vec3(std::stof(attrList[2]), std::stof(attrList[3]), std::stof(attrList[4]));
	glm::vec3 scale = glm::vec3(std::stof(attrList[5]), std::stof(attrList[6]), std::stof(attrList[7]));
	glm::vec3 color = glm::vec3(std::stof(attrList[8]), std::stof(attrList[9]), std::stof(attrList[10]));
       
	// other	
	ambient  = std::stof(attrList[11]);
	diffuse = std::stof(attrList[12]);
	specular = std::stof(attrList[13]);
	Kr = std::stof(attrList[14]);
	n = std::stof(attrList[15]);
	
	Sphere tempSphere = Sphere(name, position, scale, color, ambient,  diffuse,  specular, Kr, n);
	spheres.push_back(tempSphere);

	
}

void createLight(std::vector<std::string> attrList) {
	
	std::string name = attrList[1];
	glm::vec3 position = glm::vec3(std::stof(attrList[2]), std::stof(attrList[3]), std::stof(attrList[4]));
	glm::vec3 color = glm::vec3(std::stof(attrList[5]), std::stof(attrList[6]), std::stof(attrList[7]));

	Light tempLight = Light(name, position, color);
	lights.push_back(tempLight);

}

