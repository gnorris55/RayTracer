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

//functions
void setImage();
void readFile(char *filename);
glm::vec3 rayTrace(glm::vec4 dirVec, glm::vec4 eye, float depth, int n);
glm::vec3 getPixelColor(Sphere currSphere, float intersection, glm::vec4 dirVec, glm::vec4 eye, int n);
int getShadowRay(glm::vec4 P, Light light);
float intersects(Sphere sphere, glm::vec4 dirVec, glm::vec4 eyePosition, float ignoreDepth);
glm::vec3 capRGB(glm::vec3 color);
float getCameraParams(char *fileLine);
void getScreenRes(char *fileLine);
void createSphere(char *ptr);
void createLight(char *ptr);

//global variables
std::vector<Sphere> spheres;
std::vector<Light> lights;
std::string outputFile;
glm::vec3 backColor;
glm::vec3 ambient;

const int MAX_DEPTH = 30;
const int NUM_REFLECTIONS = 3;
const glm::vec4 EYE_POSITION = glm::vec4(0.0, 0.0, 0.0, 1.0);

int screenRes[2];


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
			glm::vec3 color = rayTrace(dirVec, EYE_POSITION, 1, NUM_REFLECTIONS+1);
		
				
			int R = (int)(color.x * 255);
			int G = (int)(color.y * 255);
			int B = (int)(color.z * 255);

			//std::cout << R << ", " << G << ", " << B << "\n";
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
		glm::vec3 reflectedColor = rayTrace(reflection, P, 0.1, n-1);

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
		float th = intersects(spheres[i], PLight, P, 0.1);
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

	// nuumber of intersections
	float nIntersect = pow(B, 2) - A*C; 
	float t1, t2;

	if (nIntersect == 0) {	
		return -(B/A);

	} else if (nIntersect > 0) {
		// roots of the conanical sphere
		t1 = -(B/A) - sqrt(pow(B,2) - A*C)/A;	
		t2 = -(B/A) + sqrt(pow(B,2) - A*C)/A;	
		
		// if the root is less than desired length then ignore it
		if (t1 < ignoreDepth) t1 = -1;
		if (t2 < ignoreDepth) t2 = -1;
		
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


	// first six lines for screen parameters
	screen.near = getCameraParams(fgets(fileLine, 100, fptr));
	screen.left = getCameraParams(fgets(fileLine, 100, fptr));
	screen.right = getCameraParams(fgets(fileLine, 100, fptr));
	screen.bottom = getCameraParams(fgets(fileLine, 100, fptr));
	screen.top = getCameraParams(fgets(fileLine, 100, fptr));

	getScreenRes(fgets(fileLine, 100, fptr));
	
	char *ptr;
	while(fgets(fileLine, 100, fptr)) {
		
		ptr = std::strtok(fileLine, " ");

		if(strcmp(ptr, "SPHERE") == 0) {
			createSphere(ptr);
		} else if(strcmp(ptr, "LIGHT") == 0) {
			createLight(ptr);
		} else {
			break;
		}	
	}
	
	// background color	
	ptr = std::strtok(NULL, " ");
	backColor.x = atof(ptr);		
	ptr = std::strtok(NULL, " ");
	backColor.y = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	backColor.z = atof(ptr);	
	
	// ambient color
	fgets(fileLine, 100, fptr);
	
	ptr = std::strtok(fileLine, " ");
	ptr = std::strtok(NULL, " ");
	ambient.x = atof(ptr);		
	ptr = std::strtok(NULL, " ");
	ambient.y = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	ambient.z = atof(ptr);	

	// output file name
	fgets(fileLine, 100, fptr);
	
	ptr = std::strtok(fileLine, " ");
	ptr = std::strtok(NULL, " ");
	outputFile = ptr;		

	fclose(fptr);
}

// gets the values of screen
float getCameraParams(char *fileLine) {
	
	char *ptr;
	ptr = std::strtok(fileLine, " ");
	ptr = std::strtok(NULL, " ");
	return atof(ptr);

}

void getScreenRes(char *fileLine) {

	char *ptr;
	ptr = std::strtok(fileLine, " ");
	
	ptr = std::strtok(NULL, " ");
	screenRes[0] = atoi(ptr);	
	
	ptr = std::strtok(NULL, " ");
	screenRes[1] = atoi(ptr);	
}

void getVec3(char *fileLine) {

	char *ptr;
	ptr = std::strtok(fileLine, " ");
	
	ptr = std::strtok(NULL, " ");
	screenRes[0] = atoi(ptr);	
	
	ptr = std::strtok(NULL, " ");
	screenRes[1] = atoi(ptr);	
	
	ptr = std::strtok(NULL, " ");
	screenRes[1] = atoi(ptr);	
}

void createSphere(char *ptr) {

	std::string name;
        glm::vec3 position;
        glm::vec3 scale;
        glm::vec3 color;
        float ambient, diffuse, specular, Kr, n;
	
	//name
	ptr = std::strtok(NULL, " ");
	name = ptr;

	//position
	ptr = std::strtok(NULL, " ");
	position.x = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	position.y = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	position.z = atof(ptr);	

	//scale
	ptr = std::strtok(NULL, " ");
	scale.x = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	scale.y = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	scale.z = atof(ptr);	

	//color
	ptr = std::strtok(NULL, " ");
	color.x = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	color.y = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	color.z = atof(ptr);	
       
	// other	
	ptr = std::strtok(NULL, " ");
	ambient  = atof(ptr);
	ptr = std::strtok(NULL, " ");
	diffuse = atof(ptr);
	ptr = std::strtok(NULL, " ");
	specular = atof(ptr);
	ptr = std::strtok(NULL, " ");
	Kr = atof(ptr);
	ptr = std::strtok(NULL, " ");
	n = atof(ptr);
	
	Sphere tempSphere = Sphere(name, position, scale, color, ambient,  diffuse,  specular, Kr, n);

	spheres.push_back(tempSphere);

	
}

void createLight(char *ptr) {
	
	char *name;
        glm::vec3 position;
        glm::vec3 color;
	
	//name
	ptr = std::strtok(NULL, " ");
	name = ptr;

	//position
	ptr = std::strtok(NULL, " ");
	position.x = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	position.y = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	position.z = atof(ptr);	

	//color
	ptr = std::strtok(NULL, " ");
	color.x = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	color.y = atof(ptr);	
	ptr = std::strtok(NULL, " ");
	color.z = atof(ptr);	
       
	
	Light tempLight = Light(name, position, color);

	lights.push_back(tempLight);

}
