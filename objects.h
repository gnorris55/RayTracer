#ifndef OBJECTS_H
#define OBJECTS_H


struct {
	float near;
	float left;
	float right;
	float bottom;
	float top;
} screen;


class Sphere {


	std::string name;
       	glm::vec3 position;
	glm::vec3 scale;
	glm::vec3 color;

	public:
	glm::mat4 A = glm::mat4(1.0);
	glm::mat4 invA = glm::mat4(1.0);
	float ka, kd, ks, kr, n;

	Sphere (std::string name, glm::vec3 position, glm::vec3 scale, glm::vec3 color, 
		float ambient, float diffuse, float specular, float kr, float n) {
		
		this->name = name;
		this->position = position;
		this->scale = scale;
		this->color = color;
		this->ka = ambient;
		this->kd = diffuse;
		this->ks = specular;
		this->kr = kr;
		this->n = n;

		createMatrix();

	}

	void createMatrix() {

		this->A = glm::translate(this->A, this->position);
		this->A = glm::scale(this->A, this->scale);
		this->invA = glm::inverse(this->A);
	}


	std::string getName() {
		return this->name;
	}

	glm::vec3 getColor() {
		return this->color;
	}



	
	void printName() {
		std::cout << this->name << "\n";
	}
};

class Light {

	std::string name;
	glm::vec3 position;
	glm::vec3 color;

	public:

	Light(std::string name, glm::vec3 position, glm::vec3 color) {
		this->name = name;
		this->position = position;
		this->color = color;
	}

	glm::vec3 getPosition() {
		return this->position;
	}

	glm::vec3 getColor() {
		return this->color;
	}

	std::string getName() {
		return this->name;
	}

};

#endif
