#ifndef MESH_H
#define MESH_H

#include <string>
#include <algorithm>
#include <vector>

#include <math.h>
#include "float.h"

#include <glm/vec3.hpp>

class Mesh{

protected:
    std::vector<glm::vec3> geom;
    std::vector<int> topo;
    std::vector<glm::vec3> normalsPoints;
    std::vector<int> normalIds;

public :

    Mesh(){}
    Mesh(std::vector<glm::vec3> listGeom, std::vector<int> listTopo): geom(listGeom), topo(listTopo)  {}

    ~Mesh(){}

    void merge(const Mesh &delta);

    std::vector<glm::vec3> getGeom() const;
    std::vector<int> getTopo() const;
    void setGeom(std::vector<glm::vec3> geom);
    void setTopo(std::vector<int> topo);
    size_t nbGeom() const;
    size_t nbTopo() const;

    void translation(const float x, const float y, const float z);
    void rescale(float scale);

    void loadFromOBJ(const glm::vec3 &center, const char* obj);

    const glm::vec3 getNormal(const glm::vec3& point) const;

    //float intersect(const Ray &ray);


protected:
};


#endif // MESH_H
