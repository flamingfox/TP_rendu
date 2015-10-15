#include "mesh.h"


/* -------------------------------------------- */
/* -------------------------------------------- */
/* -------------------------------------------- */
/*            Fonctions utilitaires             */
/* -------------------------------------------- */
/* -------------------------------------------- */
/* -------------------------------------------- */

void Mesh::merge(const Mesh &delta)
{
    if(&delta != this){
        int taille = geom.size();
        geom.reserve(taille+delta.nbGeom());

        for(size_t i=0; i< delta.geom.size(); i++){
            geom.push_back(delta.geom[i]);
        }

        topo.reserve(this->nbTopo()+delta.nbTopo());
        for(size_t i=0; i< delta.topo.size(); i++){
            topo.push_back(delta.topo[i] + taille );
        }
    }
}

/* -------------------------------------------- */
/* -------------------------------------------- */
/* -------------------------------------------- */
/*            Fonctions accessors               */
/* -------------------------------------------- */
/* -------------------------------------------- */
/* -------------------------------------------- */

std::vector<glm::vec3> Mesh::getGeom() const{
    return geom;
}

std::vector<int> Mesh::getTopo() const{
    return topo;
}

size_t Mesh::nbGeom() const
{
    return this->geom.size();
}

size_t Mesh::nbTopo() const
{
    return this->topo.size();
}

void Mesh::translation(const float x, const float y, const float z)
{
    for(glm::vec3& p: geom){
        p.x+=x;
        p.y+=y;
        p.z+=z;
    }
}

void Mesh::rescale(float scale)
{
    for(glm::vec3& p: geom){
        p *= scale;
    }
}

void Mesh::loadFromOBJ(const glm::vec3 &center, const char* obj){
    glm::vec3 minVal(1E100, 1E100, 1E100), maxVal(-1E100, -1E100, -1E100);
    FILE* f = fopen(obj, "r");
    while (!feof(f)) {
        char line[255];
        fgets(line, 255, f);
        if (line[0]=='v' && line[1]==' ') {
            glm::vec3 vec;
            sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);
            vec[2] = -vec[2];
            glm::vec3 p = vec*50.f + center;
            geom.push_back(p);
            maxVal[0] = std::max(maxVal[0], p[0]);
            maxVal[1] = std::max(maxVal[1], p[1]);
            maxVal[2] = std::max(maxVal[2], p[2]);
            minVal[0] = std::min(minVal[0], p[0]);
            minVal[1] = std::min(minVal[1], p[1]);
            minVal[2] = std::min(minVal[2], p[2]);
        }
        if (line[0]=='v' && line[1]=='n') {
            glm::vec3 vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);
            vec[2] = -vec[2];
            normalsPoints.push_back(vec);
        }
        if (line[0]=='f') {
            int i0, i1, i2;
            int j0,j1,j2;
            int k0,k1,k2;

            int n = 0;
            //sscanf(line, "f %u/%u/%u %u/%u/%u %u/%u/%u\n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2 );

            //count of '/' occuration in the line
            for (int i=0; i < sizeof(line); i++){
                if(line[i]=='/')
                    n++;
            }

            if(n==0){
                sscanf(line, "f %u %u %u\n", &i0, &i1, &i2);

                topo.push_back(i0-1);
                topo.push_back(i1-1);
                topo.push_back(i2-1);

            }
            else if(n==3){
                sscanf(line, "f %u/%u %u/%u %u/%u\n", &i0, &k0, &i1, &k1, &i2, &k2 );

                topo.push_back(i0-1);
                topo.push_back(i1-1);
                topo.push_back(i2-1);
                normalIds.push_back(k0-1);
                normalIds.push_back(k1-1);
                normalIds.push_back(k2-1);

            }
            else if(n==6){
                sscanf(line, "f %u/%u/%u %u/%u/%u %u/%u/%u\n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2 );

                topo.push_back(i0-1);
                topo.push_back(i1-1);
                topo.push_back(i2-1);
                normalIds.push_back(k0-1);
                normalIds.push_back(k1-1);
                normalIds.push_back(k2-1);
            }
        }
    }

    /*boundingSphere.C = 0.5*(minVal+maxVal);
    boundingSphere.R = sqrt((maxVal-minVal).sqrNorm())*0.5;*/

    fclose(f);
}

const glm::vec3 Mesh::getNormal(const glm::vec3& point) const{
    return glm::vec3(0,1,1);
    //return glm::normalize(glm::cross(v1-v0, v2-v0));
}

/*float Mesh::intersect(const Ray &ray)
{

    for(int i=0; i<topo.size(); i+=3){

        auto e1 = geom[i+1] - geom[i];
        auto e2 = geom[i+2] - geom[i];

        auto h = glm::cross(ray.direction, e2);
        auto a = glm::dot(e1, h);

        auto f = 1.f / a;
        auto s = ray.origin - triangle.v0;

        auto u = f * glm::dot(s, h);
        auto q = glm::cross(s, e1);
        auto v = f * glm::dot(ray.direction, q);
        auto t = f * glm::dot(e2, q);

        if(std::abs(a) < 0.00001)
            return noIntersect;
        if(u < 0 || u > 1)
            return noIntersect;
        if(v < 0 || (u+v) > 1)
            return noIntersect;
        if(t < 0)
            return noIntersect;
    }
    return t;

}*/



