// This code is highly based on smallpt
// http://www.kevinbeason.com/smallpt/
#include <cmath>
#include <algorithm>
#include <cassert>
#include <random>
#include <memory>
#include <fstream>
#include <iostream>
#include <float.h>

// GLM (vector / matrix)
#define GLM_FORCE_RADIANS

#include <glm/vec4.hpp>
#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>

//#include <mesh.h>

const float pi = 3.1415927f;
const float noIntersect = std::numeric_limits<float>::infinity();

const int N_RECURSION_RADIANCE_MAX = 4;
const float LUX = 3000;

//déclarations
struct Ray;
struct Triangle;
struct Sphere;
struct Box;
struct Mesh;

glm::vec3 radiance(const Ray &r, int nRecursion);
bool refract(glm::vec3 i, glm::vec3 n, float ior, glm::vec3 &wo);
float fresnelR(const glm::vec3 i, const glm::vec3 n, const float ior);
glm::vec3 sample_cos(const float u, const float v, const glm::vec3 n);
float random_u();

float intersect(const Ray & ray, const Triangle &triangle);
float intersect(const Ray & ray, const Sphere &sphere);
float intersect(const Ray &r, const Box& box);
float intersect(const Ray & ray, const Mesh &mesh);


bool isIntersect(float t)
{
    return t < noIntersect;
}

struct Ray
{
    const glm::vec3 origin, direction;
};

struct Sphere
{
    const float radius;
    const glm::vec3 center;

    const glm::vec3 getNormal(const glm::vec3& point) const{
        return glm::normalize(point-center);
    }
};

struct Triangle
{
    glm::vec3 v0, v1, v2;

    const glm::vec3 getNormal(const glm::vec3& point) const{
        return glm::normalize(glm::cross(v1-v0, v2-v0));
    }
};

struct Box
{
    glm::vec3 minBox, maxBox;

    bool inOut(const glm::vec3& p) const
    {
        if(p.x < minBox.x)
            return false;
        if(p.x > maxBox.x)
            return false;

        if(p.y < minBox.y)
            return false;
        if(p.y > maxBox.y)
            return false;

        if(p.z < minBox.z)
            return false;
        if(p.z > maxBox.z)
            return false;

        return true;
    }

    float intersectIn(const Ray& r) const
    {
        float tmax, tymax, tzmax;

        if(r.direction.x == 0)
            tmax = noIntersect;
        else if(r.direction.x > 0)
            tmax = (maxBox.x - r.origin.x) / r.direction.x;
        else
            tmax = (minBox.x - r.origin.x) / r.direction.x;

        if(r.direction.y == 0)
            tymax = noIntersect;
        else if(r.direction.y >= 0)
            tymax = (maxBox.y - r.origin.y) / r.direction.y;
        else
            tymax = (minBox.y - r.origin.y) / r.direction.y;

        if(tymax < tmax)
            tmax = tymax;


        if(r.direction.z == 0)
            return tmax;
        else if(r.direction.z > 0)
            tzmax = (maxBox.z - r.origin.z) / r.direction.z;
        else
            tzmax = (minBox.z - r.origin.z) / r.direction.z;

        if(tzmax < tmax)
            return tzmax;
        return tmax;
    }

    glm::vec3 getNormal(const glm::vec3& p) const{
        glm::vec3 cote = maxBox-minBox; //taille des coté
        glm::vec3 centre(minBox+cote*0.5f);
        glm::vec3 n(p-centre);
        n /= cote;  //la normal point dans la direction du cote le plus proche du point
        glm::vec3 na = glm::vec3(abs(n.x), abs(n.y), abs(n.z));

        int dir = 0;    //la normale est dans la direction de l'axe X
        if(na.x < na.y)
        {
            if(na.y < na.z)
                dir = 2;    //la normale est dans la direction de l'axe Z
            else
                dir = 1;    //la normale est dans la direction de l'axe Y
        }
        else if(na.x < na.z)
            dir = 2;        //la normale est dans la direction de l'axe Z

        switch(dir)
        {
        case 0: return (n.x < 0   ?   glm::vec3(-1,0,0)    :   glm::vec3(1,0,0));
        case 1: return (n.y < 0   ?   glm::vec3(0,-1,0)    :   glm::vec3(0,1,0));
        case 2: return (n.z < 0   ?   glm::vec3(0,0,-1)    :   glm::vec3(0,0,1));
        }
        return glm::vec3(0,0,0);
    }
};

struct BoiteEnglobante
{
    Box MaBoite;
    std::vector<BoiteEnglobante> boiteFilles;
    std::vector<Triangle*> triangleFilles;

    bool intersectBox(const Ray &r) const{
        return intersect(r, MaBoite) != noIntersect;
    }

    Triangle getTriangle(const Ray &r, float& distanceImpact) const{
        if(!triangleFilles.empty()){
            int indiceTriangle;
            for(int i=0; i < triangleFilles.size(); i++){
                if(intersect(r, *(triangleFilles[i]) ) < distanceImpact){
                    indiceTriangle = i;
                    distanceImpact = intersect(r, *(triangleFilles[i]) );
                }
            }
            return *triangleFilles[indiceTriangle];
        }

        if(!boiteFilles.empty()){
            float distanceImpact;
            int indiceBoite;
            for(int i=0; i<boiteFilles.size(); i++){
                if(distanceImpact > intersect(r, boiteFilles[i].MaBoite)){
                    boiteFilles[i].getTriangle(r, distanceImpact);
                    indiceBoite = i;
                }
            }
            return boiteFilles[indiceBoite].getTriangle(r, distanceImpact);
        }
    }

    void addTriangle(Triangle &t){
        triangleFilles.push_back(&t);
    }

    void subDivision(){

    }
};

struct Mesh
{
    std::vector<glm::vec3> geom;
    std::vector<int> topo;
    std::vector<glm::vec3> normalsPoints;
    std::vector<int> normalIds;
    std::vector<Triangle> listTriangle;

    BoiteEnglobante boite;
    glm::vec3 position;

    void translation(const float x, const float y, const float z){
        for(glm::vec3& p: geom){
            p.x+=x;
            p.y+=y;
            p.z+=z;
        }

        glm::vec3 transla(x, y, z);
        for(Triangle tri : listTriangle){
            tri.v0 += transla;
            tri.v1 += transla;
            tri.v2 += transla;
        }

        boite.MaBoite.minBox+= transla;
        boite.MaBoite.maxBox+= transla;
    }

    void rescale(float scale){

        translation(-position.x, -position.y, -position.z);

        for(glm::vec3& p: geom){
            p *= scale;
        }

        for(Triangle tri : listTriangle){
            tri.v0 *= scale;
            tri.v1 *= scale;
            tri.v2 *= scale;
        }

        boite.MaBoite.minBox *= scale;
        boite.MaBoite.maxBox *= scale;

        translation(position.x, position.y, position.z);
    }

    void rotate(const float& rX, const float& rY, const float& rZ){
        for (glm::vec3 vertex : geom) {
            if(rX != 0)
                rotateAboutAxis(vertex, rX, glm::vec3(1,0,0) );
            if(rY != 0)
                rotateAboutAxis(vertex, rY, glm::vec3(0,1,0) );
            if(rZ != 0)
                rotateAboutAxis(vertex, rZ, glm::vec3(0,0,1) );
        }

        for(Triangle tri : listTriangle){
            if(rX != 0){
                rotateAboutAxis(tri.v0, rX, glm::vec3(1,0,0) );
                rotateAboutAxis(tri.v1, rX, glm::vec3(1,0,0) );
                rotateAboutAxis(tri.v2, rX, glm::vec3(1,0,0) );
            }
            if(rY != 0){
                rotateAboutAxis(tri.v0, rY, glm::vec3(0,1,0) );
                rotateAboutAxis(tri.v1, rY, glm::vec3(0,1,0) );
                rotateAboutAxis(tri.v2, rY, glm::vec3(0,1,0) );
            }
            if(rZ != 0){
                rotateAboutAxis(tri.v0, rZ, glm::vec3(0,0,1) );
                rotateAboutAxis(tri.v1, rZ, glm::vec3(0,0,1) );
                rotateAboutAxis(tri.v2, rZ, glm::vec3(0,0,1) );
            }
        }
    }

    void rotateAboutAxis(glm::vec3& vertex, const float& angle, const glm::vec3& axis){
        float s = sinf(angle);
        float c = cosf(angle);
        float k = 1.0F - c;

        float nx = vertex.x * (c + k * axis.x * axis.x) + vertex.y * (k * axis.x * axis.y - s * axis.z)
                + vertex.z * (k * axis.x * axis.z + s * axis.y);
        float ny = vertex.x * (k * axis.x * axis.y + s * axis.z) + vertex.y * (c + k * axis.y * axis.y)
                + vertex.z * (k * axis.y * axis.z - s * axis.x);
        float nz = vertex.x * (k * axis.x * axis.z - s * axis.y) + vertex.y * (k * axis.y * axis.z + s * axis.x)
                + vertex.z * (c + k * axis.z * axis.z);

        vertex.x = nx;
        vertex.y = ny;
        vertex.z = nz;
    }

    void loadFromOBJ(const glm::vec3 &center, const char* obj){
        boite.MaBoite.minBox = glm::vec3(1E100, 1E100, 1E100);
        boite.MaBoite.maxBox = glm::vec3(-1E100, -1E100, -1E100);

        position = center;

        FILE* f = fopen(obj, "r");

        while (!feof(f)) {
            char line[255];
            fgets(line, 255, f);
            if (line[0]=='v' && line[1]==' ') {
                glm::vec3 vec;
                sscanf(line, "v %f %f %f\n", &vec[0], &vec[2], &vec[1]);
                vec[2] = -vec[2];
                glm::vec3 p = vec*50.f + center;
                geom.push_back(p);
                boite.MaBoite.maxBox[0] = std::max(boite.MaBoite.maxBox[0], p[0]);
                boite.MaBoite.maxBox[1] = std::max(boite.MaBoite.maxBox[1], p[1]);
                boite.MaBoite.maxBox[2] = std::max(boite.MaBoite.maxBox[2], p[2]);
                boite.MaBoite.minBox[0] = std::min(boite.MaBoite.minBox[0], p[0]);
                boite.MaBoite.minBox[1] = std::min(boite.MaBoite.minBox[1], p[1]);
                boite.MaBoite.minBox[2] = std::min(boite.MaBoite.minBox[2], p[2]);
            }
            if (line[0]=='v' && line[1]=='n') {
                glm::vec3 vec;
                sscanf(line, "vn %f %f %f\n", &vec[0], &vec[2], &vec[1]);
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
                for (unsigned int i=0; i < sizeof(line); i++){
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

        for(unsigned int i=0; i<topo.size(); i+=3){
            listTriangle.push_back( Triangle {geom[topo[i]], geom[topo[i+1]], geom[topo[i+2]]} );
        }

        for(Triangle trian : listTriangle){
            boite.addTriangle(trian);
        }
    }

    const glm::vec3 getNormal(const glm::vec3& point) const{
        return glm::vec3(0,1,1);
        //return glm::normalize(glm::cross(v1-v0, v2-v0));
    }
};


// WARRING: works only if r.d is normalized
float intersect (const Ray & ray, const Sphere &sphere)
{				// returns distance, 0 if nohit
    glm::vec3 op = sphere.center - ray.origin;		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    float t, b = glm:: dot(ray.direction, op), det =
            b * b - glm::dot(op, op) + sphere.radius * sphere.radius;
    if (det < 0)
        return noIntersect;
    else
        det = std::sqrt (det);
    return (t = b - det) >= 0 ? t : ((t = b + det) >= 0 ? t : noIntersect);
}

float intersect(const Ray & ray, const Triangle &triangle)
{
    auto e1 = triangle.v1 - triangle.v0;
    auto e2 = triangle.v2 - triangle.v0;

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

    return t;
}

float intersect(const Ray &r, const Box& box)
{
    if(box.inOut(r.origin))  {
        return 0;
    }

    float txmin, txmax, tymin, tymax, tzmin, tzmax;
    float div;

    if(r.direction.x == 0)    {
        txmin = FLT_MIN;
        txmax = FLT_MAX;
    }
    else if(r.direction.x > 0)    {
        div = 1 / r.direction.x;
        txmin = (box.minBox.x - r.origin.x) * div;
        txmax = box.maxBox.x * div;
    }
    else    {
        div = 1 / r.direction.x;
        txmin = (box.maxBox.x - r.origin.x) * div;
        txmax = (box.minBox.x - r.origin.x) * div;
    }

    if(r.direction.y == 0)    {
        tymin = FLT_MIN;
        tymax = FLT_MAX;
    }
    else if(r.direction.y >= 0)    {
        div = 1 / r.direction.y;
        tymin = box.minBox.y * div;
        tymax = (box.maxBox.y - r.origin.y) * div;
    }
    else    {
        div = 1 / r.direction.y;
        tymin = (box.maxBox.y - r.origin.y) * div;
        tymax = (box.minBox.y - r.origin.y) * div;
    }

    if( (txmin > tymax) || (tymin > txmax) )
        return noIntersect;

    if(r.direction.z == 0)    {
        tzmin = FLT_MIN;
        tzmax = noIntersect;
    }
    else if(r.direction.z > 0)    {
        div = 1 / r.direction.z;
        tzmin = (box.minBox.z - r.origin.z) * div;
        tzmax = (box.maxBox.z - r.origin.z) * div;
    }
    else    {
        div = 1 / r.direction.z;
        tzmin = (box.maxBox.z - r.origin.z) * div;
        tzmax = (box.minBox.z - r.origin.z) * div;
    }

    if( (txmin > tzmax) || (tymin > tzmax) || (tzmin > txmax) || (tzmin > tymax) )
        return noIntersect;

    return glm::length(glm::vec3(txmin, tymin, tzmin));
}

float intersect(const Ray & ray, const Mesh &mesh)
{
    if(intersect(ray, mesh.boite.MaBoite) == noIntersect) return noIntersect;

    float tFin = noIntersect, tTemp = noIntersect;

    for(Triangle tri : mesh.listTriangle){

        tTemp = intersect(ray, tri);

        if(tTemp != noIntersect){

            if(tFin == noIntersect){
                tFin = tTemp;
            }
            else if(tTemp < tFin){
                tFin = tTemp;
            }
        }
    }
    return tFin;
}

struct Diffuse
{
    const glm::vec3 color;
    const float coefficientSpeculaire;
    const glm::vec3 albedo() const { return color; }

    const glm::vec3 BSDF_Direct(const glm::vec3& c, const glm::vec3& n, const glm::vec3& l) const{

        return albedo();
    }

    const glm::vec3 BSDF_Indirect(const glm::vec3& c, const glm::vec3& n) const{

        return albedo();
    }

    glm::vec3 direct(const Ray& c, const glm::vec3& n, const Ray& l, const float& distanceLuxCarre) const {
        //direct = V(p, lampe) * BSDF_direct() * couleurLampe
        if( (glm::dot(n, l.direction))*glm::dot(n, c.direction) >= 0){
            float coeffLux = fabs(glm::dot(n, l.direction)/pi) * LUX / distanceLuxCarre;
            coeffLux += pow( fabs(glm::dot( reflect(l.direction, n), -c.direction )), coefficientSpeculaire );
            return coeffLux*BSDF_Direct(c.direction, n, l.direction);
        }
        return glm::vec3(0,0,0);
    }

    glm::vec3 indirect(const Ray& c, const glm::vec3& n, const glm::vec3& p, const int& nReccursion) const {

        if(nReccursion < N_RECURSION_RADIANCE_MAX){

            glm::vec3 colorSomme(0,0,0);
            int nMax = 1;

            for(int i=0; i < nMax; i++){

                glm::vec3 dSample;
                if(glm::dot(n, -c.direction) < 0)
                    dSample = sample_cos(random_u(), random_u(), -n);
                else
                    dSample = sample_cos(random_u(), random_u(), n);

                Ray sample{p+0.002f*dSample, dSample};
                colorSomme += BSDF_Indirect(c.direction, n) * radiance(sample, nReccursion+1);
            }

            return colorSomme/(float)nMax;
        }

        return glm::vec3(0,0,0);
    }


};

struct Glass
{
    const glm::vec3 color;
    const float coefficientSpeculaire;

    const glm::vec3 albedo() const { return color; }

    const glm::vec3 BSDF_Direct(const glm::vec3& c, const glm::vec3& n, const glm::vec3& l) const{

        return albedo();
    }

    const glm::vec3 BSDF_Indirect(const glm::vec3& c, const glm::vec3& n) const{

        return albedo();
    }

    glm::vec3 direct(const Ray& c, const glm::vec3& n, const Ray& l, const float& distanceLuxCarre) const {
        //direct = V(p, lampe) * BSDF_direct() * couleurLampe
        /*return BSDF_Direct(c.direction, n, l.direction) *
                (float)pow( fabs(glm::dot( reflect(l.direction, n), -c.direction )), coefficientSpeculaire ) *
                LUX / distanceLuxCarre;*/

        return glm::vec3(0,0,0);
    }

    glm::vec3 indirect(const Ray& c, const glm::vec3& n, const glm::vec3& p, const int& nReccursion) const {
        //indirect = alpha * BSDF_indirect * (radiance_reflet() + radiance_travers())

        if(nReccursion < N_RECURSION_RADIANCE_MAX){

            if(nReccursion <= 1){
                glm::vec3 dNewRayReflet = glm::normalize( reflect(c.direction, n) );
                Ray newRayReflet { p +(float)0.02*dNewRayReflet, dNewRayReflet};

                glm::vec3 reflectColor = radiance( newRayReflet, nReccursion+1);
                glm::vec3 refractColor(0,0,0);


                float fresnel = fresnelR(-c.direction, n, 1.33);
                glm::vec3 dNewRayTravers;

                // indice de refraction verre/air = 1.33
                if( refract(-c.direction, n, 1.33, dNewRayTravers) ){
                    Ray newRayTravers { p +(float)0.02*dNewRayTravers, dNewRayTravers};
                    refractColor = radiance( newRayTravers, nReccursion+1);
                }

                return BSDF_Indirect(c.direction, n) * ( fresnel*reflectColor + (1-fresnel)*refractColor) ;
            }
            else{

                float ran = random_u();
                float fresnel = fresnelR(-c.direction, n, 1.33);

                if(ran < fresnel){
                    glm::vec3 dNewRayReflet = glm::normalize( reflect(c.direction, n) );
                    Ray newRayReflet { p +(float)0.02*dNewRayReflet, dNewRayReflet};

                    return BSDF_Indirect(c.direction, n) * radiance( newRayReflet, nReccursion+1) * fresnel;

                }
                else{
                    glm::vec3 dNewRayTravers;

                    // indice de refraction verre/air = 1.33
                    refract(-c.direction, n, 1.33, dNewRayTravers);
                    Ray newRayTravers { p +(float)0.02*dNewRayTravers, dNewRayTravers};

                    return BSDF_Indirect(c.direction, n) * radiance( newRayTravers, nReccursion+1) * (1-fresnel);
                }
            }
        }

        return glm::vec3(0,0,0);
    }

};

struct Mirror
{
    const glm::vec3 color;
    const float coefficientSpeculaire;

    const glm::vec3 albedo() const { return color; }

    const glm::vec3 BSDF_Direct(const glm::vec3& c, const glm::vec3& n, const glm::vec3& l) const {

        /*
         * Direct :
         * BSDF_mirroir (C,Np,L) = 0
         */

        return albedo();
    }

    const glm::vec3 BSDF_Indirect(const glm::vec3& c, const glm::vec3& n) const {

        /*
        * Indirect :
        * BSDF_mirroir (C,Np,L) =
        */

        return albedo();
    }

    glm::vec3 direct(const Ray& c, const glm::vec3& n, const Ray& l, const float& distanceLuxCarre) const {
        //direct = V(p, lampe) * BSDF_direct() * couleurLampe
        /*return BSDF_Direct(c.direction, n, l.direction);
        (float)pow( fabs(glm::dot( reflect(l.direction, n), -c.direction )), coefficientSpeculaire ) *
                LUX / distanceLuxCarre;*/

        return glm::vec3(0,0,0);
    }

    glm::vec3 indirect(const Ray& c, const glm::vec3& n, const glm::vec3& p, const int& nReccursion) const {
        //indirect = alpha * BSDF_indirect * radiance()
        glm::vec3 retour(0,0,0);

        if(nReccursion < N_RECURSION_RADIANCE_MAX){
            glm::vec3 dNewRay = glm::normalize( reflect(c.direction, n) );
            Ray newRay { p +(float)0.02*dNewRay, dNewRay};

            retour = BSDF_Indirect(-c.direction, n) * radiance( newRay, nReccursion+1 );
        }

        return retour;
    }

};

struct Object
{
    virtual float intersect(const Ray &r) const = 0;
    virtual const glm::vec3 albedo() const = 0;
    virtual const glm::vec3 normal(const glm::vec3 &point) const = 0;
    virtual const glm::vec3 luxDirect(const Ray& c, const glm::vec3& n, const Ray& p, const float& distanceLuxCarre) const = 0;
    virtual const glm::vec3 luxIndirect(const Ray& c, const glm::vec3& n, const glm::vec3& p, const int& nReccursion) const = 0;

};

template<typename P, typename M>
struct ObjectTpl final : Object
{
    ObjectTpl(const P &_p, const M &_m)
        :primitive(_p), material(_m)
    {}

    float intersect(const Ray &ray) const
    {
        return ::intersect(ray, primitive);
    }

    const glm::vec3 albedo() const{
        return material.albedo();
    }

    const glm::vec3 normal(const glm::vec3 &point) const{
        return primitive.getNormal(point);
    }

    const glm::vec3 luxDirect(const Ray& c, const glm::vec3& n, const Ray& p, const float& distanceLuxCarre) const{
        return material.direct(c,n,p, distanceLuxCarre);
    }

    const glm::vec3 luxIndirect(const Ray& c, const glm::vec3& n, const glm::vec3& p, const int& nReccursion) const{
        return material.indirect(c,n,p, nReccursion);
    }

    const P &primitive;
    const M &material;
};


template<typename P, typename M>
std::unique_ptr<Object> makeObject(const P&p, const M&m)
{
    return std::unique_ptr<Object>(new ObjectTpl<P, M>{p, m});
}

// Scene
namespace scene
{
// Primitives

// Left Wall
const Triangle leftWallA{{0, 0, 0}, {0, 100, 0}, {0, 0, 150}};
const Triangle leftWallB{{0, 100, 150}, {0, 100, 0}, {0, 0, 150}};

// Right Wall
const Triangle rightWallA{{100, 0, 0}, {100, 100, 0}, {100, 0, 150}};
const Triangle rightWallB{{100, 100, 150}, {100, 100, 0}, {100, 0, 150}};

// Back wall
const Triangle backWallA{{0, 0, 0}, {100, 0, 0}, {100, 100, 0}};
const Triangle backWallB{{0, 0, 0}, {0, 100, 0}, {100, 100, 0}};

// Bottom Floor
const Triangle bottomWallA{{0, 0, 0}, {100, 0, 0}, {100, 0, 150}};
const Triangle bottomWallB{{0, 0, 0}, {0, 0, 150}, {100, 0, 150}};

// Top Ceiling
const Triangle topWallA{{0, 100, 0}, {100, 100, 0}, {0, 100, 150}};
const Triangle topWallB{{100, 100, 150}, {100, 100, 0}, {0, 100, 150}};

const Sphere leftSphere{16.5, glm::vec3 {27, 16.5, 47}};
const Sphere rightSphere{16.5, glm::vec3 {73, 16.5, 78}};

const glm::vec3 light{50, 70, 81.6};

// Materials
const Diffuse white{{.75, .75, .75}, 10};
const Diffuse red{{.75, .25, .25}, 10};
const Diffuse blue{{.25, .25, .75}, 10};

const Glass glass{{.9, .7, .9}, 20};
const Mirror mirror{{.9, .9, .5}, 0};

Mesh m, m2;//, m3;

// Objects
// Note: this is a rather convoluted way of initialising a vector of unique_ptr ;)
const std::vector<std::unique_ptr<Object>> objects = [] (){
    std::vector<std::unique_ptr<Object>> ret;
    //ret.push_back(makeObject(backWallA, white));
    //ret.push_back(makeObject(backWallB, white));
    //ret.push_back(makeObject(topWallA, white));
    //ret.push_back(makeObject(topWallB, white));
    //ret.push_back(makeObject(bottomWallA, white));
    //ret.push_back(makeObject(bottomWallB, white));
    //ret.push_back(makeObject(rightWallA, blue));
    //ret.push_back(makeObject(rightWallB, blue));
    //ret.push_back(makeObject(leftWallA, red));
    //ret.push_back(makeObject(leftWallB, red));

    //ret.push_back(makeObject(leftSphere, mirror));
    //ret.push_back(makeObject(rightSphere, glass));

    ret.push_back(makeObject(m, blue));
    ret.push_back(makeObject(m2, red));
    //ret.push_back(makeObject(m3, mirror));

    return ret;
}();
}

thread_local std::default_random_engine generator;
thread_local std::uniform_real_distribution<float> distribution(0.0,1.0);

float random_u()
{
    return distribution(generator);
}

glm::vec3 sample_cos(const float u, const float v, const glm::vec3 n)
{
    // Ugly: create an ornthogonal base
    glm::vec3 basex, basey, basez;

    basez = n;
    basey = glm::vec3(n.y, n.z, n.x);

    basex = glm::cross(basez, basey);
    basex = glm::normalize(basex);

    basey = glm::cross(basez, basex);

    // cosinus sampling. Pdf = cosinus
    return  basex * (std::cos(2.f * pi * u) * std::sqrt(1.f - v)) +
            basey * (std::sin(2.f * pi * u) * std::sqrt(1.f - v)) +
            basez * std::sqrt(v);
}

int toInt (const float x)
{
    return int (std::pow (glm::clamp (x, 0.f, 1.f), 1.f / 2.2f) * 255 + .5);
}

// WARNING: ASSUME NORMALIZED RAY
// Compute the intersection ray / scene.
// Returns true if intersection
// t is defined as the abscisce along the ray (i.e
//             p = r.o + t * r.d
// id is the id of the intersected object
Object* intersect (const Ray & r, float &t)
{
    t = noIntersect;
    Object *ret = nullptr;

    for(auto &object : scene::objects)
    {
        float d = object->intersect(r);
        if (isIntersect(d) && d < t)
        {
            t = d;
            ret = object.get();
        }
    }

    return ret;
}

// Reflect the ray i along the normal.
// i should be oriented as "leaving the surface"
glm::vec3 reflect(const glm::vec3 i, const glm::vec3 n)
{
    return n * (glm::dot(n, i)) * 2.f - i;
}

float sin2cos (const float x)
{
    return std::sqrt(std::max(0.0f, 1.0f-x*x));
}

// Fresnel coeficient of transmission.
// Normal point outside the surface
// ior is n0 / n1 where n0 is inside and n1 is outside
float fresnelR(const glm::vec3 i, const glm::vec3 n, const float ior)
{
    if(glm::dot(n, i) < 0)
        return fresnelR(i, n * -1.f, 1.f / ior);

    float R0 = (ior - 1.f) / (ior + 1.f);
    R0 *= R0;

    return R0 + (1.f - R0) * std::pow(1.f - glm::dot(i, n), 5.f);
}

// compute refraction vector.
// return true if refraction is possible.
// i and n are normalized
// output wo, the refracted vector (normalized)
// n point oitside the surface.
// ior is n00 / n1 where n0 is inside and n1 is outside
//
// i point outside of the surface
bool refract(glm::vec3 i, glm::vec3 n, float ior, glm::vec3 &wo)
{
    i = i * -1.f;

    if(glm::dot(n, i) > 0)
    {
        n = n * -1.f;
    }
    else
    {
        ior = 1.f / ior;
    }

    float k = 1.f - ior * ior * (1.f - glm::dot(n, i) * glm::dot(n, i));
    if (k < 0.)
        return false;

    wo = i * ior - n * (ior * glm::dot(n, i) + std::sqrt(k));

    return true;
}

glm::vec3 sample_sphere(const float r, const float u, const float v, float &pdf, const glm::vec3 normal)
{
    pdf = 1.f / (pi * r * r);
    glm::vec3 sample_p = sample_cos(u, v, normal);

    float cos = glm::dot(sample_p, normal);

    pdf *= cos;
    return sample_p * r;
}

glm::vec3 radiance (const Ray & r, int nRecursion = 0)
{
    float t = noIntersect;

    Object* obj = intersect(r, t);

    if(t != noIntersect){

        //*** Detection ombre ***//

        glm::vec3 pImpact = r.origin + t*r.direction;
        glm::vec3 n = glm::normalize(obj->normal(pImpact));


        glm::vec3 directColor(0,0,0);
        float pdf;
        int nbTestOmbre=1;

        for(int i=0; i < nbTestOmbre; i++){

            glm::vec3 pLight = scene::light + sample_sphere(9, random_u(), random_u(), pdf, glm::normalize(pImpact - scene::light));

            glm::vec3 dirOmbre = glm::normalize(pLight - pImpact);

            Ray rOmbre = {pImpact+ (float)0.018*dirOmbre, dirOmbre}; // "correction" de l'imprécision de position

            //Ray rOmbre = {pImpact, glm::normalize(scene::light - pImpact)};

            float tOmbre;
            float distanceLampeCarre =  glm::dot(pLight - pImpact, pLight - pImpact);

            intersect(rOmbre, tOmbre);
            if(tOmbre == noIntersect || tOmbre*tOmbre > distanceLampeCarre){

                directColor += obj->luxDirect(r, n, rOmbre, distanceLampeCarre);
            }
        }

        directColor /= (float)nbTestOmbre;

        /*
        float tOmbre;
        intersect(rOmbre, tOmbre);

        glm::vec3 dirOmbre = glm::normalize(scene::light - pImpact);
        Ray rOmbre = {pImpact+ (float)0.018*dirOmbre, dirOmbre}; // "correction" de l'imprécision de position

        float distanceLampeCarre =  glm::dot(scene::light - pImpact, scene::light - pImpact);

        if(tOmbre == noIntersect || tOmbre*tOmbre > distanceLampeCarre)
            directColor = obj->luxDirect(r, n, rOmbre, distanceLampeCarre);*/

        return directColor + obj->luxIndirect(r, n, pImpact, nRecursion);

    }else{
        return glm::vec3(0,0,0);
    }
}

int main (int, char **)
{
    scene::m.loadFromOBJ(glm::vec3(70,0,50), "Beautiful_Girl.obj");
    scene::m2 = scene::m;
    //scene::m3 = scene::m;

    scene::m.rescale(0.75);

    scene::m2.rescale(0.75);
    scene::m2.translation(-30,0,0);
    scene::m2.rotate(0,90,0);


    //scene::m3.translation(-60,0,0);


    int w = 768, h = 768;
    std::vector<glm::vec3> colors(w * h, glm::vec3{0.f, 0.f, 0.f});

    Ray cam {{50, 52, 295.6}, glm::normalize(glm::vec3{0, -0.042612, -1})};	// cam pos, dir
    float near = 1.f;
    float far = 10000.f;

    glm::mat4 camera =
            glm::scale(glm::mat4(1.f), glm::vec3(float(w), float(h), 1.f))
            * glm::translate(glm::mat4(1.f), glm::vec3(0.5, 0.5, 0.f))
            * glm::perspective(float(54.5f * pi / 180.f), float(w) / float(h), near, far)
            * glm::lookAt(cam.origin, cam.origin + cam.direction, glm::vec3(0, 1, 0))
            ;

    glm::mat4 screenToRay = glm::inverse(camera);

    for (int y = 0; y < h; y++)
    {
        std::cerr << "\rRendering: " << 100 * y / (h - 1) << "%";

        const int nAnti = 10;

#pragma omp parallel for schedule(dynamic, 1)
        for (unsigned short x = 0; x < w; x++)
        {
            glm::vec3 r(0,0,0);

            //antialiasing + //réduction bruit dù au indirect de diffus
            for(int i = 0; i < nAnti; i++){
                float u = random_u(), v = random_u();

                float rs = sqrt(-2*log(u));

                float delX = rs*cos(2*pi*v), delY = rs*sin(2*pi*v);

                glm::vec4 p0 = screenToRay * glm::vec4{float(x), float(h - y), 0.f, 1.f};
                glm::vec4 p1 = screenToRay * glm::vec4{float(x+delX*0.5), float(h - y+delY*0.5), 1.f, 1.f};


                glm::vec3 pp0 = glm::vec3(p0 / p0.w);
                glm::vec3 pp1 = glm::vec3(p1 / p1.w);


                glm::vec3 d = glm::normalize(pp1 - pp0);

                r += radiance (Ray{pp0, d});
            }
            r/=(float)nAnti;


            colors[y * w + x] += glm::clamp(r, glm::vec3(0.f, 0.f, 0.f), glm::vec3(1.f, 1.f, 1.f));
        }
    }

    {
        std::fstream f("image.ppm", std::fstream::out);
        f << "P3\n" << w << " " << h << std::endl << "255" << std::endl;

        for (auto c : colors)
            f << toInt(c.x) << " " << toInt(c.y) << " " << toInt(c.z) << " ";
        //f << toInt(pow(c.x, 1/2.2)) << " " << toInt(pow(c.y, 1/2.2)) << " " << toInt(pow(c.z, 1/2.2)) << " ";
    }
}
