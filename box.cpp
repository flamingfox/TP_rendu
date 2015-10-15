#include "box.h"


Box::Box()  : min(vec3(FLT_MAX,FLT_MAX,FLT_MAX)), max(vec3(FLT_MIN,FLT_MIN,FLT_MIN))
{
}

Box::Box(const vec3& _min, const vec3& _max)    :   min(_min),  max(_max)
{
}

Box::Box(const std::vector<vec3>& points)
{
    parcourtPoints(points);
}


void Box::updatePoint(const vec3& p)
{
    update(p);
}

bool Box::inOut(const vec3& p) const
{
    for(int i = 0;  i < 3;  i++)
    {
        if(p[i] < min[i])
            return false;
        if(p[i] > max[i])
            return false;
    }
    return true;
}


bool Box::intersect(const Rayon &r, float &distanceMin, float &distanceMax) const
{
    if(this->inOut(r.getOrigine()))  {
        distanceMin = 0;
        distanceMax = intersectIn(r);
        return true;
    }

    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    float div;

    if(r.getDirection().x == 0)    {
        tmin = FLT_MIN;
        tmax = FLT_MAX;
    }
    else if(r.getDirection().x > 0)    {
        div = 1 / r.getDirection().x;
        tmin = (min.x - r.getOrigine().x) * div;
        tmax = (max.x - r.getOrigine().x) * div;
    }
    else    {
        div = 1 / r.getDirection().x;
        tmin = (max.x - r.getOrigine().x) * div;
        tmax = (min.x - r.getOrigine().x) * div;
    }

    if(r.getDirection().y == 0)    {
        tymin = FLT_MIN;
        tymax = FLT_MAX;
    }
    else if(r.getDirection().y >= 0)    {
        div = 1 / r.getDirection().y;
        tymin = (min.y - r.getOrigine().y) * div;
        tymax = (max.y - r.getOrigine().y) * div;
    }
    else    {
        div = 1 / r.getDirection().y;
        tymin = (max.y - r.getOrigine().y) * div;
        tymax = (min.y - r.getOrigine().y) * div;
    }

    if( (tmin > tymax) || (tymin > tmax) )
        return false;

    if(tymin > tmin)
        tmin = tymin;

    if(tymax < tmax)
        tmax = tymax;


    if(r.getDirection().z == 0)    {
        tzmin = FLT_MIN;
        tzmax = FLT_MAX;
    }
    else if(r.getDirection().z > 0)    {
        div = 1 / r.getDirection().z;
        tzmin = (min.z - r.getOrigine().z) * div;
        tzmax = (max.z - r.getOrigine().z) * div;
    }
    else    {
        div = 1 / r.getDirection().z;
        tzmin = (max.z - r.getOrigine().z) * div;
        tzmax = (min.z - r.getOrigine().z) * div;
    }

    if( (tmin > tzmax) || (tzmin > tmax) )
        return false;

    if(tzmin > tmin)
        tmin = tzmin;

    if(tzmax < tmax)
        tmax = tzmax;

    if(tmin>=0)
        distanceMin = tmin;
    //else
    //    return false; //inutile apparament
    //distanceMin += 0.002;

    if(tmax>0)
        distanceMax = tmax;

    return true;
}



inline float Box::intersectIn(const Rayon& r) const
{
    float tmax, tymax, tzmax;

    if(r.getDirection().x == 0)
        tmax = FLT_MAX;
    else if(r.getDirection().x > 0)
        tmax = (max.x - r.getOrigine().x) / r.getDirection().x;
    else
        tmax = (min.x - r.getOrigine().x) / r.getDirection().x;

    if(r.getDirection().y == 0)
        tymax = FLT_MAX;
    else if(r.getDirection().y >= 0)
        tymax = (max.y - r.getOrigine().y) / r.getDirection().y;
    else
        tymax = (min.y - r.getOrigine().y) / r.getDirection().y;

    if(tymax < tmax)
        tmax = tymax;


    if(r.getDirection().z == 0)
        return tmax;
    else if(r.getDirection().z > 0)
        tzmax = (max.z - r.getOrigine().z) / r.getDirection().z;
    else
        tzmax = (min.z - r.getOrigine().z) / r.getDirection().z;

    if(tzmax < tmax)
        return tzmax;
    return tmax;
}


/**********************inline********************/

inline void Box::setDefaultBox()
{
    min = vec3(FLT_MAX,FLT_MAX,FLT_MAX);
    max = vec3(FLT_MIN,FLT_MIN,FLT_MIN);
}

inline void Box::update(const vec3& p)
{
    updateMin(p);
    updateMax(p);
}

inline void Box::updateMin(const vec3& p)
{
    for(int i = 0;  i < 3;  i++)
        if(p[i] < min[i])
            min[i] = p[i];

}

inline void Box::updateMax(const vec3& p)
{
    for(int i = 0;  i < 3;  i++)
        if(p[i] > max[i])
            max[i] = p[i];
}

inline void Box::parcourtPoints(const std::vector<vec3>& points)
{
    if(points.empty())
        setDefaultBox();
    else if(points.size() == 1)
    {
        min = points[0];
        max = points[0];
    }
    else
    {
        std::vector<vec3>::const_iterator it = points.begin();
        min = *it;
        max = points[points.size()-1];

        ++it;
        for(;  it != points.end()-1; ++it)
            update(*it);
    }
}

float Box::diffZ() const
{
    return max.z-min.z;
}

void Box::merge(const Box& box2)
{
    updateMin(box2.min);
    updateMax(box2.max);
}


void Box::operator+=(const vec3& t)
{
    min += t;
    max += t;
}
