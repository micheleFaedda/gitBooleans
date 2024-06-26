//
// Created by Michele on 20/06/24.
//

#ifndef GITBOOLEANS_INTERSECT_POINT_RATIONALS_H
#define GITBOOLEANS_INTERSECT_POINT_RATIONALS_H

#include <numerics.h>

class IntersectionPointRationals {

private:
    bigrational x,y,z;
    uint tri_id;

public:

    IntersectionPointRationals(bigrational x, bigrational y, bigrational z, uint tri_id){
        this->x = x;
        this->y = y;
        this->z = z;
        this->tri_id = tri_id;
    }

    bigrational getX() const {
        return x;
    }

    bigrational getY() const {
        return y;
    }

    bigrational getZ() const {
        return z;
    }

    uint getTriId() const {
        return tri_id;
    }

    void setX(bigrational x) {
        this->x = x;
    }

    void setY(bigrational y) {
        this->y = y;
    }

    void setZ(bigrational z) {
        this->z = z;
    }

    void setTriId(uint tri_id) {
        this->tri_id = tri_id;
    }

    bool operator==(const IntersectionPointRationals &p) const {
        return x == p.x && y == p.y && z == p.z;
    }

    bool operator!=(const IntersectionPointRationals &p) const {
        return x != p.x || y != p.y || z != p.z;
    }

    bool operator<(const IntersectionPointRationals &p) const {
        if(x < p.x) return true;
        if(x > p.x) return false;
        if(y < p.y) return true;
        if(y > p.y) return false;
        return z < p.z;
    }

    bool operator>(const IntersectionPointRationals &p) const {
        if(x > p.x) return true;
        if(x < p.x) return false;
        if(y > p.y) return true;
        if(y < p.y) return false;
        return z > p.z;
    }

    bool operator<=(const IntersectionPointRationals &p) const {
        return *this < p || *this == p;
    }

    bool operator>=(const IntersectionPointRationals &p) const {
        return *this > p || *this == p;
    }

    bool lessThanX(const IntersectionPointRationals &p) const {
        return x < p.x;
    }

    bool lessThanY(const IntersectionPointRationals &p) const {
        return y < p.y;
    }

    bool lessThanZ(const IntersectionPointRationals &p) const {
        return z < p.z;
    }

    bool greaterThanX(const IntersectionPointRationals &p) const {
        return x > p.x;
    }

    bool greaterThanY(const IntersectionPointRationals &p) const {
        return y > p.y;
    }

    bool greaterThanZ(const IntersectionPointRationals &p) const {
        return z > p.z;
    }


};

#endif //GITBOOLEANS_INTERSECT_POINT_RATIONALS_H
