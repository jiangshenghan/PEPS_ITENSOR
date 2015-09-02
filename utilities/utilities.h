
#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "predef.h"

inline std::ostream &operator<<(std::ostream &s, const Coordinate &coord)
{
    return s << "(" << coord[0] << "," << coord[1] << "," << coord[2] << ")";
}

#endif
