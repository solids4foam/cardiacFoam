#ifndef UINIT_H
#define UINIT_H

#include <cmath>

// Define F(x, y, z) depending on dimension
inline double F(double x, double y, double z, int dimension)
{
    if (dimension == 1) {
        return std::cos(M_PI * x);
        
    } else if (dimension == 2) {
        return std::cos(M_PI * x) * std::cos(2.0 * M_PI * y);
    } else { // 3D
        return std::cos(M_PI * x) * std::cos(2.0 * M_PI * y) * std::cos(3.0 * M_PI * z);
    }
}

// Define G(x, y, z) depending on dimension
inline double G(double x, double y, double z, int dimension)
{
    if (dimension == 1) {
        return 1.0 + x;
    } else if (dimension == 2) {
        return 1.0 + x * y * y;
    } else { // 3D
        return 1.0 + x * y * y * z * z * z;
    }
}

// Initialization function for any dimension
inline void u_init(double x, double y, double z, double &V, double &u1, double &u2, double &u3, int dimension)
{
    double f = F(x, y, z, dimension);
    double g = G(x, y, z, dimension);

    V  = f;
    u1 = g + f;
    u2 = std::pow(g, -0.5);
    u3 = 0.0;
}

#endif
