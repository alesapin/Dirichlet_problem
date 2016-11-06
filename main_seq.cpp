#include <iostream>
#include <cmath>
#include "mesh.h"
#include "iter.h"

double F(PointD p) {
    double x = p.first;
    double y = p.second;
    double mult = sqrt(4 + x*y);
    return (x*x + y*y) / (4*mult*mult*mult);
}
double phi(PointD p){
    double x = p.first;
    double y = p.second;
    return sqrt(4 + x*y);
}

int main() {
    double A = 4, B = 4;
    int rows = 2000, cols = 2000;
    Mesh result(PointUi(0,0), PointUi(A,B), rows, cols, rows, cols);
    int iters = 100;
    initMeshBoundaries(result, phi);
    CGMPuassonIterator cgm(F, result);
    for (int i = 0; i < iters; ++i) {
        std::cerr << "CGMIteration: " << i << " Error: " << cgm.iterate(result) << "\n";
    }
    return 0;
}
