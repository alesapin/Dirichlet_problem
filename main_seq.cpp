#include <iostream>
#include <cmath>
#include "mesh.h"
#include "iter.h"


double F(PointD p) {
    double x = p.first;
    double y = p.second;
    return (x*x + y*y) / (4*std::pow(4+x*y,1.5));
}
double phi(PointD p){
    double x = p.first;
    double y = p.second;
    return sqrt(4 + x*y);
}

double evaluation(Function p, const Mesh &m) {
    double diff = 0;
    for (long i = 0; i < m.getRows(); ++i){
        for(long j = 0; j < m.getColumns(); ++j) {
            diff += std::fabs(p(PointD(m.getPoint(i,j))) - m(i,j));
        }
    }
    return diff;
}

int main() {
    double A = 4, B = 4;
    int rows = 20, cols = 20;
    Mesh result(PointUi(0,0), PointUi(A,B), rows, cols, rows, cols);

    initMeshBoundaries(result, phi);
//    std::cerr << "(" <<result.getPoint(1,1).first << ","<< result.getPoint(1,1).second  << ")"<< "\n";
//    std::cerr << "(" <<result.getPoint(2,1).first << ","<< result.getPoint(2,1).second  << ")"<< "\n";
//    std::cerr << "(" <<result.getPoint(1,2).first << ","<< result.getPoint(1,2).second  << ")"<< "\n";
    //std::cerr << fiveDotScheme(result,1,1) << "\n";
    int iters = 1000;
    Mesh real(PointUi(0,0), PointUi(A,B), rows, cols, rows, cols);
    for (long i = 0; i < real.getRows(); ++i) {
        for (long j = 0; j < real.getColumns(); ++j) {
            real(i,j) = phi(real.getPoint(i,j));
        }
    }
    std::cerr << "Real:\n";
    std::cerr << real << "\n";
//    std::cerr << fiveDotScheme(real, 1, 1) << "\n";
    initMeshBoundaries(result, phi);
    CGMPuassonIterator cgm(F, result);
    for (int i = 0; i < iters; ++i) {
        cgm.iterate(result);
        //std::cerr << "RealError: " << i << " Error: " << evaluation(phi, result) << "\n";
    }

    std::cerr << "result:\n";
    std::cerr << result << "\n";
    return 0;
}
