#include <iostream>
#include <cmath>
#include "mesh.h"
#include "iter.h"
#include <getopt.h>
#include <cstdlib>
#include <fstream>
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
struct Params {
    long a;
    long b;
    long rows;
    long cols;
    std::string fname;
};

Params getParams(int argc, char** argv) {
      int c;
      char *opt_val = NULL;
      Params result;
      while ((c = getopt (argc, argv, "a:b:r:c:f:")) != -1) {
          switch(c) {
            case 'a':
                opt_val = optarg;
                result.a = strtol(opt_val, NULL, 10);
                break;
            case 'b':
                opt_val = optarg;
                result.b = strtol(opt_val, NULL, 10);
                break;
            case 'r':
                opt_val = optarg;
                result.rows = strtol(opt_val, NULL, 10);
                break;
            case 'c':
                opt_val = optarg;
                result.cols = strtol(opt_val, NULL, 10);
                break;
            case 'f':
                opt_val = optarg;
                result.fname = std::string(opt_val);
                break;
          }
      }
      return result;
}

const double EPSILON = 0.0001;

int main(int argc, char **argv) {
    Params pars = getParams(argc, argv);
    long A = pars.a, B = pars.b;
    int totalRows = pars.cols, totalCols = pars.rows;
    Mesh result(PointUi(0,0), PointUi(A,B), totalRows, totalCols, totalRows, totalCols);
    initMeshBoundaries(result, phi);
    CGMPuassonIterator cgm(F, result);
    time_t start,end;
    time(&start);
    double error = cgm.iterate(result);
    int iter = 1;
    while(error > EPSILON){
        std::cerr << "RealError: " << (error=cgm.iterate(result)) << " iter: " << iter++ <<"\n";
    }
    time (&end);
    double diff = difftime (end,start);
    std::ofstream ofs(pars.fname.c_str());
    ofs << diff << "\n";
    return 0;
}
