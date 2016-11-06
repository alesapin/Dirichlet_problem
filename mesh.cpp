#include "mesh.h"
#include <iostream>
#include <iomanip>

namespace {
    PointUi splitFunction(int N0, int N1, int p) {
        double n0, n1;
        int p0, i;

        n0 = (double) N0; n1 = (double) N1;
        p0 = 0;

        for(i = 0; i < p; i++) {
            if(n0 > n1) {
                n0 = n0 / 2.0;
                ++p0;
            } else {
                n1 = n1 / 2.0;
            }
        }
        return PointUi(p0, p-p0);
    }

}

const double Mesh::COEFF = 1.5;
const double Mesh::COEFF_DENOM = 1.82842712474;
Mesh::Mesh(PointUi leftDownCorner, PointUi rightUpCorner, long rows, long cols, long parentRows, long parentCols, long rowsShift, long colsShitf):
    leftDownCorner(leftDownCorner),
    rightUpCorner(rightUpCorner),
    rowsShift(rowsShift),
    colsShitf(colsShitf),
    parentRows(parentRows),
    parentCols(parentCols),
    data(rows, cols, 0){
}

Mesh::Mesh(PointUi leftDownCorner, PointUi rightUpCorner, long rows, long cols,  long parentRows, long parentCols, double *data, long rowsShift, long colsShitf):
    leftDownCorner(leftDownCorner),
    rightUpCorner(rightUpCorner),
    rowsShift(rowsShift),
    colsShitf(colsShitf),
    parentRows(parentRows),
    parentCols(parentCols),
    data(data, rows, cols) {
}
Mesh::Mesh(PointUi leftDownCorner, PointUi rightUpCorner, long parentRows, long parentCols, const Mat &data, long rowsShift, long colsShitf):
    leftDownCorner(leftDownCorner),
    rightUpCorner(rightUpCorner),
    rowsShift(rowsShift),
    colsShitf(colsShitf),
    parentRows(parentRows),
    parentCols(parentCols),
    data(data) {
}
PointD Mesh::getPoint(long i, long j) const {
    double dblI = static_cast<double>(i) + rowsShift;
    double dblJ = static_cast<double>(j) + colsShitf;
    double xcoeff = dblI / getParentRows();
    double ycoeff = dblJ / getParentCols();
    double x = rightUpCorner.first * f(xcoeff) + leftDownCorner.first * (1 - f(xcoeff));
    double y = rightUpCorner.second * f(ycoeff) + leftDownCorner.second * (1 - f(ycoeff));
    return PointD(x,y);
}

double fiveDotScheme(const Mesh &mesh, long i, long j) {
    PointD prevPoint = mesh.getPoint(i-1, j-1);
    PointD curPoint = mesh.getPoint(i,j);
    PointD nextPoint = mesh.getPoint(i+1, j+1);
    double h1 = nextPoint.first - curPoint.first;
    double h2 = nextPoint.second - curPoint.second;
    double hprev1 = curPoint.first - prevPoint.first;
    double hprev2 = curPoint.second - prevPoint.second;
    double hs1 = (h1 + hprev1) / 2;
    double hs2 = (h2 + hprev2) / 2;

    double leftPoint = j - 1 < 0 ? mesh.left(i):mesh(i, j-1);
    double rightPoint = j + 1 >= mesh.getColumns() ? mesh.right(i) : mesh(i, j+1);
    double downPoint = i + 1 >= mesh.getRows() ? mesh.down(j):mesh(i+1, j);
    double upPoint = i - 1 < 0 ? mesh.up(j):mesh(i-1, j);
    double ypart = 1 / hs1 * ((mesh(i,j) - downPoint)/hprev1 - (upPoint - mesh(i,j)/h1));
    double xpart = 1 / hs2 * ((mesh(i,j) - leftPoint)/hprev2 - (rightPoint - mesh(i,j)/h2));
    return xpart + ypart;
}


std::ostream &operator<< (std::ostream &os, const Mesh& m) {
    for(int i = 0; i<m.getRows(); ++i){
        for(int j = 0; j < m.getColumns(); ++j){
            os <<std::setw(10) << m(i,j) << "\t";
        }
        os << "\n";
    }
    return os;
}
void initMeshBoundaries(Mesh &mesh,Function f) {
    for (std::size_t i = 0; i < mesh.getRows() ; ++i ){
        mesh(i,0) = f(mesh.getPoint(i,0));
        mesh(i, mesh.getColumns() - 1) = f(mesh.getPoint(i, mesh.getColumns() - 1));
    } for (std::size_t j = 0; j < mesh.getColumns(); ++j) { mesh(0,j) = f(mesh.getPoint(0,j));
        mesh(mesh.getRows() - 1,j) = f(mesh.getPoint(mesh.getRows() - 1, j));
    }
}

bool checkInGlobalBorder(const Mesh &mesh, long i, long j) {
    return (
        i + mesh.rowsShift >= 0
        && i + mesh.rowsShift < mesh.parentRows
        && j + mesh.colsShitf >= 0
        && j + mesh.colsShitf < mesh.parentCols
    );
}

std::map<int, Mesh> splitMesh(const Mesh &mesh, long procPower) {
    PointUi size = splitFunction(mesh.getRows(), mesh.getColumns(), procPower);
    long rows = 1 << size.first;
    long cols = 1 << size.second;
    std::map<int, Mesh> result;
    std::map<std::pair<int,int>, Mat> submats = split(mesh.data, rows, cols);
    int procnum = 0;
    long shiftRows = -1, shiftCols = -1;
    for(std::map<std::pair<int,int>,Mat>::iterator it = submats.begin(); it != submats.end(); ++it) {
       int r = it->first.first;
       int c = it->first.second;
       if(shiftRows == -1){
           shiftRows = it->second.rowsCount();
           shiftCols = it->second.colsCount();
       }
       result[procnum++] = Mesh(mesh.leftDownCorner,
               mesh.rightUpCorner, mesh.getRows(), mesh.getColumns(), it->second, r*shiftRows, c*shiftCols);
    }
    return result;
}
