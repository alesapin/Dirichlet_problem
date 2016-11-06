#ifndef _ITER_H_
#define _ITER_H_
#include "mesh.h"
#include <mpi.h>

class IIterator {
protected:
    const Function func;
    long itercount;
public:
    IIterator(const Function &func, const Mesh& m):func(func),itercount(0) {}
    virtual double iterate(Mesh &mesh) = 0;
};

class CGMPuassonIterator : public IIterator {
private:
    Mesh rMesh;
    Mesh gMesh;
    double zeroIteration(Mesh &pMesh);
    double tau;
public:
    CGMPuassonIterator(Function func, const Mesh& m):
        IIterator(func,m),
        gMesh(m.getLeftDownCorner(), m.getRightUpCorner(), m.getRows(),m.getColumns(),m.getParentRows(), m.getParentCols(),m.getRowsShift(), m.getColumnsShift()) ,
        rMesh(m.getLeftDownCorner(), m.getRightUpCorner(), m.getRows(),m.getColumns(),m.getParentRows(), m.getParentCols(),m.getRowsShift(), m.getColumnsShift())  {}
    virtual double iterate(Mesh &pMesh);
};

class CGMPuassonIteratorMPI : public IIterator {
private:
    int left, right, up, down;
    Mesh rMesh;
    Mesh gMesh;
    double zeroIteration(Mesh &pMesh);
    double tau;
    void getMeshBorders(Mesh &mesh);
    int rank, size;
    inline bool checkBorder(const Mesh &mesh, long i, long j) {
        return (
                mesh.getRowsShift() + i - 1 >= 0
                && mesh.getRowsShift() + i + 1 < mesh.getParentRows()
                && mesh.getColumnsShift() + j - 1 >= 0
                && mesh.getColumnsShift() + j + 1 < mesh.getParentCols()
               );
    }

public:
    CGMPuassonIteratorMPI(Function func, const Mesh& m, int rank, int left, int right, int up, int down, int size):
        IIterator(func,m),
        gMesh(m.getLeftDownCorner(), m.getRightUpCorner(), m.getRows(),m.getColumns(), m.getParentRows(), m.getParentCols(), m.getRowsShift(), m.getColumnsShift()),
        rMesh(m.getLeftDownCorner(), m.getRightUpCorner(), m.getRows(),m.getColumns(), m.getParentRows(), m.getParentCols(), m.getRowsShift(), m.getColumnsShift()),
        rank(rank), left(left), right(right), up(up), down(down), size(size)
    {}
    virtual double iterate(Mesh &pMesh);
};
#endif
