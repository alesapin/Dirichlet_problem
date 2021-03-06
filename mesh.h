#ifndef _MESH_H_
#define _MESH_H_
#include <cmath>
#include <map>
#include <cfloat>
#include <vector>
#include "mat.h"

typedef std::pair<std::size_t, std::size_t> PointUi;
typedef std::pair<double, double> PointD;
static const PointD NAN_POINT(DBL_MAX, DBL_MAX);
typedef double (*Function)(PointD);

class Mesh {
    private:
        Mat data;
        Vec up;
        Vec down;
        Vec left;
        Vec right;
        PointUi leftDownCorner;
        PointUi rightUpCorner;
        std::vector<std::vector<PointD> >  pointCache;
        long rowsShift, colsShitf, parentRows, parentCols;
        static const double COEFF;
        static const double COEFF_DENOM;
        inline static double f(double t) {
            return (std::pow(1+t, COEFF) - 1) / COEFF_DENOM;
        }
        PointD countPoint(long i, long j) const;
        void initPointCache();
    public:
        Mesh(PointUi leftDownCorner, PointUi rightUpCorner, long rows, long cols, long parentRows, long parentCols, long rowsShift=0, long colsShitf=0);
        Mesh(PointUi leftDownCorner, PointUi rightUpCorner, long rows, long cols, long parentRows, long parentCols, double *data, long rowsShift=0, long colsShitf=0);
        Mesh(PointUi leftDownCorner, PointUi rightUpCorner, long parentRows, long parentCols, const Mat &data, long rowsShift =0, long colsShitf = 0);
        Mesh():leftDownCorner(0,0), rightUpCorner(0,0), rowsShift(0), colsShitf(0), parentRows(0),parentCols(0){}
        Mesh(const Mesh &other):
            data(other.data),
            leftDownCorner(other.leftDownCorner), 
            rightUpCorner(other.rightUpCorner),
            rowsShift(other.rowsShift),
            colsShitf(other.colsShitf),
            parentRows(other.parentRows),
            parentCols(other.parentCols),
            pointCache(other.pointCache)
        {
        }
        Mesh &operator=(const Mesh &other) {
            data = other.data;
            leftDownCorner = other.leftDownCorner; 
            rightUpCorner = other.rightUpCorner;
            rowsShift = other.rowsShift;
            colsShitf = other.colsShitf;
            parentRows = other.parentRows;
            parentCols = other.parentCols;
            pointCache = other.pointCache;
            return *this;
        }
        PointD getPoint(long i, long j) const;
        PointD getHShtr(long i, long j) const;
        double operator()(long i, long j) const {
            return data(i,j);
        }
        double &operator()(long i, long j) { return data(i,j); }
        long getRows() const { return data.rowsCount(); }
        long getColumns() const { return data.colsCount(); }
        void addUp(const Vec &up) {
            this->up = up;
        };
        void addDown(const Vec &down) {
            this->down = down;
        };
        void addLeft(const Vec &left) {
            this->left = left;
        };
        void addRight(const Vec &right) {
            this->right = right;
        }
        Vec getUpRow() {return data.getRow(0);}
        Vec getDownRow() {return data.getRow(data.rowsCount() - 1);}
        Vec getLeftCol() {return data.getCol(0);}
        Vec getRightCol() {return data.getCol(data.colsCount() - 1);}
        double *getData() {return data.barememptr();}
        friend double fiveDotScheme(const Mesh &mesh, long i, long j);
        friend bool checkInGlobalBorder(const Mesh &mesh, long i, long j);
        friend void initMeshBoundaries(Mesh &mesh, Function phi);
        friend std::map<int, Mesh> splitMesh(const Mesh &mesh, long procnum);
        friend Mesh collectMesh(const std::map<int, Mesh> &submeshs);
        friend double getEuError(const Mesh &mesh);
        PointUi getLeftDownCorner() const {return leftDownCorner; }
        PointUi getRightUpCorner() const {return rightUpCorner; }
        long getRowsShift() const { return rowsShift; }
        long getColumnsShift() const {return colsShitf; }
        long getParentRows() const {return parentRows; }
        long getParentCols() const {return parentCols; }
        friend std::ostream &operator<< (std::ostream &os, const Mesh& m);
        friend std::ostream &dropToStream(std::ostream& os, const Mesh &mesh);
};
#endif
