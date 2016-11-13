#include "iter.h"
double CGMPuassonIterator::zeroIteration(Mesh& pMesh) {
    for (long i = 1; i < pMesh.getRows() - 1; ++i) {
        for (long j = 1; j < pMesh.getColumns() - 1; ++j) {
            rMesh(i,j) = fiveDotScheme(pMesh, i, j) - func(pMesh.getPoint(i,j));
            gMesh(i,j) = rMesh(i,j);
        }
    }

    double numerator = 0, denumerator = 0;
    int counter = 0;
    for (std::size_t i = 1; i < pMesh.getRows() - 1; ++i) {
        for (std::size_t j = 1; j < pMesh.getColumns() - 1; ++j) {
           numerator += rMesh(i,j) * rMesh(i,j);
           denumerator += fiveDotScheme(rMesh,i,j) * rMesh(i,j);
           counter++;
        }
    }
    tau = numerator / denumerator;
    return 100000;
}

double CGMPuassonIteratorMPI::zeroIteration(Mesh &pMesh) {

    getMeshBorders(pMesh);

    for (long i = 0; i < pMesh.getRows(); ++i) {
        for (long j = 0; j < pMesh.getColumns(); ++j) {
            if(checkBorder(rMesh, i, j)) {
                rMesh(i,j) = fiveDotScheme(pMesh, i, j) - func(pMesh.getPoint(i,j));
                gMesh(i,j) = rMesh(i,j);
            }
        }
    }

    getMeshBorders(rMesh);
    getMeshBorders(gMesh);
    double numerator = 0, denumerator = 0;
    int counter = 0;
    for (std::size_t i = 0; i < pMesh.getRows(); ++i) {
        for (std::size_t j = 0; j < pMesh.getColumns(); ++j) {
            if (checkBorder(pMesh, i, j)) {
                numerator += rMesh(i,j) * rMesh(i,j);
                denumerator += fiveDotScheme(rMesh,i,j) * rMesh(i,j);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double allNumerator, allDenumerator;
    MPI_Allreduce(&numerator, &allNumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&denumerator, &allDenumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    tau = allNumerator/allDenumerator;
    return 100000;
}


void CGMPuassonIteratorMPI::getMeshBorders(Mesh &mesh) {

    MPI_Status status[4];
    MPI_Request send[4];
    MPI_Request recv[4];
    std::tr1::shared_ptr<double> upBuf, downBuf, rightBuf, leftBuf;
    if (up >= 0 && up < size) {
        upBuf = std::tr1::shared_ptr<double>(new double[mesh.getColumns()]);
        Vec upCur = mesh.getUpRow();
        MPI_Isend(upCur.memptr(), upCur.size(), MPI_DOUBLE, up, 0, MPI_COMM_WORLD, &send[0]);
        MPI_Irecv(upBuf.get(), mesh.getColumns(), MPI_DOUBLE, up, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[0]);

        MPI_Wait(&send[0],&status[0]);
        MPI_Wait(&recv[0],&status[0]);
    }
    if(down >= 0 && down < size) {
        downBuf = std::tr1::shared_ptr<double> (new double[mesh.getColumns()]);
        Vec downCur = mesh.getDownRow();
        MPI_Isend(downCur.memptr(), downCur.size(), MPI_DOUBLE, down, 0, MPI_COMM_WORLD, &send[1]);
        MPI_Irecv(downBuf.get(), mesh.getColumns(), MPI_DOUBLE, down, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[1]);

        MPI_Wait(&send[1],&status[1]);
        MPI_Wait(&recv[1],&status[1]);
    }
    if (right >= 0 && right < size) {
        rightBuf = std::tr1::shared_ptr<double> (new double[mesh.getRows()]);
        Vec rightCur = mesh.getRightCol();
        MPI_Isend(rightCur.memptr(), rightCur.size(), MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &send[2]);
        MPI_Irecv(rightBuf.get(), mesh.getRows(), MPI_DOUBLE, right, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[2]);

        MPI_Wait(&send[2],&status[2]);
        MPI_Wait(&recv[2],&status[2]);
    }
    if(left >= 0 && left < size) {
        leftBuf = std::tr1::shared_ptr<double> (new double[mesh.getRows()]);
        Vec leftCur = mesh.getLeftCol();
        MPI_Isend(leftCur.memptr(), leftCur.size(), MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &send[3]);
        MPI_Irecv(leftBuf.get(), mesh.getRows(), MPI_DOUBLE, left, MPI_ANY_TAG, MPI_COMM_WORLD, &recv[3]);

        MPI_Wait(&send[3],&status[3]);
        MPI_Wait(&recv[3],&status[3]);
    }
    if (upBuf) {
        Vec upR(upBuf, mesh.getColumns());
        mesh.addUp(upR);
    }
    if (downBuf) {
        Vec downR(downBuf, mesh.getColumns());
        mesh.addDown(downR);
    }
    if (rightBuf) {
        Vec rightR(rightBuf, mesh.getRows());
        mesh.addRight(rightR);
    }
    if (leftBuf) {
        Vec leftR(leftBuf, mesh.getRows());
        mesh.addLeft(leftR);
    }
}

double CGMPuassonIteratorMPI::iterate(Mesh &pMesh) {
    if(itercount == 0){
        itercount++;
        return zeroIteration(pMesh);
    }

    long rows = pMesh.getRows();
    long cols = pMesh.getColumns();
    double err = 0;
    for (std::size_t i = 0; i < rows ; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            if (checkBorder(rMesh, i, j)) {
                double val = pMesh(i,j) - tau*gMesh(i,j);
                double errV = pMesh(i,j) - val;
                PointD hs = pMesh.getHShtr(i,j);
                err += errV*errV*hs.first*hs.second;
                pMesh(i,j) = val;
            }
        }
    }

    getMeshBorders(pMesh);
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if (checkBorder(rMesh, i, j)) {
                rMesh(i,j) = fiveDotScheme(pMesh, i, j) - func(pMesh.getPoint(i,j));
            }
        }
    }
    getMeshBorders(rMesh);
    double anumerator = 0, adenumerator = 0;
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(rMesh, i, j)) {
                anumerator += fiveDotScheme(rMesh, i, j) * gMesh(i,j);
                adenumerator += fiveDotScheme(gMesh, i, j) * gMesh(i,j);
            }
        }
    }
    double allAnumerator, allAdenumerator;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&anumerator, &allAnumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&adenumerator, &allAdenumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double alpha = allAnumerator / allAdenumerator;

    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(gMesh, i, j)) {
                gMesh(i,j) = rMesh(i,j) - alpha*gMesh(i,j);
            }
        }
    }

    getMeshBorders(gMesh);
    double tnumerator = 0, tdenumerator = 0;
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < cols; ++j) {
            if(checkBorder(gMesh, i, j)) {
                tnumerator += gMesh(i,j) * rMesh(i,j);
                tdenumerator += fiveDotScheme(gMesh, i, j) * gMesh(i, j);
            }
        }
    }

    double allTnumerator, allTdenumerator;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(&tnumerator, &allTnumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&tdenumerator, &allTdenumerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    tau = allTnumerator/ allTdenumerator;
    itercount++;
    int total_size;
    MPI_Comm_size(MPI_COMM_WORLD, &total_size);
    double total_error = 0;
    MPI_Allreduce(&err, &total_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(total_error);
}
double CGMPuassonIterator::iterate(Mesh &pMesh) {
    if(itercount == 0) {
        itercount++;
        return zeroIteration(pMesh);
    }
    long rows = pMesh.getRows();
    long cols = pMesh.getColumns();
    double err;
    for (std::size_t i = 1; i < rows - 1; ++i) {
        for (std::size_t j = 1; j < cols - 1; ++j) {
            double val = pMesh(i,j) - tau*gMesh(i,j);
            err += (std::fabs(pMesh(i,j) - val));
            pMesh(i,j) = val;
        }
    }
    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {
            rMesh(i,j) = fiveDotScheme(pMesh, i, j) - func(pMesh.getPoint(i,j));
        }
    }
    double anumerator = 0, adenumerator = 0;

    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {

            anumerator += fiveDotScheme(rMesh, i, j) * gMesh(i,j);
            adenumerator += fiveDotScheme(gMesh, i, j) * gMesh(i,j);
        }
    }
    double alpha = anumerator / adenumerator;

    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {
            gMesh(i,j) = rMesh(i,j) - alpha*gMesh(i,j);
        }
    }
    double tnumerator = 0, tdenumerator = 0;
    for (long i = 1; i < rows - 1; ++i) {
        for (long j = 1; j < cols - 1; ++j) {
            tnumerator += gMesh(i,j) * rMesh(i,j);
            tdenumerator += fiveDotScheme(gMesh, i, j) * gMesh(i, j);
        }
    }
    tau = tnumerator / tdenumerator;
    itercount++;
    return err;
}


