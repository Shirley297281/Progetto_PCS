#ifndef Utils_H
#define Utils_H
#include "Eigen/Eigen"
#include <iostream>
#include "FracturesTraces.hpp"

using namespace std;
using namespace Eigen;

namespace GeometryLibrary {

bool ImportMesh(Fractures& mesh);


bool ImportFR(const string &filename,
                   Fractures& fracture);


}

// compute euclidean distance between two points 3-dimensional
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices);

// compute the barycenter of a fracture with "num_vertices" vertices
vector <double> barycenter(const MatrixXd& m, unsigned int num_vertices);

#endif
