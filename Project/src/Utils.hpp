#ifndef Utils_H
#define Utils_H
#include "Eigen/Eigen"
#include <iostream>
#include "FracturesTraces.hpp"

using namespace std;
using namespace Eigen;

namespace GeometryLibrary {


bool ImportFR(const string &filename,
                   Fractures& fracture);

void CalcoloTracce(Fractures& fracture, Traces& trace);

}

// compute euclidean distance between two points 3-dimensional
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices);

// compute the barycenter of a fracture with "num_vertices" vertices
vector <double> barycenter(const MatrixXd& m, unsigned int num_vertices);

// rough intersection control of fractures
bool check_sphere(const double bar1, const double bar2, const double l1, const double l2);

// calculate the normal vector of a fracture
array <double, 3> normal_vector(const MatrixXd& m);

//vector <double> planes_intersection (const piano1, const piano2);

#endif
