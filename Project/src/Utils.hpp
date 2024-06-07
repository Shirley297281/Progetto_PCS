#ifndef Utils_H
#define Utils_H
#include "Eigen/Eigen"
#include "inline.hpp"
#include <algorithm> // qui si trova std::max_element


using namespace std;
using namespace Eigen;


extern double tolDefault;

// compute the max eedge in fracture
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices);

// compute the barycenter of a fracture with "num_vertices" vertices
array <double,3> barycenter(const MatrixXd& m, unsigned int num_vertices);

// rough intersection control of fractures
bool check_sphere(const array<double,3> bar1, const array<double,3> bar2, const double l1, const double l2);

// calculate the normal vector of a fracture
Vector3d normal_vector(const MatrixXd& m);

// function to implement solution of systems 3x3 to find planes intersections
bool system_solution (Vector3d& n1,
                     Vector3d& n2,
                     array <double, 3>& b1,
                     array <double, 3>& b2,
                     Vector3d& t,
                     Vector3d& Point);

// troviamo il punto di intersezione tra la retta (Point, t)  e la retta di prolungamento del segmento V1V2
bool soluzione_sistema3x2 (Vector3d& t,
                          Vector3d& V1,
                          Vector3d& V2,
                          Vector3d& Point,
                          Vector3d& Punto0);


template<typename T>
void BubbleSort_mod(vector<array<T,2>>& data)
{
    size_t rem_size = data.size();
    size_t last_seen = rem_size;
    bool swapped = true;

    while (swapped) {
        swapped = false;
        for (size_t i = 1; i < rem_size; i++) {
            if (data[i-1][1] > data[i][1]) {
                swap(data[i-1], data[i]);
                swapped = true;
                last_seen = i;
            }
        }

        rem_size = last_seen;
    }
}



#endif
