#ifndef Utils_H
#define Utils_H
#include "Eigen/Eigen"
#include <algorithm> // qui si trova std::max_element


using namespace std;
using namespace Eigen;




// compute the max eedge in fracture
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices);

// compute the barycenter of a fracture with "num_vertices" vertices
array <double,3> barycenter(const MatrixXd& m, unsigned int num_vertices);

// rough intersection control of fractures
bool check_sphere(const array<double,3> bar1, const array<double,3> bar2, const double l1, const double l2);

// calculate the normal vector of a fracture
Vector3d normal_vector(const MatrixXd& m);


template<typename T>
void BubbleSort(vector<array<T,2>>& data)
{
    size_t rem_size = data.size();
    size_t last_seen = rem_size;
    bool swapped = true;

    while (swapped) {
        swapped = false;
        for (size_t i = 1; i < rem_size; i++) {
            if (data[i-1][1] > data[i][1]) {
                std::swap(data[i-1], data[i]);
                swapped = true;
                last_seen = i;
            }
        }

        rem_size = last_seen;
    }
}



#endif
