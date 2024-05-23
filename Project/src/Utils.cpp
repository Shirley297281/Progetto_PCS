#include "Utils.hpp"
//#include "FracturesTraces.hpp"
#include "inline.hpp"
#include <vector>
#include "Eigen/Eigen"
#include <cmath> // per sqrt
#include <vector>


using namespace std;


double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices)
{


    vector<double> distance_vector(num_vertices, 0.0);
    // Itero sulle colonne
    for(unsigned int i = 0; i < num_vertices - 1; ++i){
        Vector3d point_a = m.col(i);
        Vector3d point_b = m.col(i + 1);
        distance_vector[i] = euclidean_distance(point_a, point_b);
    }

    // Calcolo la distanza tra l'ultimo punto e il primo
    Vector3d last_point = m.col(num_vertices - 1);
    Vector3d first_point = m.col(0);
    distance_vector[num_vertices - 1] = euclidean_distance(last_point, first_point);

    // Trovo la distanza massima
    auto it_max_distance = max_element(distance_vector.begin(), distance_vector.end());
    double max_distance = *it_max_distance;

    return max_distance;

}


//baricentro
array <double,3> barycenter(const MatrixXd& m, unsigned int num_vertices)
{

    array<double,3> barycenter_coord ={};

    for(unsigned int i = 0; i < 3; ++i){
        for(unsigned int j = 0; j < num_vertices ; ++j){
            barycenter_coord[i] += m(i,j);
        }
        barycenter_coord[i] /= num_vertices;
    }
    return barycenter_coord;
}

//funzione sfere
bool check_sphere(const array<double,3> bar1, const array<double,3> bar2, const double l1, const double l2)
{
    //controlliamo la distanza tra i due baricentri
    double distance_bar = 0.0;
    distance_bar = sqrt( (bar1[0] - bar2[0])*(bar1[0] - bar2[0]) +
                         (bar1[1] - bar2[1])*(bar1[1] - bar2[1]) +
                         (bar1[2] - bar2[2])*(bar1[2] - bar2[2]) );

    double max_distance = 0.0;
    max_distance = (l1 + l2) ;///h

    if (distance_bar > max_distance){
        return false; // le fratture non si intersecano
    }

    return true;
}

//normale al piano
Vector3d normal_vector(const MatrixXd& m)
{

    Vector3d v1 = {}; // vettore1
    Vector3d v2 = {}; // vettore2
    Vector3d v = {};

    v1[0] = m(0,1) - m(0,0);
    v1[1] = m(1,1) - m(1,0);
    v1[2] = m(2,1) - m(2,0);

    v2[0] = m(0,1) - m(0,2);
    v2[1] = m(1,1) - m(1,2);
    v2[2] = m(2,1) - m(2,2);

    // prodotto vettoriale
    v = vec_product(v1, v2);

    return v;

}




