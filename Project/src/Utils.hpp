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
array <double,3> barycenter(const MatrixXd& m, unsigned int num_vertices);

// rough intersection control of fractures
bool check_sphere(const double bar1, const double bar2, const double l1, const double l2);

// calculate the normal vector of a fracture
array <double, 3> normal_vector(const MatrixXd& m);

// inline function to implement vectorial product
inline array <double,3> vec_product(array <double, 3>& v1, array <double, 3>& v2){

    array <double, 3> v = {};
    v[0] = v1[1]*v2[2] - v2[1]*v1[2];
    v[1] = - v1[0]*v2[2] + v1[2]*v2[0];
    v[2] = v1[0]*v2[1] - v2[0]*v1[1];

    return v;
}

// inline function to implement scalar product
inline Vector3d system_solution (array <double, 3>& n1,
                                        array <double,3>& n2,
                                        array <double, 3>& b1,
                                        array <double, 3>& b2,
                                        array <double, 3>& t){
    double d1 = 0.0; //termine noto piano 1
    double d2 = 0.0; //termine noto piano 2
     //termine noto piano d'intersezione = 0;
    d1 = n1[0] * b1[0] + n1[1] * b1[1]  + n1[2] * b1[2];
    d2 = n2[0] * b2[0] + n2[1] * b2[1]  + n2[2] * b2[2];

    // Convert arrays to Eigen vectors
    Vector3d n1_eigen(n1.data());
    Vector3d n2_eigen(n2.data());
    Vector3d t_eigen(t.data());

    // Create the matrix A and initialize it with vectors n1, n2, and t
    Matrix3d A;
    A << n1_eigen, n2_eigen, t_eigen;

    Vector3d P ={};

    //risolviamo il sistema
    if ( A.determinant() != 0 ){ //aggiungere non complanarity?!
        //risolvo il prodotto AP = d

        Vector3d B;
        B << d1, d2, 0;

        // Risoluzione del sistema lineare
        Vector3d P = A.colPivHouseholderQr().solve(B);

        // Stampa delle soluzioni
        cout << "Soluzioni per P:" << endl;
        cout << "P0: " << P(0) << endl;
        cout << "P1: " << P(1) << endl;
        cout << "P2: " << P(2) << endl;

    }else{
        //there is no intersection??
    }

    return P;
}

//vector <double> planes_intersection (const piano1, const piano2);

#endif
