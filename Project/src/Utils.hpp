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

// inline function to implement solution of systems 3x3 to find planes intersections
inline bool system_solution (array <double, 3>& n1,
                            array <double,3>& n2,
                            array <double, 3>& b1,
                            array <double, 3>& b2,
                            array <double, 3>& t,
                            Vector3d& Point){
    double d1 = 0.0; //termine noto piano 1
    double d2 = 0.0; //termine noto piano 2
    //termine noto piano d'intersezione = 0 come prodotto scalare tra vettore normale
    // al piano e punto generico della frattura (abbiamo scelto il baricentro)
    d1 = n1[0] * b1[0] + n1[1] * b1[1]  + n1[2] * b1[2];
    d2 = n2[0] * b2[0] + n2[1] * b2[1]  + n2[2] * b2[2];

    // Convert arrays to Eigen vectors
    Vector3d n1_eigen(n1.data());
    Vector3d n2_eigen(n2.data());
    Vector3d t_eigen(t.data());

    // Create the matrix A and initialize it with vectors n1, n2, and t
    Matrix3d A;
    A.row(0) = n1_eigen;
    A.row(1) = n2_eigen;
    A.row(2) = t_eigen;

    //risolviamo il sistema
    // det(A) != 0 : per verificare che il istema abbia soluzione = i due piani si intersecano
    // aggiungere non complanarity! controllo che il modulo di t sia diverso da 0 per evitare prodotto vettoriale != 0

    // Verifica della possibilità di intersezione
    if (A.determinant() != 0 && t_eigen.dot(t_eigen) != 0) {
        // Risoluzione del sistema lineare per trovare il punto di intersezione
        Vector3d B;
        B << d1, d2, 0;
        Point = A.colPivHouseholderQr().solve(B);
        return true;

    } else {
        // I piani sono paralleli o coincidenti
        if (A.determinant() == 0) {
            if (abs(d1 - d2) < 1e-6) {
                // I piani sono coincidenti
                cout << "I piani sono coincidenti." << endl;
                return false;
            } else {
                // I piani sono paralleli ma non coincidenti
                cout << "I piani sono paralleli.";
                return false;
            }
        } else {
            // I piani sono complanari ma non paralleli
            cout << "I piani sono complanari." << endl;
            return false;
        }
    }

}

//troviamo il punto di intersezione tra la retta (Point, t)  e la retta di prolungamento del segmento V1V2
inline bool soluzione_sistema3x2 (array <double, 3>& t,
                                     Vector3d& V1,
                                     Vector3d& V2,
                                     Vector3d& Point,
                                     Vector3d& Punto0)
{

    Vector3d t_eigen(t.data());
    Vector3d vettoreDirezioneI = V1 - V2;
    VectorXd termine_noto = V1 - Point;

    MatrixXd A(3,2);
    A.col(0) << t_eigen;
    A.col(1) << vettoreDirezioneI;

    // Calcola il rango utilizzando la decomposizione LU
    int rank = A.fullPivLu().rank();
    //cout << "\n\nIl rango della matrice A è: " << rank << std::endl;
    if (rank == 2) {
        VectorXd sol;
        if (A.rows() >= A.cols()) {

            /*std::chrono::steady_clock::time_point t_begin= chrono::steady_clock::now();
            sol = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(termine_noto);
            std::chrono::steady_clock::time_point t_end= chrono::steady_clock::now();
            double duration = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();
            cout<<"\tjacobisvd.solve time: "<<duration<<" microseconds\n" <<endl;*/

            //chrono::steady_clock::time_point t_begin= chrono::steady_clock::now();
            sol = (A.transpose() * A).ldlt().solve(A.transpose() * termine_noto);
            //chrono::steady_clock::time_point t_end= chrono::steady_clock::now();
            //double duration = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();
            //cout<<"\tother time: "<<duration<<" microseconds\n" <<endl;

            cout << "Soluzione:\n" << sol << std::endl;
        } else {
            std::cout << "Impossibile risolvere: il sistema è sottodeterminato." << std::endl;
        }
    }

    return true;

}



#endif
