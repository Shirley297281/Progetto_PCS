#ifndef Utils_H
#define Utils_H
#include "Eigen/Eigen"
#include <algorithm> // qui si trova std::max_element e find


using namespace std;
using namespace Eigen;


// Funzione che calcola il lato di lunghezza massimo
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices);

// Funzione che calcola il baricentro
array <double,3> barycenter(const MatrixXd& m, unsigned int num_vertices);

// Funzione che controlla se c'è possibilità che due fratture si intersechino
bool check_sphere(const array<double,3> bar1, const array<double,3> bar2, const double l1, const double l2);

// Funzione che calcola il vettore normale al piano dove giace la frattura
Vector3d normal_vector(const MatrixXd& m);

// Funzione che calcola intersezione tra piani
bool intersezione_piani(Vector3d& n1,
                     Vector3d& n2,
                     array <double, 3>& b1,
                     array <double, 3>& b2,
                     Vector3d& t,
                     Vector3d& Point);

// Funzione che calcola intersezione tra rette
bool intersezione_rette (Vector3d& t,
                          Vector3d& V1,
                          Vector3d& V2,
                          Vector3d& Point,
                          Vector3d& Punto0);


// Funzione creata a partire dall'algoritmo di ordinamento BubbleSort la quale dato un vettore di array
// ordina il vettore in base alla seconda componente dell'array
template<typename T>
void BubbleSort_mod(vector<array<T,2>>& data)
{
    size_t r_size = data.size();
    size_t prec = r_size;
    bool swapped = true;

    while (swapped) {
        swapped = false;
        for (size_t i = 1; i < r_size; i++) {
            if (data[i-1][1] > data[i][1]) {
                swap(data[i-1], data[i]);
                swapped = true;
                prec = i;
            }
        }

        r_size = prec;
    }
}



#endif
