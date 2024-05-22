#ifndef Utils_H
#define Utils_H
#include "Eigen/Eigen"
#include <iostream>
#include <algorithm> // qui si trova std::max_element
#include "FracturesTraces.hpp"

using namespace std;
using namespace Eigen;

namespace GeometryLibrary {


bool ImportFR(const string &filename,
                   Fractures& fracture);

void CalcoloTracce(Fractures& fracture, Traces& trace);

bool Tips_Shy(Fractures fracture, Traces trace,const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ, const int i, const int j ,map<double, Vector3d>& dizfreeParToVec,
              double freeParP1,
              double freeParP2,
              double freeParP3,
              double freeParP4);

int Controllo_tracce2(Fractures fracture, Traces trace, const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ,
                      const Vector3d Point, Vector3d t, unsigned int i, unsigned int j);

//calculate free parameter of points
unsigned int Calcolo_par(Vector3d& t, Vector3d& Point,  int i, vector<Vector3d>& vecI, Fractures& fracture);

void inserimento_map(double& idpar1, double& idpar2, Traces& trace);

}




// compute the max eedge in fracture
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices);

// compute the barycenter of a fracture with "num_vertices" vertices
array <double,3> barycenter(const MatrixXd& m, unsigned int num_vertices);

// rough intersection control of fractures
bool check_sphere(const array<double,3> bar1, const array<double,3> bar2, const double l1, const double l2);

// calculate the normal vector of a fracture
Vector3d normal_vector(const MatrixXd& m);





// inline function to implement vectorial product
inline Vector3d vec_product(Vector3d& v1, Vector3d& v2){

    Vector3d v = {};
    v[0] = v1[1]*v2[2] - v2[1]*v1[2];
    v[1] = - v1[0]*v2[2] + v1[2]*v2[0];
    v[2] = v1[0]*v2[1] - v2[0]*v1[1];

    return v;
}

// inline function to implement solution of systems 3x3 to find planes intersections
inline bool system_solution (Vector3d& n1,
                            Vector3d& n2,
                            array <double, 3>& b1,
                            array <double, 3>& b2,
                            Vector3d& t,
                            Vector3d& Point){
    double d1 = 0.0; //termine noto piano 1
    double d2 = 0.0; //termine noto piano 2
    //termine noto piano d'intersezione = 0 come prodotto scalare tra vettore normale
    // al piano e punto generico della frattura (abbiamo scelto il baricentro)
    d1 = n1[0] * b1[0] + n1[1] * b1[1]  + n1[2] * b1[2];
    d2 = n2[0] * b2[0] + n2[1] * b2[1]  + n2[2] * b2[2];


    // Create the matrix A and initialize it with vectors n1, n2, and t
    Matrix3d A;
    A.row(0) = n1;
    A.row(1) = n2;
    A.row(2) = t;

    //risolviamo il sistema
    // det(A) != 0 : per verificare che il istema abbia soluzione = i due piani si intersecano
    // aggiungere non complanarity! controllo che il modulo di t sia diverso da 0 per evitare prodotto vettoriale != 0

    // Verifica della possibilità di intersezione
    if (abs(A.determinant()) > 1e-14 && t.dot(t) > 1e-14) {
        // Risoluzione del sistema lineare per trovare il punto di intersezione
        Vector3d B;
        B << d1, d2, 0;
        Point = A.colPivHouseholderQr().solve(B);
        return true;

    } else {
        // I piani sono paralleli o coincidenti
        return false;
    }

}

//troviamo il punto di intersezione tra la retta (Point, t)  e la retta di prolungamento del segmento V1V2
inline bool soluzione_sistema3x2 (Vector3d& t,
                                     Vector3d& V1,
                                     Vector3d& V2,
                                     Vector3d& Point,
                                     Vector3d& Punto0)
{

    Vector3d vettoreDirezioneI = V1 - V2;
    VectorXd termine_noto = V1 - Point;

    MatrixXd A(3,2);
    A.col(0) << t;
    A.col(1) << vettoreDirezioneI;

    MatrixXd Completa(3,3); //è A con ultima colonna termine noto
    Completa.col(0) << t;
    Completa.col(1) << vettoreDirezioneI;
    Completa.col(2) << termine_noto;

    // Calcola il rango utilizzando la decomposizione LU
    int rankA = A.fullPivLu().rank();

    int rankC = Completa.fullPivLu().rank();


    //cout << "\n\nIl rango della matrice A è: " << rank << std::endl;
    if (rankA == rankC && rankA == 2) {
        VectorXd sol;
        if (A.rows() >= A.cols()) {

            /*std::chrono::steady_clock::time_point t_begin= chrono::steady_clock::now();
            sol = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(termine_noto);
            std::chrono::steady_clock::time_point t_end= chrono::steady_clock::now();
            double duration = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();
            cout<<"\tjacobisvd.solve time: "<<duration<<" microseconds\n" <<endl;
            */

            //chrono::steady_clock::time_point t_begin2= chrono::steady_clock::now();
            sol = (A.transpose() * A).ldlt().solve(A.transpose() * termine_noto);
            //chrono::steady_clock::time_point t_end2= chrono::steady_clock::now();
            //double duration2 = chrono::duration_cast<chrono::microseconds>(t_end2 - t_begin2).count();
            //cout<<"\tother time: "<<duration2<<" microseconds\n" <<endl;

            //ora abbiamo trovato i coeffieiente s e aplha

            // Calcola il punto di intersezione Punto0
            Punto0 = Point + sol[0] * t ;
            return true;
        }
    }
    return false;



}

// Funzione per calcolare la distanza euclidea senza usare .norm()
inline double euclidean_distance(const Vector3d& a, const Vector3d& b) {
    return sqrt((a(0) - b(0)) * (a(0) - b(0)) +
                (a(1) - b(1)) * (a(1) - b(1)) +
                (a(2) - b(2)) * (a(2) - b(2)));
}


inline void inserimento_map(double pass, unsigned int idpar, GeometryLibrary::Traces& trace) {
    if (pass==0) {

        // Inserisci nella mappa delle tracce passanti
        auto ret = trace.TraceIdsPassxFracture.insert({idpar, {(trace.numTraces - 1)}});
        if (!ret.second) {
            ret.first->second.push_back(trace.numTraces - 1);
        }
    } else {
        // Inserisci nella mappa delle tracce non passanti
        auto ret = trace.TraceIdsNoPassxFracture.insert({idpar, {(trace.numTraces - 1)}});
        if (!ret.second) {
            ret.first->second.push_back(trace.numTraces - 1);
        }
    }
}




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
