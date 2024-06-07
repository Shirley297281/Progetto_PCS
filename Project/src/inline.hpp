#ifndef Inline_H
#define Inline_H

#include "Eigen/Eigen"
#include "FracturesTraces.hpp"


using namespace std;
using namespace Eigen;

// inline function to implement vectorial product
inline Vector3d vec_product(Vector3d& v1, Vector3d& v2){

    Vector3d v = {};
    v[0] = v1[1]*v2[2] - v2[1]*v1[2];
    v[1] = - v1[0]*v2[2] + v1[2]*v2[0];
    v[2] = v1[0]*v2[1] - v2[0]*v1[1];

    return v;
}


// Funzione per calcolare la distanza euclidea senza usare .norm()
inline double euclidean_distance(const Vector3d& a, const Vector3d& b) {
    return sqrt((a(0) - b(0)) * (a(0) - b(0)) +
                (a(1) - b(1)) * (a(1) - b(1)) +
                (a(2) - b(2)) * (a(2) - b(2)));
}


inline void inserimento_map(int pass, unsigned int idpar, GeometryLibrary::Traces& trace) {
    if (pass == 0) {
        // Assicurati che il vettore interno esista
        if (idpar >= trace.TraceIdsPassxFracture.size()) {
            trace.TraceIdsPassxFracture.resize(idpar + 1);
        }
        // Inserisci l'ID della traccia nel vector TraceIdsPassxFracture nella posizione idpar
        trace.TraceIdsPassxFracture[idpar].push_back(trace.numTraces - 1);

    } else {
        // Assicurati che il vettore interno esista
        if (idpar >= trace.TraceIdsNoPassxFracture.size()) {
            trace.TraceIdsNoPassxFracture.resize(idpar + 1);
        }
        // Inserisci l'ID della traccia nel vector TraceIdsNoPassxFracture nella posizione idpar
        trace.TraceIdsNoPassxFracture[idpar].push_back(trace.numTraces - 1);
    }
}

// Funzione inline per calcolare che un punto sia combinazione convessa di altri due
inline bool combinazione_convessa(Vector3d v1, Vector3d v2, Vector3d p){
    // risolvo per alpha in (0,1)
    // p[0] = (1-alpha)*v1[0] + alpha*v2[0]
    // p[1] = (1-alpha)*v1[1] + alpha*v2[1]
    // p[2] = (1-alpha)*v1[2] + alpha*v2[2]

    double alpha_x = (p[0] - v1[0])/(v2[0] - v1[0]);
    double alpha_y = (p[1] - v1[1])/(v2[1] - v1[1]);
    double alpha_z = (p[2] - v1[2])/(v2[2] - v1[2]);

    // verifico che i valori di alhpa calcolati nelle coordinate x, y, z siano approssivativamente uguali e che alpha sia in (0,1)
    if(abs(alpha_x - alpha_y) < 1e-6 && abs(alpha_y - alpha_z) < 1e-6 && alpha_x >= 0 && alpha_x <= 1){
        return true;
    }
    return false;
}

#endif
