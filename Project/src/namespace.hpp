#ifndef Namespace_H
#define Namespace_H
#include "Eigen/Eigen"
#include "FracturesTracesPolygons.hpp"

using namespace std;
using namespace Eigen;


namespace GeometryLibrary {

// Funzione per l'importazione dei dati sulle fratture da file assegnati
bool ImportFR(const string &filename,
              Fractures& fracture);

// Funzione per l'individuazione di tracce
void CalcoloTracce(Fractures& fracture, Traces& trace);

// Funzione per discriminare tracce passanti e non
int distinzioneTipoTraccia1(Traces& trace, const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ,
                      const Vector3d& Point, Vector3d& t, unsigned int i, unsigned int j);

// Funzione per discriminare tracce passanti e non
void Tips_Shy(Fractures& fracture, Traces& trace, const int i, const int j );

// Funzione sfruttata in Calcolo tracce per il controllo dell'intersezione tra traccia e frattura
unsigned int Calcolo_par(Vector3d& t, Vector3d& Point,  int i, vector<Vector3d>& vec, Fractures& fracture);

// Funzione per esportare i dati riguardo le tracce nel primo file di export
bool exportFR1(const string &filename, const Traces& trace);

// Funzione per esportare i dati riguardo le tracce nel secondo file di export
// e per ogni fratture ordina le tracce passanti e non in ordine di lunghezza
bool secondoOutput(const string &filename, const Fractures& fracture, Traces& trace);


}

#endif
