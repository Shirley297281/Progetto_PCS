#ifndef Namespace_H
#define Namespace_H
#include "Eigen/Eigen"
#include "FracturesTraces.hpp"

using namespace std;
using namespace Eigen;


namespace GeometryLibrary {


bool ImportFR(const string &filename,
              Fractures& fracture);

void CalcoloTracce(Fractures& fracture, Traces& trace);

int Controllo_tracce2(Fractures& fracture, Traces& trace, const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ,
                      const Vector3d& Point, Vector3d& t, unsigned int i, unsigned int j);

bool Tips_Shy(Fractures& fracture, Traces& trace, const int i, const int j );


//calculate free parameter of points
unsigned int Calcolo_par(Vector3d& t, Vector3d& Point,  int i, vector<Vector3d>& vec, Fractures& fracture);


bool exportFR1(const string &filename, const Traces& trace);

bool secondoOutput(const string &filename, const Fractures& fracture, Traces& trace);


}

#endif
