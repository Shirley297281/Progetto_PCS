#ifndef Namespace_H
#define Namespace_H
#include "Eigen/Eigen"
#include "FracturesTracesPolygons.hpp"

using namespace std;
using namespace Eigen;


namespace GeometryLibrary {


bool ImportFR(const string &filename,
              Fractures& fracture);

void CalcoloTracce(Fractures& fracture, Traces& trace);

int distinzioneTipoTraccia1(Traces& trace, const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ,
                      const Vector3d& Point, Vector3d& t, unsigned int i, unsigned int j);

void Tips_Shy(Fractures& fracture, Traces& trace, const int i, const int j );

unsigned int Calcolo_par(Vector3d& t, Vector3d& Point,  int i, vector<Vector3d>& vec, Fractures& fracture);


bool exportFR1(const string &filename, const Traces& trace);

bool secondoOutput(const string &filename, const Fractures& fracture, Traces& trace);


}

#endif
