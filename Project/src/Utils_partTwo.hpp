#ifndef Utils_2
#define Utils_2
#include "Eigen/Eigen"
#include "FracturesTraces.hpp"

using namespace std;
using namespace Eigen;

namespace Vettore{
bool operator==(const VectorXd &v1, const VectorXd &v2);
}


namespace GeometryLibrary {


void MemorizzaVertici_Cell0Ds(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z);


void Creazioni_Sequenze(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z);







}




#endif
