#ifndef Utils_2
#define Utils_2
#include "Eigen/Eigen"
#include "FracturesTracesPolygons.hpp"

using namespace std;
using namespace Eigen;

namespace Vettore{
    bool operator==(const VectorXd &v1, const VectorXd &v2);
}


namespace GeometryLibrary {


void MemorizzaVerticiPassanti_Cell0Ds(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z);
void MemorizzaVerticiNonPassanti_Cell0Ds (const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z);


void Creazioni_Sequenze_Passanti(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z);


void Creo_sottopoligono(unsigned int num_fracture, unsigned int num_sottopoligono, list<unsigned int> listaIdVertici, Polygons& poligoni, Fractures& fracture);




}




#endif
