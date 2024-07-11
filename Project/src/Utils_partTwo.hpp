#ifndef Utils_2
#define Utils_2
#include "Eigen/Eigen"
#include "FracturesTracesPolygons.hpp"

using namespace std;
using namespace Eigen;

// Estensione operatore di uguaglianza per vettori di Eigen
namespace Vettore{
    bool operator==(const VectorXd &v1, const VectorXd &v2);
}

namespace GeometryLibrary {

// Funzione che memorizza i punti appartenenti alla Mesh dovuti ai vertici della frattura di partenza, intersezioni lato-traccia passante e intersezioni tra tracce passanti
void MemorizzaVerticiPassanti_Cell0Ds(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z);

// Funzione che memorizza i punti appartenenti alla Mesh dovuti a intersezioni tracce non passanti - passanti e non passanti - non passanti già iterate
void MemorizzaVerticiNonPassanti_Cell0Ds (const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z, vector<Matrix<double, 3, 4>>& NuoviEstremi);

// Funzione che assegna delle sequenze ai punti della mesh in base alla loro posizione rispetto alle tracce passanti
void Creazioni_Sequenze_Passanti(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z);

// Funzione che assegna delle sequenze ai punti della mesh in base alla loro posizione rispetto alle tracce non passanti
void Creazioni_Sequenze_NONPassanti(const Fractures& fracture, Polygons& sottoPoligono, unsigned int z, vector<Matrix<double, 3, 4>>& NuoviEstremi);

// Funzione che dati gli insiemi di punti di un sottopoligono li ordina in ordine antiorario e aggiorna Cell2D
void Creo_sottopoligono(unsigned int num_fracture, unsigned int num_sottopoligono, list<unsigned int> listaIdVertici, Polygons& poligoni, Fractures& fracture);

}




#endif
