#include "src/Utils.hpp"
#include "src/FracturesTraces.hpp"
#include "src/utils_partTwo.hpp"
#include <sstream>
#include <iostream>
#include "src/namespace.hpp"


using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;

int main()
{
    Fractures fracture;
    string filename = "FR82_data.txt";

    if( !ImportFR(filename, fracture) )
        return 1;

    // richiamo funzione calcola tracce
    Traces trace;
    CalcoloTracce(fracture, trace);

    cout<<"Tolleranza di defalut: "<<tolDefault<<endl;


    if (!exportFR1(filename, trace))
    {
        cerr<<"errore nella scrittura del file1";
    }

    if (!secondoOutput(filename, fracture, trace))
    {
        cerr<<"errore nella scrittura del file2";
    }

    ///Seconda Parte
    /// riferimento a utils_partTwo
    ///
    // DA QUI
    // scelta da utente l'id della frattura di cui voglio vedere la sottopoligonazione
    Polygons sottoPoligono;
    unsigned int z = 0;
    // INIZIO SALVATAGGIO PUNTI
    // (Parte di shyrl)
    Memo_Cello0Ds(fracture, trace, sottoPoligono, z);
    sottoPoligono.NumberCell0D = 0;
    // inserimento in Cell0Ds dei vertici della frattura Id (numeri naturali): in Cell0DId
    // inserimento in delle loro coordinate: in Cell0DCoordinates
    // incrementare NumberCell0D
    unsigned int numVerticiFrattZ = fracture.numVertices[z];
    MatrixXd insiemeVerticiFrattZ = fracture.CoordinatesVertice[z];
    // ciclo su fratture (cambiare TraceIdsPassxFracture con quello ordinato)
    unsigned int numTraccePassantiInZ = trace.TraceIdsPassxFracture[z].size();
    Vector3d vecNormaleAfratt = fracture.vettoreNormalePiano[z]; //mi serve dopo



    return 0;
}

