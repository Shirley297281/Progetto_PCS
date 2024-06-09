#include "Utils.hpp"
#include "FracturesTraces.hpp"
#include "inline.hpp"
#include <vector>
#include "Eigen/Eigen"

#include <vector>


using namespace std;
using namespace Eigen;

namespace Vettore{
bool operator==(const VectorXd &v1, const VectorXd &v2) // operatore per uguaglianza tra vettori in std
{
    if(v1.size() == 0 || v1.size()!= v2.size())
        return false;
    for (unsigned int i = 0; i < v1.size(); i++) {
        if (v1[i] != v2[i]) {
            return false;
        }
    }
    return true;
}
}



namespace GeometryLibrary{

void Memo_Cello0Ds(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z){
    // salvo i punti coincidenti con i vertici della frattura
    MatrixXd insiemeVerticiFrattZ = fracture.CoordinatesVertice[z];//estraggo la matrice di vertici di z
    unsigned int numVerticiFrattZ = fracture.numVertices[z];
    map<unsigned int, list<unsigned int>>& markerDiz = sottoPoligono.Cell1DMarkers; // passo in referenza così lo modifico direttamente senza farne una copia
    for (unsigned int i = 0; i<numVerticiFrattZ; i++) //ciclo sui vertici
    {
        Vector3d vertice = insiemeVerticiFrattZ.col(i);
        unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
        sottoPoligono.NumberCell0D = NumPuntiFinora  + 1;
        // visto che push_back alloca memoria e inserisce dovrebbe funzionare comunque (reserve non sapevo come farlo)
        sottoPoligono.Cell0DId.push_back(NumPuntiFinora); //d'ora in poi userò come id dei vertici il numero di punti salvati
        sottoPoligono.Cell0DCoordinates.push_back(vertice);
        MatrixXd M(0,0);
        sottoPoligono.SequenzeXpunto.push_back(M);
        markerDiz[0].push_back(NumPuntiFinora); // marker con chiave 0 se vertici
    }
    // salvo punti che derivano da intersezione lato-traccia
    /// ciclo su fratture (cambiare TraceIdsPassxFracture con quello ordinato)
    unsigned int numTraccePassantiInZ = trace.TraceIdsPassxFracture[z].size();
    Vector3d vecNormaleAfratt = fracture.vettoreNormalePiano[z]; //mi serve dopo
    for (unsigned int i = 0; i<numTraccePassantiInZ ; i++ )
    {
        unsigned int idTraccia = trace.TraceIdsPassxFracture[z][i];
        Vector3d Estremo1Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(0); //estraggo estremo 1 della traccia i
        Vector3d Estremo2Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(1); //estraggo estremo 2 della traccia i
        Vector3d t = Estremo2Traccia - Estremo1Traccia; //serve dopo
        // ciclo su lati della frattura
        for (unsigned int j = 0; j < numVerticiFrattZ ; j++)
        {
            // estraggo un vertice della frattura e il successivo per individuare un lato della frattura
            Vector3d Vertice1 = insiemeVerticiFrattZ.col(j);
            Vector3d Vertice2 = {0,0,0};
            if (j == numVerticiFrattZ -1)
            {
                Vector3d Vertice2 = insiemeVerticiFrattZ.col(0);
            }
            else
            {
                Vector3d Vertice2 = insiemeVerticiFrattZ.col(j+1);
            }

            // calcolo intersezione retta su cui giace il lato e retta su cui giace la traccia
            Vector3d Punto0 = {0, 0, 0};
            if (!soluzione_sistema3x2(t, Vertice1,Vertice2,Estremo1Traccia, Punto0))
            {
                continue;
            }
            // controllo se il punto è compreso tra gli estremi della traccia in modo da
            // individuare solo punti all'interno della frattura (infatti essendo queste passanti avrenno estremi su lati)
            // inline per combinazione convessa
            if (!combinazione_convessa(Estremo1Traccia, Estremo2Traccia, Punto0))
            {
                continue;
            }
            // controllo che punto non sia già stato inserito nelle strutture (in caso di intersezioni coincidenti)
            vector<Vector3d> VettoreCoordinateIn0D = sottoPoligono.Cell0DCoordinates;
            // CheckInserimento returna true se l'inserimento non è ancora avvenuto e false se è già avvenuto (inline fun)
            if (!checkInserimento(Punto0, VettoreCoordinateIn0D)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
            {
                continue;
            }
            unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
            sottoPoligono.NumberCell0D = NumPuntiFinora  + 1;
            /// visto che push_back alloca memoria e inserisce dovrebbe funzionare comunque (reserve non sapevo come farlo)
            sottoPoligono.Cell0DId.push_back(NumPuntiFinora); // es: se ho già 2 punti questi hanno identificativo 0,1. Quando aggiungo il terzo questo avrà id = 2
            sottoPoligono.Cell0DCoordinates.push_back(Punto0);
            // non so se sia necessario questo però sto allocando memoria
            MatrixXd M(0,0); //matrice vuota di dimensione
            sottoPoligono.SequenzeXpunto.push_back(M);
            markerDiz[1].push_back(NumPuntiFinora); //marker con chiave 1 per punti sui lati
        }
        // salvo punti che derivano da intersezione di due tracce
        // ciclo sulle fratture che rimangono (prendo sempre quelle successive ma prima controllo di non sforare con l'iteratore)
        if (i == trace.TraceIdsPassxFracture[z].size()-1)
        {
            continue;
        }
        for (unsigned int j = i + 1; i< trace.TraceIdsPassxFracture[z].size() ; j++ )
        {
            unsigned int idTraccia2 = trace.TraceIdsPassxFracture[z][j];
            Vector3d Estremo1Traccia2 = trace.CoordinatesEstremiTraces[idTraccia2].col(0); //estraggo estremo 1 della traccia j
            Vector3d Estremo2Traccia2 = trace.CoordinatesEstremiTraces[idTraccia2].col(1); //estraggo estremo 2 della traccia j
            // calcolo intersezione retta su cui giace il lato e retta su cui giace la traccia
            Vector3d Punto02 = {0, 0, 0};
            if (!soluzione_sistema3x2(t, Estremo1Traccia2, Estremo2Traccia2,Estremo1Traccia, Punto02))
            {
                continue;
            }
            // controllo se il punto è compreso tra gli estremi di almeno una delle due traccia (in teoria basta perchè i punti di una traccia per def appartengono alla frattura)
            // ho scelto senza un motivo (è uguale) la traccia i
            // inline per combinazione convessa
            if (!combinazione_convessa(Estremo1Traccia, Estremo2Traccia, Punto02))
            {
                continue;
            }
            // controllo che punto non sia già stato inserito nelle strutture (in caso di intersezioni coincidenti)
            vector<Vector3d> VettoreCoordinateIn0D = sottoPoligono.Cell0DCoordinates;
            // CheckInserimento returna true se l'inserimento non è ancora avvenuto e false se è già avvenuto (inline fun)
            if (!checkInserimento(Punto02, VettoreCoordinateIn0D)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
            {
                continue;
            }
            unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
            sottoPoligono.NumberCell0D = NumPuntiFinora  + 1;
            /// visto che push_back alloca memoria e inserisce dovrebbe funzionare comunque (reserve non sapevo come farlo)
            sottoPoligono.Cell0DId.push_back(NumPuntiFinora); // es: se ho già 2 punti questi hanno identificativo 0,1. Quando aggiungo il terzo questo avrà id = 2
            sottoPoligono.Cell0DCoordinates.push_back(Punto02);
            // non so se sia necessario questo però sto allocando memoria
            MatrixXd M(0,0); //matrice vuota di dimensione
            sottoPoligono.SequenzeXpunto.push_back(M); //pensare a un resize eventuale
            markerDiz[2].push_back(NumPuntiFinora); //marker con chiave 2 per punti interni
        }
    }

}















}

