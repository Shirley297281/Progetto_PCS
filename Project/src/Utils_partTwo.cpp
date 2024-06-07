#include "Utils.hpp"
#include "FracturesTraces.hpp"
#include "inline.hpp"
#include <vector>
#include "Eigen/Eigen"

#include <vector>


using namespace std;
using namespace Eigen;


namespace GeometryLibrary{

void Memo_Cello0Ds(const Fractures& fracture, const Traces& trace, Polygons& poligono, unsigned int z){
    //aggiornamento temporaneo del numero di Punti utili
    int numVertices_forthisfracture=fracture.numVertices[z];
    //vale solo per le passanti:
    //poligono.NumberCell0D = numVertices_forthisfracture + trace.NumTraccePerFrattura[z][0] * 2; per questo un controllo serve sull'uguag coord

    for(int i = 0; i<numVertices_forthisfracture; i++){
        poligono.Cell0DCoordinates.push_back(fracture.CoordinatesVertice[z].col(i));
        poligono.NumberCell0D++;
    }

    //crush
    /*//for sugli id delle tracce passanti
    for (const auto& idt : trace.TraceIdsPassxFracture[z] ){
        //for che itera sulle coordinate già inserite
        for (int i = 0; i < poligono.Cell0DCoordinates.size(); ++i) {
            //controllo il primo estremo della traccia
            if (poligono.Cell0DCoordinates[i] == trace.CoordinatesEstremiTraces[idt].col(0)) {
                //vertice già presente in celle0d
            }else{
                //aggiungere il vertice al vector
                poligono.Cell0DCoordinates.push_back(trace.CoordinatesEstremiTraces[idt].col(0));
                //e incremento il numero di Cell0D
                poligono.NumberCell0D ++;
            }

            //controllo il secondo estremo della traccia
            if (poligono.Cell0DCoordinates[i] == trace.CoordinatesEstremiTraces[idt].col(1)) {
                //vertice già presente in celle0d
            }else{
                //aggiungere il vertice al vector
                poligono.Cell0DCoordinates.push_back(trace.CoordinatesEstremiTraces[idt].col(1));
                //e incremento il numero di Cell0D
                poligono.NumberCell0D ++;
            }
        }
    }*/



    // // inserimento in delle loro coordinate: in Cell0DCoordinates
    // // incrementare NumberCell0D
    // unsigned int numVerticiFrattZ = fracture.numVertices[z];
    // MatrixXd insiemeVerticiFrattZ = fracture.CoordinatesVertice[z];
    // // ciclo su fratture (cambiare TraceIdsPassxFracture con quello ordinato)
    // unsigned int numTraccePassantiInZ = trace.TraceIdsPassxFracture[z].size();
    // Vector3d vecNormaleAfratt = fracture.vettoreNormalePiano[z]; //mi serve dopo

}















}

