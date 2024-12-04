#ifndef __FRACTURESTRACESPOLYGONS_H
#define __FRACTURESTRACESPOLYGONS_H

#include "Eigen/Eigen"
#include <array>
#include <vector>

using namespace Eigen;
using namespace std;

// Definizione della variabile globale
extern double tolDefault; //= 100 * std::numeric_limits<double>::epsilon();

namespace GeometryLibrary {

// 1^ PARTE
struct Fractures{
    unsigned int NumFractures = 0;
    vector<unsigned int> IdFractures ; //vettore con Id fratture
    vector<MatrixXd> CoordinatesVertice ; // matrice con coordinate dei vertici
    vector<double> lenghtMaxEdges = {} ; // vettore di lato con lunghezza massima per ogni frattura
    vector<array<double, 3>> baricentro = {}; // vettore baricentro per ogni frattura
    vector<Vector3d> vettoreNormalePiano = {}; // vettore di vettori normali al piano di giacenza per ogni frattura
    vector<unsigned int> numVertices = {}; //numero vertici per ogni frattura

    Fractures() = default;
    Fractures(const vector<unsigned int> IdFractures,
              const  vector<MatrixXd> CoordinatesVertice):
        IdFractures(IdFractures),
        CoordinatesVertice(CoordinatesVertice)
    {}
};


struct Traces{
    unsigned int numTraces = 0; // numero tracce totali
    vector<unsigned int> IdTraces = {};  // vettore con Id tracce
    vector<Matrix<double, 3, 2>> CoordinatesEstremiTraces = {}; // vettore con matrice((x,y,z), estremo1 x estremo2)
    vector<double> lengthTraces = {};  // vettore lunghezza tracce

    //Per ogni traccia memorizzo i due id delle fratture
    vector<array<unsigned int,2>> IdsFractures;

    // Per ogni frattura, elenco delle tracce passanti
    vector<vector<unsigned int>> TraceIdsPassxFracture;
    // Per ogni frattura, elenco delle tracce non passanti
    vector<vector<unsigned int>> TraceIdsNoPassxFracture;

    Traces() = default;
    Traces(const vector<vector<unsigned int>>& TraceIdsPassxFracture,
          const vector<vector<unsigned int>>& TraceIdsNoPassxFracture):
         TraceIdsPassxFracture( TraceIdsPassxFracture),
         TraceIdsNoPassxFracture(TraceIdsNoPassxFracture)
    {}

};

// 2^PARTE
struct Polygons{

    /// Cell0D
    unsigned int NumberCell0D = 0;
    // Id of each point in mesh
    vector<unsigned int> Cell0DId = {};
    // coordinates of each point
    vector<Vector3d> Cell0DCoordinates = {};

    // map for linking each marker (key) with the list of vertices Id associated with that marker
    map<unsigned int, list<unsigned int>> Cell0DMarkers = {};


    /// Cell1D
    unsigned int NumberCell1D = 0;
    vector<unsigned int> Cell1DId = {};
    vector<Vector2i> Cell1DVertices = {};  // per ogni lato coppia id vertici


    /// Cell2D
    unsigned int NumberCell2D= 0;
    vector<unsigned int> Cell2DId = {};
    list<unsigned int> NumberVertices = {};  // numero vertici per poligono
    vector<vector<unsigned int>> Cell2DVertices = {};       // Per ogni poligono gli id dei vertici
    list<unsigned int> NumberEdges  = {};   // numero lati per poligono

    vector<vector<unsigned int>> Cell2DEdges = {};  // Per ogni poligono gli id dei lati

    // strutture di supporto utilizzate nella parte 2
    vector<MatrixXd> SequenzeXpunto;
    vector<Vector3d> CoordinatesPunto;


    void GedimInterface(vector<vector<unsigned int>>& triangles,
                        VectorXi& materials);
};
}

#endif
