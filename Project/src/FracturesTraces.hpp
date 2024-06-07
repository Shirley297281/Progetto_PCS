#ifndef __FRACTURESTRACES_H
#define __FRACTURESTRACES_H

#include "Eigen/Eigen"
#include <array>
#include <vector>

using namespace Eigen;
using namespace std;


namespace GeometryLibrary {

// 1^ PARTE
struct Fractures{
    unsigned int NumFractures = 0;
    vector<unsigned int> IdFractures ; //vettore con Id fratture
    vector<MatrixXd> CoordinatesVertice ; // matrice con coordinate dei vertici (ogni colonna individua un vertice)
    vector<double> lenghtMaxEdges = {} ; // vettore di lunghezza massima per ogni frattura
    vector<array<double, 3>> baricentro = {}; // vettore baricentro per ogni frattura
    vector<Vector3d> vettoreNormalePiano = {}; // vettore di versori per ogni frattura (i,j,k)
    vector<unsigned int> numVertices = {}; //numero vertici per ogni poligono

    // (list o vector dipende da cosa ci serve dopo)
    Fractures() = default; //DA CHIEDERE
    Fractures(const vector<unsigned int> IdFractures,
              const  vector<MatrixXd> CoordinatesVertice):
        IdFractures(IdFractures),
        CoordinatesVertice(CoordinatesVertice)
    {}
    // void importPolygonsList(const string& filepath, Polygons& polygons); NEL UTILS
    // void GedimInterface(vector<vector<unsigned int>>& triangles, VectorXi& materials); PER PARAVIEW
    // altri metodi specifici per la struttura (calcolo lunghezza max e baricentro??)
};
struct Traces{
    unsigned int numTraces = 0; // numero tracce totali
    vector<unsigned int> IdTraces = {};  // vettore con Id tracce
    vector<Matrix<double, 3, 2>> CoordinatesEstremiTraces = {}; // vettore con matrice((x,y,z), estremo1 x estremo2)
    vector<double> lengthTraces = {};  // vettore lunghezza tracce

    //Per ogni traccia memorizzo i due id delle fratture
    vector<array<unsigned int,2>> IdsFractures;

    // Per ogni frattura, elenco delle tracce passanti assegnate
    vector<vector<unsigned int>> TraceIdsPassxFracture;
    // Per ogni frattura, elenco delle tracce non passanti assegnate
    vector<vector<unsigned int>> TraceIdsNoPassxFracture;
    // vector<vector<vector<unsigned int>>> TriangulatePolygons();
    // vector<double> computePolygonsArea();


    Traces() = default;
    Traces(const vector<vector<unsigned int>>& TraceIdsPassxFracture,
          const vector<vector<unsigned int>>& TraceIdsNoPassxFracture):
         TraceIdsPassxFracture( TraceIdsPassxFracture),
         TraceIdsNoPassxFracture(TraceIdsNoPassxFracture)
    {}

};

// 2^PARTE
struct Polygons{

    // from Cell0
    // number of cell0D
    unsigned int NumberCell0D = 0;
    // Id of each point in mesh
    vector<unsigned int> Cell0DId = {};
    // coordinates of each point
    vector<Vector3d> Cell0DCoordinates = {};
    /*
    // map for linking each marker (key) with the list of vertices Id associated with that marker
    map<unsigned int, list<unsigned int>> Cell0DMarkers = {};
    */

    // from Cell1D
    // number of cell1D
    unsigned int NumberCell1D = 0;
    // Id of each edge in mesh
    vector<unsigned int> Cell1DId = {};
    // two vertices Id of each edge (origin and end)
    vector<Vector2i> Cell1DVertices = {};
    /*
    // map for linking each marker (key) with the list of edges Id associated with that marker
    map<unsigned int, list<unsigned int>> Cell1DMarkers = {};
    */

    // from Cell2D
    // number of Cell2D
    unsigned int NumberCell2D= 0;
    // Id of each polygon in mesh
    vector<unsigned int> Cell2DId = {};
    // number of polygon's vertices
    list<unsigned int> NumberVertices = {};
    // vertices Id of each polygon
    vector<vector<unsigned int>> Cell2DVertices = {};
    // number of polygon's edges
    list<unsigned int> NumberEdges  = {};
    // edges Id of each polygon
    vector<vector<unsigned int>> Cell2DEdges = {};
    /*
    // map for linking each marker (key) with the list of polygon Id associated with that marker
    map<unsigned int, list<unsigned int>> Cell2DMarkers = {};
    */

    // strutture di supporto
    vector<MatrixXd> SequenzeXpunto;
    /// numCol = numSequenze (=num sottopoligoni)
    /// numRow = numTraces
    vector<Vector3d> CoordinatesPunto;
    /// mi servono per valutare se ogni punto è a destra o sinistra
};
}

///parentesi fine geometryNamespace
//double tol=numeric_limits<double>::epsilon() * 100;


#endif
