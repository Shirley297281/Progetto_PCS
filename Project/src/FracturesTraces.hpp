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
    vector<unsigned int> IdFractures ; //vettore con Id fratture
    vector<MatrixXd> CoordinatesVertice ; // matrice con coordinate dei vertici (ogni colonna individua un vertice)
    vector<double> lenghtMaxEdges = {} ; // vettore di lunghezza massima per ogni frattura
    vector<vector<double>> baricentro = {}; // vettore baricentro per ogni frattura
    vector<array<double, 3>> versoreNormalePiano = {}; // vettore di versori per ogni frattura (i,j,k)
    vector<unsigned int> numVertices = {}; //numero vertici per ogni poligono
    map<unsigned int, list<unsigned int>> TraceIdsxFracture = {}; //per ogni frattura elenco tracce assegnate
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
    vector<bool> vectorTips = {} ;  // vettore di bool (se passante o no)
    // vector<vector<vector<unsigned int>>> TriangulatePolygons();
    // vector<double> computePolygonsArea();

};

// 2^PARTE
struct Polygons{
    //...
};
///parentesi fine geometryNamespace
}

#endif
