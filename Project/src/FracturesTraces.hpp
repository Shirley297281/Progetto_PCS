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

};

// 2^PARTE
struct Polygons{
    //...
};
///parentesi fine geometryNamespace
//double tolDef=numeric_limits<double>::epsilon() * 100;
}


#endif
