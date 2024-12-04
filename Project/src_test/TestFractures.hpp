#ifndef _TestFractures_test_HPP_
#define _TestFractures_test_HPP_

#include <fstream>
#include <iostream>
#include "FracturesTracesPolygons.hpp"
#include "Utils.hpp"
#include <gtest/gtest.h>
#include "Eigen/Eigen"
#include <math.h>
#include "namespace.hpp"
#include "inline.hpp"
#include "Utils_partTwo.hpp"


using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;

// test su fratture e tracce
TEST(TestFractures, TestDFNReading)
{
    ifstream file("FR2_datatest.txt");
    if (!file.is_open()) {
        cerr << "Error opening file FR2_datatest.txt" << endl;
        exit(1);
    }

    Fractures fractures;

    // Read the fractures from the file
    bool importSuccess = ImportFR("FR2_datatest.txt", fractures);
    EXPECT_TRUE(importSuccess);
    if(!importSuccess){
        cerr << "Error importing fractures from file." << endl;
        return;
    }
    // Check that the number of fractures is correct
    EXPECT_EQ(fractures.NumFractures, 2);

    // Check the properties of the first fracture
    if(fractures.CoordinatesVertice.size() > 0){
        MatrixXd coord1 = fractures.CoordinatesVertice[0];
        EXPECT_EQ(fractures.IdFractures[0], 0);
        EXPECT_EQ(fractures.numVertices[0], 3);
        EXPECT_DOUBLE_EQ(coord1(0,0), 0.0);
        EXPECT_DOUBLE_EQ(coord1(1,0), 0.0);
        EXPECT_DOUBLE_EQ(coord1(2,0), 0.0);
        EXPECT_DOUBLE_EQ(coord1(0,1), 1.0);
        EXPECT_DOUBLE_EQ(coord1(1,1), 0.0);
        EXPECT_DOUBLE_EQ(coord1(2,1), 0.0);
        EXPECT_DOUBLE_EQ(coord1(0,2), 1.0);
        EXPECT_DOUBLE_EQ(coord1(1,2), 1.0);
        EXPECT_DOUBLE_EQ(coord1(2,2), 0.0);
    }else{
        cerr << "No coordinates found for the first fracture" << endl;
    }


    // Check the properties of the second fracture
    if(fractures.CoordinatesVertice.size() > 0){
        MatrixXd coord2 = fractures.CoordinatesVertice[1];
        EXPECT_EQ(fractures.IdFractures[1], 1);
        EXPECT_EQ(fractures.numVertices[1], 4);
        EXPECT_DOUBLE_EQ(coord2(0,0), 0.8);
        EXPECT_DOUBLE_EQ(coord2(1,0), 0.0);
        EXPECT_DOUBLE_EQ(coord2(2,0), -0.1);
        EXPECT_DOUBLE_EQ(coord2(0,1), 0.8);
        EXPECT_DOUBLE_EQ(coord2(1,1), 0.0);
        EXPECT_DOUBLE_EQ(coord2(2,1), 0.29999999999999999);
        EXPECT_DOUBLE_EQ(coord2(0,2), 0.8);
        EXPECT_DOUBLE_EQ(coord2(1,2), 1.0);
        EXPECT_DOUBLE_EQ(coord2(2,2), 0.29999999999999999);
        EXPECT_DOUBLE_EQ(coord2(0,3), 0.8);
        EXPECT_DOUBLE_EQ(coord2(1,3), 1.0);
        EXPECT_DOUBLE_EQ(coord2(2,3), -0.1);
    }else{
        cerr << "No coordinates found for the second fracture" << endl;
    }
}


// test sulle function contenute in Finding_Traces.cpp

// test sulla function distinzioneTipoTraccia1 (5 casistiche)
TEST(TestdistinzioneTipoTraccia1, Test_1pass_1notpass) // 1^ caso
{
    // candidati estremi delle tracce
    Vector3d v1(1.0, 0.0, 1.0);
    Vector3d v2(-1.0, 0.0, 1.0);
    Vector3d v3(3.0, 0.0, 1.0);
    Vector3d v4(-3.0, 0.0, 1.0);

    Fractures fracture;
    Traces trace;
    vector<Vector3d> vecI;
    vecI.reserve(2);
    vecI.push_back(v1);
    vecI.push_back(v2);

    vector<Vector3d> vecJ;
    vecJ.reserve(2);
    vecJ.push_back(v3);
    vecJ.push_back(v4);

    Vector3d Point = v1;
    Vector3d t(2, 0, 0);
    unsigned int i = 0;
    unsigned int j = 1;
    int sol;
    sol = distinzioneTipoTraccia1( trace, vecI, vecJ, Point, t, i, j);
    EXPECT_EQ(sol, 4);
}

TEST(TestdistinzioneTipoTraccia1, Test_1pass_1notpass2) // 2^ caso
{
    // candidati estremi delle tracce
    Vector3d v1(1.0, 0.0, 1.0);
    Vector3d v2(-1.0, 0.0, 1.0);
    Vector3d v3(3.0, 0.0, 1.0);

    Fractures fracture;
    Traces trace;
    vector<Vector3d> vecI;
    vecI.reserve(2);
    vecI.push_back(v1);
    vecI.push_back(v2);

    vector<Vector3d> vecJ;
    vecJ.reserve(2);
    vecJ.push_back(v2);
    vecJ.push_back(v3);

    Vector3d Point = v1;
    Vector3d t(2, 0, 0);
    unsigned int i = 0;
    unsigned int j = 1;
    int sol;
    sol = distinzioneTipoTraccia1(trace, vecI, vecJ, Point, t, i, j);
    EXPECT_EQ(sol, 4);
}

TEST(TestdistinzioneTipoTraccia1, Test_Notrace) // 3^ caso
{
    // candidati estremi delle tracce
    Vector3d v1(5.0, 0.0, 1.0);
    Vector3d v2(3.0, 0.0, 1.0);
    Vector3d v3(1.0, 0.0, 1.0);
    Vector3d v4(-1.0, 0.0, 1.0);

    Fractures fracture;
    Traces trace;
    vector<Vector3d> vecI;
    vecI.reserve(2);
    vecI.push_back(v1);
    vecI.push_back(v2);

    vector<Vector3d> vecJ;
    vecJ.reserve(2);
    vecJ.push_back(v3);
    vecJ.push_back(v4);

    Vector3d Point = v1;
    Vector3d t(2, 0, 0);
    unsigned int i = 0;
    unsigned int j = 1;
    int sol;
    sol = distinzioneTipoTraccia1(trace, vecI, vecJ, Point, t, i, j);
    EXPECT_EQ(sol, 1);
}

TEST(TestdistinzioneTipoTraccia1, Test_2coindenteExtremes) // 4^ caso
{
    // candidati estremi delle tracce
    Vector3d v1(1.0, 0.0, 1.0);
    Vector3d v2(-1.0, 0.0, 1.0);
    Fractures fracture;
    Traces trace;
    vector<Vector3d> vecI;
    vecI.reserve(2);
    vecI.push_back(v1);
    vecI.push_back(v2);

    vector<Vector3d> vecJ;
    vecJ.reserve(2);
    vecJ.push_back(v1);
    vecJ.push_back(v2);

    Vector3d Point = v2;
    Vector3d t(-2, 0, 0);
    unsigned int i = 0;
    unsigned int j = 1;
    int sol;
    sol = distinzioneTipoTraccia1( trace, vecI, vecJ, Point, t, i, j);
    EXPECT_EQ(sol, 3);
}

TEST(TestdistinzioneTipoTraccia1, Test_2notpass) // 5^ caso
{
    // candidati estremi delle tracce
    Vector3d v1(1.0, 0.0, 1.0);
    Vector3d v2(-1.0, 0.0, 1.0);
    Vector3d v3(3.0, 0.0, 1.0);
    Vector3d v4(0.0, 0.0, 1.0);

    Fractures fracture;
    Traces trace;
    vector<Vector3d> vecI;
    vecI.reserve(2);
    vecI.push_back(v1);
    vecI.push_back(v2);

    vector<Vector3d> vecJ;
    vecJ.reserve(2);
    vecJ.push_back(v3);
    vecJ.push_back(v4);

    Vector3d Point = v1;
    Vector3d t(2, 0, 0);
    unsigned int i = 0;
    unsigned int j = 1;
    int sol;
    sol = distinzioneTipoTraccia1(trace, vecI, vecJ, Point, t, i, j);
    EXPECT_EQ(sol, 6);
}

// test su Calcolo_par
// t (direzione della retta di intersezione fra piani) coincide con uno dei lati della frattura
TEST(TestCalcolo_par, Test_t_coincident){
    Vector3d a(-2.0, 0.0, 0.0);
    Vector3d b(0.0, 0.0, 2.0);
    Vector3d c(0.0, 0.0, 0.0);
    Vector3d d(2.0, 0.0, 2.0);
    MatrixXd A(3,4);
    A.col(0) << d;
    A.col(1) << b;
    A.col(2) << a;
    A.col(3) << c;

    Vector3d t(2.0, 0.0, 2.0);
    Vector3d Point = a;
    int i = 0;
    vector<Vector3d> vec = {};

    Fractures fracture;
    fracture.CoordinatesVertice.push_back(A);
    unsigned int par;
    par = Calcolo_par(t, Point, i, vec, fracture);
    EXPECT_EQ(par, 2);
}

// t è esterna alla frattura
TEST(TestCalcolo_par, Test_t_external){
    Vector3d a(-2.0, 0.0, 0.0);
    Vector3d b(0.0, 0.0, 2.0);
    Vector3d c(0.0, 0.0, 0.0);
    Vector3d d(2.0, 0.0, 2.0);
    MatrixXd A(3,4);
    A.col(0) << d;
    A.col(1) << b;
    A.col(2) << a;
    A.col(3) << c;

    Vector3d t(2.0, 0.0, 2.0); // sto trovando una retta parallela al lato ab
    Vector3d Point(-3.0, 0.0, 0.0);
    int i = 0;
    vector<Vector3d> vec = {};

    Fractures fracture;
    fracture.CoordinatesVertice.push_back(A);
    unsigned int par;
    par = Calcolo_par(t, Point, i, vec, fracture);
    EXPECT_EQ(par, 0);
}

// t interseca due lati della frattura
TEST(TestCalcolo_par, Test_t_intersect_2edges){
    Vector3d a(-2.0, 0.0, 0.0);
    Vector3d b(0.0, 0.0, 2.0);
    Vector3d c(0.0, 0.0, 0.0);
    Vector3d d(2.0, 0.0, 2.0);
    MatrixXd A(3,4);
    A.col(0) << d;
    A.col(1) << b;
    A.col(2) << a;
    A.col(3) << c;

    Vector3d t(2.0, 0.0, 2.0);
    Vector3d Point(1.0, 0.0, 2.0);
    int i = 0;
    vector<Vector3d> vec = {};

    Fractures fracture;
    fracture.CoordinatesVertice.push_back(A);
    unsigned int par;
    par = Calcolo_par(t, Point, i, vec, fracture);
    EXPECT_EQ(par, 2);
}

// t interseca un vertice della frattura
TEST(TestCalcolo_par, Test_t_intersect_1vertice){
    Vector3d a(-2.0, 0.0, 0.0);
    Vector3d b(0.0, 0.0, 2.0);
    Vector3d c(0.0, 0.0, 0.0);
    Vector3d d(2.0, 0.0, 2.0);
    MatrixXd A(3,4);
    A.col(0) << d; // la riempio così perchè i vertici devono essere ordinati e consecutivi
    A.col(1) << b;
    A.col(2) << a;
    A.col(3) << c;

    Vector3d t(-2.0, 0.0, 1.0);
    Vector3d Point = a;
    int i = 0;
    vector<Vector3d> vec = {};

    Fractures fracture;
    fracture.CoordinatesVertice.push_back(A);
    unsigned int par;
    par = Calcolo_par(t, Point, i, vec, fracture);
    EXPECT_EQ(par, 1);
}

TEST(TestCalcolo_par, Test_t_intersect_2vertices){
    Vector3d a(-2.0, 0.0, 0.0);
    Vector3d b(0.0, 0.0, 2.0);
    Vector3d c(0.0, 0.0, 0.0);
    Vector3d d(2.0, 0.0, 2.0);
    MatrixXd A(3,4);
    A.col(0) << d;
    A.col(1) << b;
    A.col(2) << a;
    A.col(3) << c;

    Vector3d t(0.0, 0.0, -2.0);
    Vector3d Point = c;
    int i = 0;
    vector<Vector3d> vec = {};

    Fractures fracture;
    fracture.CoordinatesVertice.push_back(A);
    unsigned int par;
    par = Calcolo_par(t, Point, i, vec, fracture);
    EXPECT_EQ(par, 2);
}


// test su function implementate nel file Utils.cpp
TEST(TestUtils, TestBarycenter)
{
    Vector3d p1(0.0, 0.0, 2.0);
    Vector3d p2(0.0, 2.0, 0.0);
    Vector3d p3(0.0, 0.0, 0.0);
    MatrixXd A(3,3);
    A.col(0) << p1;
    A.col(1) << p2;
    A.col(2) << p3;

    array<double,3> barycenter_coord = barycenter(A, 3);
    EXPECT_DOUBLE_EQ(barycenter_coord[0], 0.0);
    EXPECT_DOUBLE_EQ(barycenter_coord[1], 2.0/3);
    EXPECT_DOUBLE_EQ(barycenter_coord[2], 2.0/3);
}

TEST(TestUtils, TestCheckSphere)
{
    array <double, 3> bar1{1.0, 1.0, 1.0};
    array <double, 3> bar2{3.0, 1.0, 2.0};
    double l1 = 2.0;
    double l2 = 4.0;
    bool check = check_sphere(bar1, bar2, l1, l2);
    ASSERT_TRUE(check); // se prendo due punti che distano di meno di l1+l2 => restituisce TRUE (la distanza vera è sqrt(5))
}

TEST(TestUtils, TestNormalVector)
{
    MatrixXd m = (Matrix3d(3,3) << -1.0, 0.0, -1.0,
                  1.0, 0.0, -2.0,
                  -2.0, 0.0, -3.0).finished() ;
    Vector3d v = normal_vector(m);
    EXPECT_DOUBLE_EQ(v[0], -7.0);
    EXPECT_DOUBLE_EQ(v[1], -1.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
}


// test sulle inline function
TEST(TestInline, TestEuclideanDistance)
{

    Vector3d p1(0.0, 0.0, 2.0);
    Vector3d p2(0.0, 2.0, 0.0);
    Vector3d p3(0.0, 0.0, 0.0);
    MatrixXd A(3,3);
    A.col(0) << p1;
    A.col(1) << p2;
    A.col(2) << p3;

    double d12 = euclidean_distance(p1, p2);
    double d13 = euclidean_distance(p1, p3);

    // verifico che la distanza non sia negativa -> sto dicendo che se la distanza è <0.0 => il programma deve restituire un errore
    ASSERT_GE(d12, 0.0);
    ASSERT_GE(d13, 0.0);

    // verifico che calcoli correttamente la distanza
    EXPECT_DOUBLE_EQ(d12, sqrt(8.0));
    EXPECT_DOUBLE_EQ(d13, 2.0);

    double max = max_euclidean_distance(A, 3);

    // verifico che selezioni correttamente la distanza massima
    EXPECT_DOUBLE_EQ(max, sqrt(8.0));

}

TEST(TestInline, TestVecProdZero)
{
    Vector3d v = {};
    Vector3d v1(1.0, 2.0, -1.0);
    Vector3d v2(3.0, 6.0, -3.0);
    v = vec_product(v1, v2);
    EXPECT_DOUBLE_EQ(v[0], 0.0);
    EXPECT_DOUBLE_EQ(v[1], 0.0);
    EXPECT_DOUBLE_EQ(v[2], 0.0);
}

TEST(TestInline, TestVecProd)
{
    Vector3d v = {};
    Vector3d v1(1.0, 2.0, -1.0);
    Vector3d v2(-3.0, 4.0, 5.0);
    v = vec_product(v1, v2);
    EXPECT_DOUBLE_EQ(v[0], 14.0);
    EXPECT_DOUBLE_EQ(v[1], -2.0);
    EXPECT_DOUBLE_EQ(v[2], 10.0);
}

// controllo che l'inserimento attraverso la function inserimento_map() funzioni correttamente
TEST(TestInline, Test_inserimento_map){
    // caso 1 passante 1 non passante di Test_1pass_1notpass
    Vector3d v1(1.0, 0.0, 1.0);
    Vector3d v2(-1.0, 0.0, 1.0);
    Vector3d v3(3.0, 0.0, 1.0);
    Vector3d v4(-3.0, 0.0, 1.0);

    Traces trace;
    vector<Vector3d> vecI;
    vecI.reserve(2);
    vecI.push_back(v1);
    vecI.push_back(v2);

    vector<Vector3d> vecJ;
    vecJ.reserve(2);
    vecJ.push_back(v3);
    vecJ.push_back(v4);

    trace.numTraces = 1; // se ho una traccia => il suo identificativo è 0
    unsigned int i = 0; // idpar della frattura i
    unsigned int j = 1; // idpar della frattura j
    int pass0 = 0;
    int pass1 = 1;
    inserimento_map(pass0, i, trace); // restituisce passante per i
    inserimento_map(pass1, j, trace); // restituisce non passante per j
    int id_pass = trace.TraceIdsPassxFracture[i][0];
    int id_nopass = trace.TraceIdsNoPassxFracture[j][0];
    ASSERT_EQ(id_pass, 0);
    ASSERT_EQ(id_nopass, 0);
}

// p non appartiene al segmento v1v2
TEST(TestInline, Test_combinazione_convessa_si){
    Vector3d v1(0.0, 0.0, 0.0);
    Vector3d v2(2.0, 2.0, 2.0);
    Vector3d p(1.0, 1.0, 1.0);
    bool sol = combinazione_convessa(v1, v2, p);
    ASSERT_TRUE(sol);
}

// p appartiene al segmento v1v2
TEST(TestInline, Test_combinazione_convessa_no){
    Vector3d v1(0.0, 0.0, 0.0);
    Vector3d v2(2.0, 2.0, 2.0);
    Vector3d p(3.0, 3.0, 3.0);
    bool sol = combinazione_convessa(v1, v2, p);
    ASSERT_FALSE(sol);
}


// test su intersezione_piani con vettori v1, v2 (vettori normali ai piani) tra loro paralleli
TEST(TestIntersezione_piani, TestVecParallel)
{
    // i vettori v1 e v2 sono linearmente dipendenti -> Non hanno un punto di intersezione unico
    Vector3d n1(1.0, 2.0, 3.0);
    Vector3d n2(2.0, 4.0, 6.0);
    array <double, 3> b1 = {0.0, 1.0, 2.0};
    array <double, 3> b2= {0.0, 0.0, 1.0};
    Vector3d t(1.0, 2.0, 0.0);
    Vector3d Point = {};
    bool sol = intersezione_piani(n1, n2, b1, b2, t, Point); // il sistema dovrebbe restituire 'false' perchè n1 e n2 sono paralleli
    ASSERT_FALSE(sol); // verifica che la variabile sol sia falsa
}

// test su intersezione_piani con vettori v1, v2 paralleli e coincidenti
TEST(TestIntersezione_piani, TestVecCoincident)
{
    Vector3d n1(1.0, 2.0, 3.0);
    Vector3d n2(1.0, 2.0, 3.0);
    array <double, 3> b1 = {0.0, 1.0, 2.0};
    array <double, 3> b2 = {0.0, 0.0, 1.0};
    Vector3d t(1.0, 2.0, 0.0);
    Vector3d Point = {};
    bool sol = intersezione_piani(n1, n2, b1, b2, t, Point); // il sistema dovrebbe restituire 'false' perchè n1 e n2 sono coincidenti
    ASSERT_FALSE(sol); // verifica che la variabile sol sia falsa
}

// test su intersezione_piani con vettori v1, v2 nè paralleli nè coincidenti
TEST(TestIntersezione_piani, TestCorrectSol)
{
    Vector3d n1(1.0, 2.0, 3.0);
    Vector3d n2(3.0, 5.0, 7.0);
    array <double, 3> b1 = {0.0, 1.0, 2.0};
    array <double, 3> b2 = {0.0, 0.0, 1.0};
    Vector3d t(1.0, 2.0, 0.0);
    Vector3d Point = {};
    bool sol = intersezione_piani(n1, n2, b1, b2, t, Point); // il sistema dovrebbe restituire 'true' perchè n1 e n2 non sono nè paralleli nè coicidenti
    ASSERT_TRUE(sol); // verifica che la variabile sol sia vera
}

// test su intersezione_rette con vettori t e V1-V2 paralleli
TEST(TestIntersezione_rette, TestSoluzioneIntersezione_rette_vecParallel)
{
    Vector3d t(0.0, -2.0, 0.0);
    Vector3d V1(1.0, 2.0, 3.0);
    Vector3d V2(1.0, 6.0, 3.0); // così il vettoreDirezione = V1 - V2 = (0,-4,0)
    Vector3d Point(1.0, 6.0, 3.0);
    Vector3d Punto0 = {};
    bool sol = intersezione_rette(t, V1, V2, Point, Punto0);
    ASSERT_FALSE(sol);
}

// test su intersezione_rette con vettori t e V1-V2 NON paralleli
TEST(TestIntersezione_rette, TestSoluzioneIntersezione_rette_vecNotParallel)
{
    Vector3d t(0.0, -2.0, 0.0);
    Vector3d V1(2.0, 3.0, 4.0);
    Vector3d V2(1.0, 1.0, 1.0); // così il vettoreDirezione = V1 - V2 = (1,2,3)
    Vector3d Point(2.0, 1.0, 5.0);
    Vector3d Punto0 = {};
    bool sol = intersezione_rette(t, V1, V2, Point, Punto0);
    ASSERT_TRUE(sol);
}


/// TEST PARTE 2
// test su creo_sottopoligoni
TEST(TestSottopoligoni, Testcreo_sottopoligoni){
    unsigned int num_fracture = 0;
    vector<list<unsigned int>> VettSequenza_Punto = {};
    Fractures fracture;
    Polygons poligono;

    poligono.Cell0DCoordinates.reserve(6);
    poligono.Cell0DCoordinates[0] = {2.0, 1.0, 0.0};
    poligono.Cell0DCoordinates[1] = {-3.0, 1.0, 0.0};
    poligono.Cell0DCoordinates[2] = {-3.0, -2.0, 0.0};
    poligono.Cell0DCoordinates[3] = {2.0, -2.0, 0.0};
    poligono.Cell0DCoordinates[4] = {-3.0, 0.0, 0.0};
    poligono.Cell0DCoordinates[5] = {2.0, 0.0, 0.0};

    poligono.Cell0DId.reserve(6);
    poligono.Cell0DId = {0, 1, 2, 3, 4, 5, 6};

    fracture.vettoreNormalePiano.reserve(1);
    fracture.vettoreNormalePiano[num_fracture] = {0.0, 0.0, 1.0};

    list<unsigned int> listaIdVertici1;
    listaIdVertici1.push_back(0);
    listaIdVertici1.push_back(1);
    listaIdVertici1.push_back(4);
    listaIdVertici1.push_back(5);
    VettSequenza_Punto.push_back(listaIdVertici1);

    list<unsigned int> listaIdVertici2;
    listaIdVertici2.push_back(2);
    listaIdVertici2.push_back(3);
    listaIdVertici2.push_back(5);
    listaIdVertici2.push_back(4);
    VettSequenza_Punto.push_back(listaIdVertici2);


    poligono.Cell2DEdges.resize(2);
    poligono.Cell2DVertices.resize(2);
    poligono.NumberCell2D = 2;
    poligono.NumberCell0D = 6;

    for (unsigned int i = 0; i < 2; i++)
    {
        list<unsigned int> listaIdVertici = VettSequenza_Punto[i];
        Creo_sottopoligono(num_fracture, i,listaIdVertici, poligono, fracture);
    }

    EXPECT_EQ(poligono.Cell2DEdges[0][0], 0);
    EXPECT_EQ(poligono.Cell2DEdges[0][1], 1);
    EXPECT_EQ(poligono.Cell2DEdges[0][2], 2);
    EXPECT_EQ(poligono.Cell2DEdges[0][3], 3);
    EXPECT_EQ(poligono.Cell2DEdges[1][0], 4);
    EXPECT_EQ(poligono.Cell2DEdges[1][1], 5);
    EXPECT_EQ(poligono.Cell2DEdges[1][2], 2);
    EXPECT_EQ(poligono.Cell2DEdges[1][3], 6);

    EXPECT_EQ(poligono.Cell2DVertices[0][0], 0);
    EXPECT_EQ(poligono.Cell2DVertices[0][1], 1);
    EXPECT_EQ(poligono.Cell2DVertices[0][2], 4);
    EXPECT_EQ(poligono.Cell2DVertices[0][3], 5);
    EXPECT_EQ(poligono.Cell2DVertices[1][0], 2);
    EXPECT_EQ(poligono.Cell2DVertices[1][1], 3);
    EXPECT_EQ(poligono.Cell2DVertices[1][2], 5);
    EXPECT_EQ(poligono.Cell2DVertices[1][3], 4);

    EXPECT_EQ(poligono.NumberEdges.front(), 4);
    EXPECT_EQ(poligono.NumberVertices.front(), 4);
    EXPECT_EQ(poligono.NumberEdges.back(), 4);
    EXPECT_EQ(poligono.NumberVertices.back(), 4);
    EXPECT_EQ(poligono.NumberCell1D, 7);

    EXPECT_EQ(poligono.Cell2DId[0], 0);
    EXPECT_EQ(poligono.Cell2DId[1], 1);

    Vector2i v0 = poligono.Cell1DVertices[0];
    Vector2i v1 = poligono.Cell1DVertices[1];
    Vector2i v2 = poligono.Cell1DVertices[2];
    Vector2i v3 = poligono.Cell1DVertices[3];
    EXPECT_EQ(v0[0], 0);
    EXPECT_EQ(v0[1], 1);
    EXPECT_EQ(v1[0], 1);
    EXPECT_EQ(v1[1], 4);
    EXPECT_EQ(v2[0], 4);
    EXPECT_EQ(v2[1], 5);
    EXPECT_EQ(v3[0], 5);
    EXPECT_EQ(v3[1], 0);

    Vector2i v4 = poligono.Cell1DVertices[4];
    Vector2i v5 = poligono.Cell1DVertices[5];
    Vector2i v6 = poligono.Cell1DVertices[6];
    EXPECT_EQ(v4[0], 2);
    EXPECT_EQ(v4[1], 3);
    EXPECT_EQ(v5[0], 3);
    EXPECT_EQ(v5[1], 5);
    EXPECT_EQ(v6[0], 4);
    EXPECT_EQ(v6[1], 2);

    EXPECT_EQ(poligono.Cell1DId[0], 0);
    EXPECT_EQ(poligono.Cell1DId[1], 1);
    EXPECT_EQ(poligono.Cell1DId[2], 2);
    EXPECT_EQ(poligono.Cell1DId[3], 3);
    EXPECT_EQ(poligono.Cell1DId[4], 4);
    EXPECT_EQ(poligono.Cell1DId[5], 5);
    EXPECT_EQ(poligono.Cell1DId[6], 6);

}


#endif
