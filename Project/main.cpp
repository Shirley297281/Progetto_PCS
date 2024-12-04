#include "src/FracturesTracesPolygons.hpp"
#include "src/Utils_partTwo.hpp"
#include <iostream>
#include "src/namespace.hpp"
#include "TestingParaview/Code/src/UCDUtilities.hpp" //per Paraview esportazione


using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;

int main()
{

    ///
    /// Prima Parte
    /// riferimento a utils
    ///
    ///
    cout<<"\n\n---PRIMA PARTE---\n\n";

    Fractures fracture;
    string filename = "FR3_data.txt";

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


    ///
    /// Seconda Parte
    /// riferimento a utils_partTwo
    ///
    ///
    cout<<"\n\n---SECONDA PARTE---\n\n";

    // scelta dell'id della frattura di si vuole creare la mesh
    Polygons sottoPoligono;
    unsigned int z = 0;

    cout << " ### ANALISI Frattura con ID "<< z <<endl;

    // INIZIO SALVATAGGIO PUNTI DA TRACCE PASSANTI
    cout << "\n -SALVATAGGIO PUNTI DA TRACCE PASSANTI-\n\n";
    MemorizzaVerticiPassanti_Cell0Ds(fracture, trace, sottoPoligono, z);
    // FINE SALVATAGGIO PUNTI DA TRACCE PASSANTI

    // INIZIO SALVATAGGIO PUNTI DA TRACCE NON PASSANTI
    cout << "\n -SALVATAGGIO PUNTI DA TRACCE NON PASSANTI-\n\n";

    vector<Matrix<double, 3, 4>> NuoviEstremi = {};   // mi salvo i nuovi estremi delle tracce non passanti
                                                      // e i vettori direttori delle rette passanti per i nuovi estremi delle tracce non passanti
    NuoviEstremi.reserve(trace.TraceIdsNoPassxFracture[z].size());

    MemorizzaVerticiNonPassanti_Cell0Ds (fracture, trace, sottoPoligono, z, NuoviEstremi);
        // FINE SALVATAGGIO PUNTI DA TRACCE NON PASSANTI

    cout << " " <<endl;
    cout << " " <<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Ciclo per stampare gli elementi
    cout << " CONTROLLO SALVATAGGIO PUNTI DELLA MESH" <<endl;
    for (size_t i = 0; i < sottoPoligono.Cell0DId.size(); ++i) {
        cout << "Id punto:    " <<  sottoPoligono.Cell0DId[i] << endl;
        cout << "Coordinate:    " ;
        cout << sottoPoligono.Cell0DCoordinates[i].transpose() << endl;
    }
    cout << " " <<endl;
    cout << " " <<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // INIZIO CREAZIONI SEQUENZE CON TRACCE PASSANTI
    Creazioni_Sequenze_Passanti(fracture, trace, sottoPoligono, z);
    // FINE CREAZIONI SEQUENZE CON TRACCE PASSANTI

    // INIZIO CREAZIONI SEQUENZA CON TRACCE NON PASSANTI
    Creazioni_Sequenze_NONPassanti(fracture, sottoPoligono, z, NuoviEstremi);
    // FINE CREAZIONI SEQUENZE CON TRACCE NON PASSANTI

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "CONTROLLO CREAZIONI SEQUENZA"<<endl;
    cout << "numero punti con sequenze = "<< sottoPoligono.SequenzeXpunto.size() << endl;
    for (size_t k = 0; k < sottoPoligono.SequenzeXpunto.size(); ++k) {
        const MatrixXd& matrice = sottoPoligono.SequenzeXpunto[k];
        int numRighe = matrice.rows();
        int numColonne = matrice.cols();

        cout << "Id: " << k << ":" << endl;
        cout << "Numero di righe: " << numRighe << endl;
        cout << "Numero di colonne: " << numColonne << endl;

        for (int row = 0; row < numRighe; ++row) {
            for (int col = 0; col < numColonne; ++col) {
                cout << matrice(row, col) << " ";
            }
            cout << endl; // Fine della riga (colonna della matrice)
        }
        cout << "Fine id " << k << endl << endl; // Fine della matrice
    }
    cout << " " <<endl;
    cout << " " <<endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////


    // INIZIO CONFRONTI SEQUENZE e RAGGRUPPAMENTI IN POLIGONI
    vector<list<unsigned int>> VettSequenza_Punto = {};
    vector<VectorXd> LinktraSequenzaELista = {};
    // la posizione nei due vettori corrisponde a un sottopoligono nella frattura
    for(unsigned int i = 0; i < sottoPoligono.NumberCell0D; i++) //ciclo per i punti. i = identificativo punto in Cell0D
    {
        MatrixXd A = sottoPoligono.SequenzeXpunto[i]; //estraggo la matrice
        for (int j = 0; j < A.cols(); j++) // ciclo sulle colonne
        {
            VectorXd SequenzaJ = A.col(j); // estraggo la sequenza j-esima del punto i
            /// Inserisco la sequenza nel vettore VettSequenza_Punto come chiave se non c'è già, altrimenti inserisco nella lista
            // controllo se non ha ancora elementi
            bool trovato = false;
            if (VettSequenza_Punto.size() == 0)
            {
                LinktraSequenzaELista = {SequenzaJ};
                VettSequenza_Punto = {{i}} ;
                continue;
            }
            // controllo se la sequenza è già inserita
            for (unsigned int z = 0; z < LinktraSequenzaELista.size(); z++)
            {
                VectorXd SequenzaAttuale = LinktraSequenzaELista[z];
                if (Vettore::operator==(SequenzaJ, SequenzaAttuale)) // ho trovato la sequenza
                {
                    VettSequenza_Punto[z].push_back(i); // aggiungo allora nella lista l'id del punto
                    trovato = true;
                    break;

                }
            }
            if (trovato == true)
            {
                continue;
            }
            else
            {
                // non ho trovato la sequenza in LinktraSequenzaELista quindi devo inserire SequenzaJ sia in LinktraSequenzaELista e in VettSequenza_Punto
                // inserisco in coda, se inserisco entrambe in coda viene mantenuto l'ordine: in LinktraSequenzaELista in posizione 4 (per esempio)
                // ho la sequenza (1,0,1,1,0) allora in posizione 4 di VettSequenza_Punto ho la lista di id dei punti associati a (1,0,1,1,0)
                LinktraSequenzaELista.push_back(SequenzaJ);
                VettSequenza_Punto.push_back({i});
            }

        }
    }
    // FINE CONFRONTI SEQUENZE e RAGGRUPPAMENTI IN POLIGONI

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "CONTROLLO RAGGRUPPAMENTO PUNTI IN SOTTOPOLIGONI"<<endl;
    for (size_t i = 0; i < LinktraSequenzaELista.size(); ++i) {
        cout << "sottopoligono: " << i << " : " << LinktraSequenzaELista[i].transpose() << endl;
        cout << "ids: ";
        for (auto it = VettSequenza_Punto[i].begin(); it != VettSequenza_Punto[i].end(); ++it) {
            cout << *it << " ";
        }
        cout << endl;
    }
    cout << " " <<endl;
    cout << " " <<endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // INIZIO RIMOZIONE DI SEQUENZE CON N°ID PUNTI <= 2
    vector<list<unsigned int>> VettSequenza_PuntoDesiderato1 = {};
    vector<VectorXd> LinktraSequenzaEListaDesiderato1 = {};
    VettSequenza_PuntoDesiderato1.reserve(VettSequenza_Punto.size());
    for (unsigned int i = 0; i < VettSequenza_Punto.size(); ++i)
    {
        if (VettSequenza_Punto[i].size() > 2)
        {
            VettSequenza_PuntoDesiderato1.push_back(VettSequenza_Punto[i]);
            LinktraSequenzaEListaDesiderato1.push_back(LinktraSequenzaELista[i]);
        }
    }
    // FINE RIMOZIONE DI SEQUENZE CON N°ID PUNTI <= 2

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "CONTROLLO RIMOZIONI SEQUENZE"<<endl;
    for (size_t i = 0; i < LinktraSequenzaEListaDesiderato1.size(); ++i) {
        cout << "sottopoligono: " << i << " : " << LinktraSequenzaEListaDesiderato1[i].transpose() << endl;
        cout << "ids: ";
        for (auto it = VettSequenza_PuntoDesiderato1[i].begin(); it != VettSequenza_PuntoDesiderato1[i].end(); ++it) {
            cout << *it << " ";
        }
        cout << endl;
    }
    cout << " " <<endl;
    cout << " " <<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // INIZIO ORDINAMENTO LATI E SALVATAGGIO IN CELL2D
    size_t numSottopoligoni = VettSequenza_PuntoDesiderato1.size(); // ogni sottopoligono è univocamente determinato da una sequenza: numSottoPol = numSequenze
    sottoPoligono.NumberCell2D = numSottopoligoni;

    sottoPoligono.Cell2DEdges.resize(numSottopoligoni);
    sottoPoligono.Cell2DVertices.resize(numSottopoligoni);

    for (unsigned int i = 0; i < numSottopoligoni; i++)
    {
        list<unsigned int> listaIdVertici = VettSequenza_PuntoDesiderato1[i];
        Creo_sottopoligono(z, i,listaIdVertici, sottoPoligono, fracture);
    }
    // FINE ORDINAMENTO LATI E SALVATAGGIO IN CELL2D per FILE DA 3 FRATTURE

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout << "CONTROLLO ORDINAMENTO per FILE da 3 FRATTURE"<<endl;
    for(unsigned int j = 0; j < sottoPoligono.Cell1DId.size(); j++){
        cout << "gli estremi del lato con id " << sottoPoligono.Cell1DId[j] << " sono: " << sottoPoligono.Cell1DVertices[j][0] << " e " << sottoPoligono.Cell1DVertices[j][1] << endl;
    }

    cout << "NumberCell1D: " << sottoPoligono.NumberCell1D << endl;
    cout << endl;
    cout << " " <<endl;
    cout << " " <<endl;
    cout << "CONTROLLO RIEMPIMENTO di Cell2DVertices"<<endl;
    //verifico che sia giusto il riempimento di Cell2DVertices -> GIUSTO
    cout << "id degli estremi del sottopoligono 0" << endl;
    for(unsigned int i = 0; i < 4; i++){
        cout << "estremo " << i << " : " << sottoPoligono.Cell2DVertices[0][i] << endl;
    }
    cout << "id degli estremi del sottopoligono 1" << endl;
    for(unsigned int i = 0; i < 4; i++){
        cout << "estremo " << i << " : " << sottoPoligono.Cell2DVertices[1][i] << endl;
    }
    cout << "id degli estremi del sottopoligono 2" << endl;
    for(unsigned int i = 0; i < 4; i++){
        cout << "estremo " << i << " : " << sottoPoligono.Cell2DVertices[2][i] << endl;
    }
    cout << " " <<endl;
    cout << " " <<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    ///
    ///EXPORTING PARAVIEW
    ///
    ///
    ///

    cout<<"\n\n--Esportazione per PARAVIEW--"<<endl;

    // RIMOZIONE 2: SOTTOPOLIGONI FORMATI DA ALTRI SOTTOPOLIGONI
    vector<list<unsigned int>> VettSequenza_PuntoDesiderato2 = {{0, 5, 6, 7}, {1, 2, 4, 5},  {3, 4 ,6, 7}} ;



    Gedim::UCDUtilities exporter;
    std::vector<std::vector<unsigned int>> triangles;
    Eigen::VectorXi materials;
    sottoPoligono.GedimInterface(triangles, materials);
    cout<<"triangolazione sottopoligoni effettuata."<<endl;


    // Controlla se il vettore di imput è vuoto
    if (sottoPoligono.Cell2DVertices.empty()) {
       throw runtime_error("Cell2DVertices is empty");
    }

    // Crea MatrixXd per aver una struttura compatibile con la funzione in UDCUtilities
    MatrixXd eigenMatrix(3, sottoPoligono.NumberCell0D);

    for (unsigned int p = 0; p < sottoPoligono.NumberCell0D; ++p) {
        Vector3d coord = sottoPoligono.Cell0DCoordinates[p];
        eigenMatrix(0, p) = coord[0];
        eigenMatrix(1, p) = coord[1];
        eigenMatrix(2, p) = coord[2];
    }


    exporter.ExportPolygons("./Polygon0_FR3.inp",
                            eigenMatrix,
                            triangles,
                            {},
                            {},
                            materials);

    cout << "stampato file 'Polygon0_FR3.inp'."<<endl;

    // Numero di linee
    const int numLines = sottoPoligono.NumberCell1D;

    // Creazione della matrice MatrixXi per le linee
    MatrixXi lines(2, numLines);

    for (int i = 0; i < numLines; ++i) {
        lines(0, i) = sottoPoligono.Cell1DVertices[i][0];
        lines(1, i) = sottoPoligono.Cell1DVertices[i][1];
    }

    exporter.ExportSegments("./Segments0_FR3.inp",
                            eigenMatrix,
                            lines,
                            {},
                            {},
                            materials);

    cout << "stampato file 'Segments0_FR3.inp'."<<endl;

    ///fine PARAVIEW
    ///
    /// */
    return 0;
}
