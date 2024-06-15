#include "src/Utils.hpp"
#include "src/FracturesTraces.hpp"
#include "src/utils_partTwo.hpp"
#include <sstream>
#include <iostream>
#include "src/namespace.hpp"
#include "src/inline.hpp"


using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;

int main()
{
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
    ///Seconda Parte
    /// riferimento a utils_partTwo
    ///
    ///
    cout<<"\n\n---SECONDA PARTE---\n\n";

    // scelta da utente l'id della frattura di cui voglio vedere la sottopoligonazione
    Polygons sottoPoligono;
    unsigned int z = 1;

    // INIZIO SALVATAGGIO PUNTI
    MemorizzaVertici_Cell0Ds(fracture, trace, sottoPoligono, z);
    // FINE SALVATAGGIO PUNTI


    // INIZIO CREAZIONI SEQUENZE

    Creazioni_Sequenze(fracture, trace, sottoPoligono, z);

    // FINE CREAZIONI SEQUENZE


    // INIZIO CONFRONTI SEQUENZE e RAGGRUPPAMENTI IN POLIGONI
    vector<list<unsigned int>> VettSequenza_Punto = {};
    vector<VectorXd> LinktraSequenzaELista = {};
    // la posizione nei due vettori corrisponde a un poligono sulla frattura
    for(unsigned int i = 0; i < sottoPoligono.NumberCell0D; i++) //ciclo per i punti. i = identificativo punto in Cell0D
    {
        MatrixXd A = sottoPoligono.SequenzeXpunto[i]; //estraggo la matrice
        for (int j = 0; j < A.cols(); j++) // ciclo sulle colonne
        {
            VectorXd SequenzaJ = A.col(j); // estraggo la sequenza j-esima del punto i
            // la inserisco nel vettore VettSequenza_Punto come chiave se non c'è già, altrimenti inserisco nella lista
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
                continue; // posso uscire ho già fatto quello che dovevo fare
            }
            else
            {
                // non ho trovato la sequenza in LinktraSequenzaELista quindi devo inserire SequenzaJ sia in LinktraSequenzaELista e in VettSequenza_Punto
                // inserisco in coda, in teoria se inserisco entrambe in coda dovrei mantenere l'ordine: in posizione 4 (per esempio) ho la (1,0,1,1,0)
                // in LinktraSequenzaELista allora in posizione 4 di VettSequenza_Punto ho la lista di id dei punti associati a (1,0,1,1,0)
                LinktraSequenzaELista.push_back(SequenzaJ);
                VettSequenza_Punto.push_back({i});
            }

        }
    }
    // FINE CONFRONTI SEQUENZE e RAGGRUPPAMENTI IN POLIGONI
    // INIZIO ORDINAMENTO LATI e SALVATAGGIO IN CELL2D
    size_t numSottopoligoni = LinktraSequenzaELista.size(); // ogni sottopoligono è univocamente determinato da una sequenza: numSottoPol = numSequenze
    sottoPoligono.NumberCell2D = numSottopoligoni;


    sottoPoligono.Cell2DEdges.resize(numSottopoligoni);
    sottoPoligono.Cell2DVertices.resize(numSottopoligoni);
    sottoPoligono.Cell1DVertices.reserve(numSottopoligoni);
    sottoPoligono.Cell1DId.reserve(numSottopoligoni*4);

    for (unsigned int i = 0; i < numSottopoligoni; i++)
    {
        list<unsigned int> listaIdVertici = VettSequenza_Punto[i];
        // dò in pasto a funzione di Ceci per ordinamento
        Creo_sottopoligono(z, i, listaIdVertici, sottoPoligono, fracture);
    }


    // controllo che stampi tutto bene

    for(unsigned int j = 0; j < sottoPoligono.Cell1DId.size(); j++){ // giusto
        cout << "gli estremi del lato con id " << sottoPoligono.Cell1DId[j] << " sono: " << sottoPoligono.Cell1DVertices[j][0] << " e " << sottoPoligono.Cell1DVertices[j][1] << endl;
    }

    cout << "NumberCell1D: " << sottoPoligono.NumberCell1D << endl; // giusto

    cout << endl;
    // verifico gli identificativi dei lati del primo sottopoligono
    //cout << "num lati per secondo sottopoligono: " << sottoPoligono.NumberEdges.back() << endl;
    for(unsigned int i=0; i < 4; i++){
        cout << "lato " << i << " : id "<< sottoPoligono.Cell2DEdges[1][i] << endl; // SBAGLIATO
    }

    //verifico che sia giusto il riempimento di Cell2DVertices -> GIUSTO
    cout << "id degli estremi del sottopoligono 0" << endl;
    for(unsigned int i = 0; i < 4; i++){
        cout << "estremo " << i << " : " << sottoPoligono.Cell2DVertices[0][i] << endl;
    }
    cout << "id degli estremi del sottopoligono 1" << endl;
    for(unsigned int i = 0; i < 4; i++){
        cout << "estremo " << i << " : " << sottoPoligono.Cell2DVertices[1][i] << endl;
    }

    // attenzione: ne viene uno in + e gli ultimi sono sballati
    cout << "Cell1DVertices[0][0]: " << sottoPoligono.Cell1DVertices[0][0] << ", Cell1DVertices[0][1]: " << sottoPoligono.Cell1DVertices[0][1] << endl;
    cout << "Cell1DVertices[1][0]: " << sottoPoligono.Cell1DVertices[1][0] << ", Cell1DVertices[1][1]: " << sottoPoligono.Cell1DVertices[1][1] << endl;
    cout << "Cell1DVertices[2][0]: " << sottoPoligono.Cell1DVertices[2][0] << ", Cell1DVertices[2][1]: " << sottoPoligono.Cell1DVertices[2][1] << endl;
    cout << "Cell1DVertices[3][0]: " << sottoPoligono.Cell1DVertices[3][0] << ", Cell1DVertices[3][1]: " << sottoPoligono.Cell1DVertices[3][1] << endl;
    cout << "Cell1DVertices[4][0]: " << sottoPoligono.Cell1DVertices[4][0] << ", Cell1DVertices[4][1]: " << sottoPoligono.Cell1DVertices[4][1] << endl;
    cout << "Cell1DVertices[5][0]: " << sottoPoligono.Cell1DVertices[5][0] << ", Cell1DVertices[5][1]: " << sottoPoligono.Cell1DVertices[5][1] << endl;
    cout << "Cell1DVertices[6][0]: " << sottoPoligono.Cell1DVertices[6][0] << ", Cell1DVertices[6][1]: " << sottoPoligono.Cell1DVertices[6][1] << endl;
    cout << endl;
    cout << "Cell1DId[0]" << sottoPoligono.Cell1DId[0] << endl;
    cout << "Cell1DId[1]" << sottoPoligono.Cell1DId[1] << endl;
    cout << "Cell1DId[2]" << sottoPoligono.Cell1DId[2] << endl;
    cout << "Cell1DId[3]" << sottoPoligono.Cell1DId[3] << endl;
    cout << "Cell1DId[4]" << sottoPoligono.Cell1DId[4] << endl;
    cout << "Cell1DId[5]" << sottoPoligono.Cell1DId[5] << endl;
    cout << "Cell1DId[6]" << sottoPoligono.Cell1DId[6] << endl;
    // FINE ORDINAMENTO LATI e SALVATAGGIO IN CELL2D (già fatto tutto in Creo_sottopoligono)




    return 0;
}
