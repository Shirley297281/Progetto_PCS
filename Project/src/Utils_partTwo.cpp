#include "Utils.hpp"
#include "FracturesTracesPolygons.hpp"
#include "inline.hpp"
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

namespace Vettore{
bool operator==(const VectorXd &v1, const VectorXd &v2) // operatore per uguaglianza tra vettori in Eigen
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


void GeometryLibrary::Polygons::GedimInterface(vector<vector<unsigned int>>& triangles,
                                               VectorXi& materials)
{
    // funzione di triangolazione
    const unsigned int numPolygons = NumberCell2D;
    vector<vector<vector<unsigned int>>> triangleList(numPolygons);

    for (unsigned int p = 0; p < numPolygons; p++)
    {

        const unsigned int numPolygonVertices = Cell2DVertices[p].size();

        for (unsigned int v = 0; v < numPolygonVertices; v++)
        {
            const unsigned int nextVertex = Cell2DVertices[p][(v + 1) % numPolygonVertices];
            const unsigned int nextNextVertex = Cell2DVertices[p][(v + 2) % numPolygonVertices];

            if ((v + 2) % numPolygonVertices == 0)
                break;

            vector<unsigned int> triangle_vertices = {Cell2DVertices[p][0], nextVertex, nextNextVertex};

            triangleList[p].push_back(triangle_vertices);
        }

    }

    //effettiva funzione GedimInterface
    unsigned int numTotalTriangles = 0;
    for (unsigned int p = 0; p < numPolygons; p++)
        numTotalTriangles += triangleList[p].size();

    triangles.reserve(numTotalTriangles);
    materials = VectorXi::Zero(numTotalTriangles);

    unsigned int count = 0;
    for (unsigned int p = 0; p < numPolygons; p++)
    {
        for (unsigned int t = 0; t < triangleList[p].size(); t++)
        {
            triangles.push_back(triangleList[p][t]);
            materials(count) = p;
            count++;
        }
    }
}



namespace GeometryLibrary{
void MemorizzaVerticiPassanti_Cell0Ds(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z){

    // salvo i punti coincidenti con i vertici della frattura
    MatrixXd insiemeVerticiFrattZ = fracture.CoordinatesVertice[z];//estraggo la matrice di vertici di z

    unsigned int numVerticiFrattZ = fracture.numVertices[z];
    map<unsigned int, list<unsigned int>>& markerDiz = sottoPoligono.Cell0DMarkers; // passo in referenza così lo modifico direttamente senza farne una copia
    cout << "Numero di vertici della frattura: " << numVerticiFrattZ << endl;

    //inserimento dei vertici della frattura nella mesh
    for (unsigned int i = 0; i<numVerticiFrattZ; i++) //ciclo sui vertici
    {
        Vector3d vertice = insiemeVerticiFrattZ.col(i);

        unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
        sottoPoligono.NumberCell0D = NumPuntiFinora  + 1;

        sottoPoligono.Cell0DId.push_back(NumPuntiFinora); // d'ora in poi userò come id dei vertici il numero di punti salvati
        sottoPoligono.Cell0DCoordinates.push_back(vertice);
        MatrixXd M(0,0);
        sottoPoligono.SequenzeXpunto.push_back(M);
        markerDiz[0].push_back(NumPuntiFinora); // marker con chiave 0 se vertici
    }

    // salvo puntiestremi di tracce passanti traccia
    /// ciclo su fratture
    unsigned int numTraccePassantiInZ = trace.TraceIdsPassxFracture[z].size();
    for (unsigned int i = 0; i<numTraccePassantiInZ ; i++ ){


        unsigned int idTraccia = trace.TraceIdsPassxFracture[z][i];

        for (int t= 0; t<2; t++){
        Vector3d EstremoTraccia = trace.CoordinatesEstremiTraces[idTraccia].col(t); //estraggo estremo 1 della traccia i

        vector<Vector3d> VettoreCoordinateIn0D = sottoPoligono.Cell0DCoordinates;
        // CheckInserimento returna true se l'inserimento non è ancora avvenuto e false se è già avvenuto (inline)
        if (checkInserimento(EstremoTraccia, VettoreCoordinateIn0D)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
        {
            continue;
        }

        cout<<"aggiungo estremi di traccia "<<idTraccia<<endl;
        addAndPrintPoint(sottoPoligono, markerDiz, EstremoTraccia, 1);

        }
        cout<<"\n\nmarker aggiornati dopo inserimento estremi di traccia: "<<endl;
        for (const auto& pair : markerDiz) {
            cout << "Marker " << pair.first << ": ";
            for (const auto& id : pair.second) {
                cout << id << " ";
            }
            cout << endl;
        }

        // salvo punti che derivano da intersezione di due tracce
        // ciclo sulle fratture che rimangono (prendo sempre quelle successive ma prima controllo di non sforare con l'iteratore)
        if (i == trace.TraceIdsPassxFracture[z].size()-1)
        {
            continue;
        }
        for (unsigned int j = i + 1; j< trace.TraceIdsPassxFracture[z].size()-1; j++ )
        {
            unsigned int idTraccia2 = trace.TraceIdsPassxFracture[z][j];
            Vector3d Estremo1Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(0); //estraggo estremo 1 della traccia i
            Vector3d Estremo2Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(1); //estraggo estremo 2 della traccia i

            Vector3d Estremo1Traccia2 = trace.CoordinatesEstremiTraces[idTraccia2].col(0); //estraggo estremo 1 della traccia j
            Vector3d Estremo2Traccia2 = trace.CoordinatesEstremiTraces[idTraccia2].col(1); //estraggo estremo 2 della traccia j
            // calcolo intersezione retta su cui giace il lato e retta su cui giace la traccia
            Vector3d Punto02 = {0, 0, 0};
            Vector3d t = Estremo2Traccia - Estremo1Traccia;
            if (!intersezione_rette(t, Estremo1Traccia2, Estremo2Traccia2,Estremo1Traccia, Punto02))
            {
                continue;
            }
            // controllo se il punto è compreso tra gli estremi di almeno una delle due traccia
            //(basta perchè i punti di una traccia per definizione appartengono alla frattura)
            if (!combinazione_convessa(Estremo1Traccia, Estremo2Traccia, Punto02))
            {
                continue;
            }
            // controllo che punto non sia già stato inserito nelle strutture (in caso di intersezioni coincidenti)
            vector<Vector3d> VettoreCoordinateIn0D = sottoPoligono.Cell0DCoordinates;
            if (checkInserimento(Punto02, VettoreCoordinateIn0D)) // se esiste già lo stesso punto in Cell0D (=false) non aggiungerlo
            {
                continue;
            }
            cout<<"C'è intersezione tra traccia "<<idTraccia<< " e "<<idTraccia2<<" interna alla frattura. "<<endl;
            // a questo punto inserisco il punto aggiornando le varie strutture dati
            addAndPrintPoint(sottoPoligono, markerDiz, Punto02, 2);
        }


    }

    cout<<"\n\nmarker aggiornati alla fine del controllo sulle passanti: "<<endl;
    for (const auto& pair : markerDiz) {
        cout << "Marker " << pair.first << ": ";
        for (const auto& id : pair.second) {
            cout << id << " ";
        }
        cout << endl;
    }
    cout <<endl;
}

void MemorizzaVerticiNonPassanti_Cell0Ds (const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z, vector<Matrix<double, 3, 4>>& NuoviEstremi)
{
    MatrixXd insiemeVerticiFrattZ = fracture.CoordinatesVertice[z];
    unsigned int numVerticiFrattZ = fracture.numVertices[z];
    unsigned int numTracceNoPassantiInZ = trace.TraceIdsNoPassxFracture[z].size();
    unsigned int numTraccePassantiInZtrace = trace.TraceIdsPassxFracture[z].size();
    /// ciclo per le tracce non passanti
    for (unsigned int i = 0; i<numTracceNoPassantiInZ ; i++ )
    {
        // 1) calcolo tutti i punti "papabili" nuovi estremi della traccia non passante
        cout << "Traccia non passante: " << i << std::endl;
        vector<Vector3d> puntiIntersPapabili = {}; //salvo i punti di intersezione tra traccia non passante i e lati/tracce
        puntiIntersPapabili.reserve(numTraccePassantiInZtrace + numVerticiFrattZ + numTracceNoPassantiInZ);
        vector<Vector3d> vettoriDirettoriRette = {}; // salvo vettori direttori rette per ogni punto in puntiIntersPapabili
        vettoriDirettoriRette.reserve(numTraccePassantiInZtrace + numVerticiFrattZ);
        // in puntiIntersPapabili e vettoriDirettoriRette sfrutto la posizione vettore direttore di
        // puntiIntersPapabili[j] = vettoriDirettoriRette[j] evitando la creazione di dizionari
        unsigned int idTraccia = trace.TraceIdsNoPassxFracture[z][i];
        Vector3d Estremo1Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(0); //estraggo estremo 1 della traccia i
        Vector3d Estremo2Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(1); //estraggo estremo 2 della traccia i
        Vector3d t = Estremo2Traccia - Estremo1Traccia;
        cout << "Estremi della traccia: (" << Estremo1Traccia.transpose() << "), (" << Estremo2Traccia.transpose() << ")" << endl;
        /// ciclo su vertici e determino intersezione lati-traccia non passante i-esima
        for (unsigned int j = 0; j < numVerticiFrattZ ; j++)
        {
            // estraggo un vertice della frattura e il successivo per individuare un lato della frattura
            Vector3d Vertice1 = insiemeVerticiFrattZ.col(j);
            Vector3d Vertice2 = {0,0,0};
            if (j == numVerticiFrattZ -1)
            {
                Vertice2 = insiemeVerticiFrattZ.col(0);
            }
            else
            {
                Vertice2 = insiemeVerticiFrattZ.col(j+1);
            }

            // calcolo intersezione retta su cui giace il lato e retta su cui giace la traccia
            Vector3d Punto0 = {0, 0, 0};
            if (!intersezione_rette(t, Vertice1,Vertice2,Estremo1Traccia, Punto0))
            {
                continue;
            }
            // controllo se il punto è compreso tra gli estremi del lato in modo da
            // individuare solo punti all'interno della frattura
            if (!combinazione_convessa(Vertice1, Vertice2, Punto0))
            {
                continue;
            }
            puntiIntersPapabili.push_back(Punto0); // superati i controlli inserisco il punto di intersezione tra i papabili
            Vector3d vettDir = Vertice2-Vertice1;
            vettoriDirettoriRette.push_back(vettDir);
        }
        cout << "Numero di punti di intersezione dopo i lati della frattura: " << puntiIntersPapabili.size() << endl;
        /// ciclo su tracce passanti e determino intersezione traccia passante - traccia non passante i -esima
        for (unsigned int j = 0; j<  numTraccePassantiInZtrace ; j++ )
        {
            unsigned int idTraccia2 = trace.TraceIdsPassxFracture[z][j];
            Vector3d Estremo1Traccia2 = trace.CoordinatesEstremiTraces[idTraccia2].col(0); //estraggo estremo 1 della traccia j
            Vector3d Estremo2Traccia2 = trace.CoordinatesEstremiTraces[idTraccia2].col(1); //estraggo estremo 2 della traccia j
            // calcolo intersezione retta su cui giace la traccia non passante e retta su cui giace la traccia passante
            Vector3d Punto02 = {0, 0, 0};
            if (!intersezione_rette(t, Estremo1Traccia2, Estremo2Traccia2,Estremo1Traccia, Punto02))
            {
                continue;
            }
            // controllo se il punto è compreso tra gli estremi di almeno una delle due traccia
            // (basta perchè gli estremi di una traccia per def appartengono alla frattura)
            if (!combinazione_convessa(Estremo1Traccia2, Estremo2Traccia2, Punto02))
            {
                continue;
            }
            puntiIntersPapabili.push_back(Punto02); // superati i controlli inserisco il punto di intersezione tra i papabili
            Vector3d vettDir = Estremo1Traccia2-Estremo2Traccia2;
            vettoriDirettoriRette.push_back(vettDir);
        }
        cout << "Numero di punti di intersezione dopo il controllo con le tracce passanti: " << puntiIntersPapabili.size() << endl;
        /// ciclo su tracce non passanti già iterate (quindi di lunghezza maggiore e contenute in NuoviEstremi)
        for (unsigned int j = 0; j< NuoviEstremi.size() ; j++ )
        {

            Vector3d Estremo1Traccia3 = NuoviEstremi[j].col(0);
            Vector3d Estremo2Traccia3 = NuoviEstremi[j].col(2);
            // calcolo intersezione retta su cui giace la traccia non passante e retta su cui giace la traccia non passante i-esima
            Vector3d Punto03 = {0, 0, 0};
            if (!intersezione_rette(t, Estremo1Traccia3, Estremo2Traccia3,Estremo1Traccia, Punto03))
            {
                continue;
            }
            // il punto di intersezione deve appartenere agli estremi della traccia non passante che sto iterando
            if (!combinazione_convessa(Estremo2Traccia3, Estremo1Traccia3, Punto03))
            {
                continue;
            }
            puntiIntersPapabili.push_back(Punto03); // superati i controlli inserisco il punto di intersezione tra i papabili
            Vector3d vettDir = Estremo1Traccia3-Estremo2Traccia3;
            vettoriDirettoriRette.push_back(vettDir);
        }
        cout << "Numero di punti di intersezione dopo le tracce non passanti: " << puntiIntersPapabili.size() << endl;
        // 2) valutazione punti papabili
        // ora valuto tra i punti papabili quali sono i punti che mi interessano:
        // ci saranno due punti che saranno deputati all'essere estremi della traccia non passante
        // (quello più grande tra il minimo e il più piccolo tra il massimo se guardo i parametri liberi)
        // e abbiamo quindi allungato la traccia non passante.
        // ci saranno altri punti (interni al segmento della traccia allungata che non saranno estremi
        // ma andranno a inseriti nella mesh (in Cell0D) in quanto sono visibili.
        vector<array<double,2>> PerEstremoSinistro = {};
        vector<array<double,2>> PerEstremoDestro = {};
        PerEstremoDestro.reserve(puntiIntersPapabili.size());
        PerEstremoSinistro.reserve(puntiIntersPapabili.size());
        // osserviamo che per come defininiamo in seguito la combinazione convessa convenzionalmente avremo:
        double alphaEstTr1 = 1;
        double alphaEstTr2 = 0;
        /// ciclo sui punti papabili e calcolo il parametro libero considerando come retta quella della traccia non passante
        for (unsigned int k = 0; k <  puntiIntersPapabili.size() ; k++ )
        {
            Vector3d punto = puntiIntersPapabili[k];
            double alpha = 0;
            // inizio calcolo parametro libero punto
            for (int s=0;s<3;s++)
            {
                if (abs(Estremo1Traccia[s]-Estremo2Traccia[s])>tolDefault)
                {
                    alpha = (punto[s] - Estremo2Traccia[s])/(Estremo1Traccia[s]-Estremo2Traccia[s]);
                    break;
                }
            }
            // fine calcolo parametro libero punto
            /// valutazione tipo punto: interno alla traccia non passante o no
            if (alpha <= alphaEstTr2)
            //controllo se alpha è minore dell'alpha minore tra i due alpha dei due estremi e raccolgo i papabili estremi da un lato
            {
                array<double,2> arraino = {double(k),alpha}; //trasformo k in double per essere conforme a bubbleSort
                PerEstremoSinistro.push_back(arraino);
            }
            else if(alpha >= alphaEstTr1)
            //controllo se alpha è maggiore dell'alpha maggiore tra i due alpha dei due estremi e raccolgo i papabili estremi dall'altro lato
            {
                array<double,2> arraino = {double(k),alpha}; //trasformo k in double per essere conforme a bubbleSort
                PerEstremoDestro.push_back(arraino);
            }
            else
            {   /// punto di intersezione in mezzo agli estremi "vecchi" della traccia non passante va già inserito
                if (!checkInserimento(punto, sottoPoligono.Cell0DCoordinates)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
                {
                    continue;
                }

                addAndPrintPoint(sottoPoligono, sottoPoligono.Cell0DMarkers, punto,3);
            }
        }
        /// trovo i nuovi estremi della traccia non passante
        BubbleSort_mod(PerEstremoDestro); // sto ordinando in base al parametro libero
        BubbleSort_mod(PerEstremoSinistro); // sto ordinando in base al parametro libero
        // sfrutto l'ordinamento per prendere il massimo e il minimo e allo stesso tempo portare dietro l'informazione su k
        array<double, 2> UltimoArray = PerEstremoSinistro.back();
        array<double, 2> PrimoArray = PerEstremoDestro.front(); // devo prendere il minimo tra i parametri liberi maggiori = al parametro libero dell'estremo della traccia non passante
        // questo lo facciamo perchè gli estremi "nuovi" della traccia non passante devono essere i punti di intersezione più vicini agli estremi originali
        double pos1Double = UltimoArray[0]; // estraggo k in modo da riuscire a estrarre nuovamente i punti "buoni"
        double pos2Double = PrimoArray[0];
        unsigned int pos1 = static_cast<unsigned int>(pos1Double);  // pos1 e pos2 sono due interi ora
        unsigned int pos2 = static_cast<unsigned int>(pos2Double);
        // inizio salvataggio in Cell0D
        Vector3d nuovoEstremo1 = puntiIntersPapabili[pos1];
        if (!checkInserimento(nuovoEstremo1, sottoPoligono.Cell0DCoordinates)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
        {
            addAndPrintPoint(sottoPoligono, sottoPoligono.Cell0DMarkers, nuovoEstremo1,3);
        }
        Vector3d nuovoEstremo2 = puntiIntersPapabili[pos2];
        if (!checkInserimento(nuovoEstremo2, sottoPoligono.Cell0DCoordinates)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
        {
            addAndPrintPoint(sottoPoligono, sottoPoligono.Cell0DMarkers, nuovoEstremo2,3);
        }
        /// inizio salvataggio in "NuoviEstremi"
        Matrix<double, 3, 4> matrice;
        Vector3d primaColonna = nuovoEstremo1;
        Vector3d secondaColonna = vettoriDirettoriRette[pos1];
        Vector3d terzaColonna = nuovoEstremo2;
        Vector3d quartaColonna = vettoriDirettoriRette[pos2];
        matrice.col(0) = primaColonna;
        matrice.col(1) = secondaColonna;
        matrice.col(2) = terzaColonna;
        matrice.col(3) = quartaColonna;
        NuoviEstremi.push_back(matrice);
        // in precedenza ho sfruttato il fatto che i pushback in puntiIntersPapabili e vettoriDirettoriRette sono stati contemporanei
        // ho quindi acceduto alle informazioni corrispondenti per posizione risparmiando l'utilizzo di un dizionario
    }
}

void Creazioni_Sequenze_Passanti(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z)
{
    Vector3d vecNormaleAfratt = fracture.vettoreNormalePiano[z];
    /// ciclo su tracce passanti
    for (unsigned int i = 0; i<trace.TraceIdsPassxFracture[z].size() ; i++ )
    {
        unsigned int idTraccia = trace.TraceIdsPassxFracture[z][i];
        Vector3d Estremo1Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(0); //estraggo estremo 1 della traccia i
        Vector3d Estremo2Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(1); //estraggo estremo 2 della traccia i
        // controllo per ogni punto se è a "destra" o "sinistra" della retta che individua la traccia
        /// ciclo su tutti i punti salvati in precedenza in Cell0D
        for(unsigned int j = 0; j < sottoPoligono.NumberCell0D; j++)
        {
            //estraggo la matrice riferita al punto con id = j dal vettore per riferimento: modificando M modifico quella nel vettore
            // evitando copie
            MatrixXd& M = sottoPoligono.SequenzeXpunto[j];
            Vector3d coordinatePuntoInCell0d = sottoPoligono.Cell0DCoordinates[j];
            Vector3d vec1 = Estremo2Traccia - Estremo1Traccia;
            Vector3d vec2 = coordinatePuntoInCell0d - Estremo1Traccia;
            Vector3d prodVett = vec1.cross(vec2);
            double prodScal = prodVett.dot(vecNormaleAfratt);
            if (abs(prodScal)< tolDefault) //prodScal = 0 se e solo se prodVett = 0 se e solo se punto appartiene alla traccia
            {
                // duplico le sequenze e assegno sia 0 che 1
                if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                {
                    MatrixXd MatriceDiSupporto(1, 2);
                    MatriceDiSupporto.row(0) << 1, 0;
                    M = MatriceDiSupporto;
                }
                else // se la matrice non è vuota dovrò duplicare le sequenza già esistenti e inserire vettori di 0 in una e vettori di 1 nell'altra
                {
                    MatrixXd MatriceDiSupporto1(M.rows() + 1, M.cols());
                    MatrixXd MatriceDiSupporto2(M.rows() + 1, M.cols());
                    RowVectorXd nuovaRiga0 = RowVectorXd::Zero(M.cols());
                    RowVectorXd nuovaRiga1 = RowVectorXd::Ones(M.cols());
                    MatriceDiSupporto1 << M, nuovaRiga0;
                    MatriceDiSupporto2 << M, nuovaRiga1;
                    MatrixXd MConcatenata(MatriceDiSupporto1.rows(), MatriceDiSupporto1.cols() + MatriceDiSupporto2.cols());
                    MConcatenata << MatriceDiSupporto1, MatriceDiSupporto2;
                    M = MConcatenata;
                }
            }
            else if( prodScal < 0)
            {
                // assegno 1 alla sequenza (convenzione)
                unsigned int numCols;
                if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                {
                    M = MatrixXd:: Ones(1, 1);
                }
                else
                {
                    numCols = M.cols();
                    RowVectorXd nuovaRiga = RowVectorXd::Ones(numCols); // creo vettore riga di tutti 1
                    MatrixXd MatriceDiSupporto(M.rows() + 1, numCols);

                    MatriceDiSupporto << M, nuovaRiga; // copio A e aggiungo vettore di 1

                    // Aggiorniamo la matrice A e quindi anche quella nel vettore di matrici
                    M = MatriceDiSupporto;
                }

            }
            else
            {
                // assegno 0 alla sequenza (convenzione)
                unsigned int numCols;
                if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                {
                    M = MatrixXd:: Zero(1, 1);
                }
                else
                {
                    numCols = M.cols();
                    RowVectorXd nuovaRiga = RowVectorXd::Zero(numCols); // creo vettore riga di tutti 0
                    MatrixXd MatriceDiSupporto(M.rows() + 1, numCols);

                    MatriceDiSupporto << M, nuovaRiga; // copio A e aggiungo vettore di 0

                    // Aggiorniamo la matrice A e quindi anche quella nel vettore di matrici
                    M = MatriceDiSupporto;
                }
            }
        }
    }
}


void Creazioni_Sequenze_NONPassanti(const Fractures& fracture, Polygons& sottoPoligono, unsigned int z,  vector<Matrix<double, 3, 4>>& NuoviEstremi) {
    Vector3d vecNormaleAfratt = fracture.vettoreNormalePiano[z]; // anche qesta variabile l'abbiamo usata in altri posti

    for (unsigned int i = 0; i<NuoviEstremi.size() ; i++ ) //ciclo sulle tracce non passanti
    {
        Vector3d Estremo1Traccia = NuoviEstremi[i].col(0);
        Vector3d VettDirTraccia1 = NuoviEstremi[i].col(1);
        Vector3d Estremo2Traccia = NuoviEstremi[i].col(2);
        Vector3d VettDirTraccia2= NuoviEstremi[i].col(3);

        // controllo per ogni punto se è a destra o sinistra
        /// ciclo sui punti in Cell0D
        for(unsigned int j = 0; j < sottoPoligono.NumberCell0D; j++)
        {
            // non riesco a spiegarla sul codice questa parte
            MatrixXd& M = sottoPoligono.SequenzeXpunto[j]; //estraggo la matrice riferita al punto con id = j dal vettore per riferimento: modificando M modifico quella nel vettore
            Vector3d coordinatePuntoInCell0d = sottoPoligono.Cell0DCoordinates[j];
            Vector3d vec1 = Estremo2Traccia - Estremo1Traccia;
            Vector3d vec2 = coordinatePuntoInCell0d -Estremo1Traccia;
            Vector3d prodVett1 = VettDirTraccia1.cross(vec1);
            Vector3d prodVett2 = VettDirTraccia1.cross(vec2);
            double prodScal1 = prodVett1.dot(vecNormaleAfratt);
            double prodScal2 = prodVett2.dot(vecNormaleAfratt);
            Vector3d vec3 = Estremo1Traccia - Estremo2Traccia;
            Vector3d vec4 = coordinatePuntoInCell0d -Estremo2Traccia;
            Vector3d prodVett3 = VettDirTraccia2.cross(vec3);
            Vector3d prodVett4 = VettDirTraccia2.cross(vec4);
            double prodScal3 = prodVett3.dot(vecNormaleAfratt);
            double prodScal4 = prodVett4.dot(vecNormaleAfratt);

            if (abs(prodScal2) < tolDefault || abs(prodScal4) < tolDefault)
            // 1) ovvero se il punto appartiene a una delle due rette individuate dai vettori direttori
            {
                // L'assegnazione della sequenza dipende dalla casistica
                Vector3d prodVett = vec1.cross(vec2);
                double prodScal = prodVett.dot(vecNormaleAfratt);
                //1.1) punto appartiene anche alla traccia non passante
                if (abs(prodScal)< 1e-14) //prodScal = 0 se e solo se prodVett = 0 se e solo se punto appartiene alla traccia non passante
                // sarà in particolare un estremo della traccia non passante
                {
                    // duplico le sequenze e assegno sia 0,1
                    if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                    {
                        MatrixXd MatriceDiSupporto(1, 2);
                        MatriceDiSupporto.row(0) << 1, 0;
                        M = MatriceDiSupporto;
                    }
                    else // se la matrice non è vuota dovrò duplicare le sequenza già esistenti e inserire vettori di 0 in una e vettori di 1 nell'altra
                    {
                        MatrixXd MatriceDiSupporto1(M.rows() + 1, M.cols());
                        MatrixXd MatriceDiSupporto2(M.rows() + 1, M.cols());
                        RowVectorXd nuovaRiga0 = RowVectorXd::Zero(M.cols());
                        RowVectorXd nuovaRiga1 = RowVectorXd::Ones(M.cols());
                        MatriceDiSupporto1 << M, nuovaRiga0;
                        MatriceDiSupporto2 << M, nuovaRiga1;
                        MatrixXd MConcatenata(MatriceDiSupporto1.rows(), MatriceDiSupporto1.cols() + MatriceDiSupporto2.cols());
                        MConcatenata << MatriceDiSupporto1, MatriceDiSupporto2;
                        M = MConcatenata;
                    }
                }
                else if( prodScal < 0)
                //1.2) ovvero il punto è da un lato rispetto alla traccia non passante
                {
                    // duplico le sequenze e assegno sia 1 che 2
                    if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                    {
                        MatrixXd MatriceDiSupporto(1, 2);
                        MatriceDiSupporto.row(0) << 2, 1;
                        M = MatriceDiSupporto;
                    }
                    else // se la matrice non è vuota dovrò duplicare le sequenza già esistenti e inserire vettori di 0 in una e vettori di 1 nell'altra
                    {
                        MatrixXd MatriceDiSupporto1(M.rows() + 1, M.cols());
                        MatrixXd MatriceDiSupporto2(M.rows() + 1, M.cols());
                        RowVectorXd nuovaRiga0 = RowVectorXd::Ones(M.cols())*2;
                        RowVectorXd nuovaRiga1 = RowVectorXd::Ones(M.cols());
                        MatriceDiSupporto1 << M, nuovaRiga0;
                        MatriceDiSupporto2 << M, nuovaRiga1;
                        MatrixXd MConcatenata(MatriceDiSupporto1.rows(), MatriceDiSupporto1.cols() + MatriceDiSupporto2.cols());
                        MConcatenata << MatriceDiSupporto1, MatriceDiSupporto2;
                        M = MConcatenata;
                    }

                }
                else
                // 1.3) ovvero il punto è dall'altro lato rispetto alla traccia non passante
                {
                    // duplico le sequenze e assegno sia 0 che 2
                    if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                    {
                        MatrixXd MatriceDiSupporto(1, 2);
                        MatriceDiSupporto.row(0) << 2, 0;
                        M = MatriceDiSupporto;
                    }
                    else // se la matrice non è vuota dovrò duplicare le sequenza già esistenti e inserire vettori di 0 in una e vettori di 1 nell'altra
                    {
                        MatrixXd MatriceDiSupporto1(M.rows() + 1, M.cols());
                        MatrixXd MatriceDiSupporto2(M.rows() + 1, M.cols());
                        RowVectorXd nuovaRiga0 = RowVectorXd::Ones(M.cols())*2;
                        RowVectorXd nuovaRiga1 = RowVectorXd::Zero(M.cols());
                        MatriceDiSupporto1 << M, nuovaRiga0;
                        MatriceDiSupporto2 << M, nuovaRiga1;
                        MatrixXd MConcatenata(MatriceDiSupporto1.rows(), MatriceDiSupporto1.cols() + MatriceDiSupporto2.cols());
                        MConcatenata << MatriceDiSupporto1, MatriceDiSupporto2;
                        M = MConcatenata;
                    }
                }
            }
            //2) il punto è fuori dall'area di influenza della traccia non passante
            else if (prodScal1*prodScal2<0 || prodScal3*prodScal4<0 )

            {
                // assegno 2 alla sequenza (convenzione)
                unsigned int numCols;
                if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                {
                    M = MatrixXd:: Ones(1, 1)*2;
                }
                else
                {
                    numCols = M.cols();
                    RowVectorXd nuovaRiga = RowVectorXd::Ones(numCols)*2; // creo vettore riga di tutti 2
                    MatrixXd MatriceDiSupporto(M.rows() + 1, numCols);

                    MatriceDiSupporto << M, nuovaRiga; // copio A e aggiungo vettore di 2

                    // Aggiorniamo la matrice A e quindi anche quella nel vettore di matrici
                    M = MatriceDiSupporto;
                }
            }
            // 3) siamo interni all'area di influenza della traccia non passante
            else
            {
                // assegno 0 o 1 o entrambi sdoppiando come già fatto nel caso passante
                Vector3d prodVett = vec1.cross(vec2);
                double prodScal = prodVett.dot(vecNormaleAfratt);
                if (abs(prodScal)< 1e-14) //prodScal = 0 se e solo se prodVett = 0 se e solo se punto appartiene alla traccia
                {
                    // duplico le sequenze e assegno sia 0 che 1
                    if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                    {
                        MatrixXd MatriceDiSupporto(1, 2);
                        MatriceDiSupporto.row(0) << 1, 0;
                        M = MatriceDiSupporto;
                    }
                    else // se la matrice non è vuota dovrò duplicare le sequenza già esistenti e inserire vettori di 0 in una e vettori di 1 nell'altra
                    {
                        MatrixXd MatriceDiSupporto1(M.rows() + 1, M.cols());
                        MatrixXd MatriceDiSupporto2(M.rows() + 1, M.cols());
                        RowVectorXd nuovaRiga0 = RowVectorXd::Zero(M.cols());
                        RowVectorXd nuovaRiga1 = RowVectorXd::Ones(M.cols());
                        MatriceDiSupporto1 << M, nuovaRiga0;
                        MatriceDiSupporto2 << M, nuovaRiga1;
                        MatrixXd MConcatenata(MatriceDiSupporto1.rows(), MatriceDiSupporto1.cols() + MatriceDiSupporto2.cols());
                        MConcatenata << MatriceDiSupporto1, MatriceDiSupporto2;
                        M = MConcatenata;
                    }
                }
                else if( prodScal < 0)
                {
                    // assegno 1 alla sequenza (convenzione)
                    unsigned int numCols;
                    if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                    {
                        M = MatrixXd:: Ones(1, 1);
                    }
                    else
                    {
                        numCols = M.cols();
                        RowVectorXd nuovaRiga = RowVectorXd::Ones(numCols); // creo vettore riga di tutti 1
                        MatrixXd MatriceDiSupporto(M.rows() + 1, numCols);

                        MatriceDiSupporto << M, nuovaRiga; // copio A e aggiungo vettore di 1

                        // Aggiorniamo la matrice A e quindi anche quella nel vettore di matrici
                        M = MatriceDiSupporto;
                    }

                }
                else
                {
                    // assegno 0 alla sequenza (convenzione)
                    unsigned int numCols;
                    if (M.cols() == 0)  //la matrice è ancora vuota: è la prima volta che "pesco" il punto
                    {
                        M = MatrixXd:: Zero(1, 1);
                    }
                    else
                    {
                        numCols = M.cols(); //conto le colonne
                        RowVectorXd nuovaRiga = RowVectorXd::Zero(numCols); // creo vettore riga di tutti 0
                        MatrixXd MatriceDiSupporto(M.rows() + 1, numCols);

                        MatriceDiSupporto << M, nuovaRiga; // copio A e aggiungo vettore di 0

                        // Aggiorniamo la matrice A e quindi anche quella nel vettore di matrici
                        M = MatriceDiSupporto;
                    }
                }
            }
            ///////////////////////////////////////////////////////////////////////////////////////////
        }
    }


}

  
void Creo_sottopoligono(unsigned int num_fracture, unsigned int num_sottopoligono,list<unsigned int> listaIdVertici, Polygons& sottopoligono, Fractures& fracture){

    vector<unsigned int> estremi(listaIdVertici.begin(), listaIdVertici.end()); // trasformo la lista in un vector
    unsigned int n = estremi.size(); // num dei vertici del sottopoligono
    vector<Vector2i> id_estremi_lato; // lati identificati dagli id degli estremi -> per Cell2DVertices
    id_estremi_lato.reserve(n); // n vertici => avrò n lati
    MatrixXd vertices(3, n); // per baricentro
    Vector3d vett_normale_frattura = fracture.vettoreNormalePiano[num_fracture];
    vector<unsigned int> id_lati; // per Cell2DEdges
    id_lati.reserve(n);


    // con i seleziono un vertice, con j il consecutivo (effettuo la verifica con k)
    unsigned int num_iterazioni = 0;
    unsigned int i = 0;
    unsigned int id_i = estremi[i];
    unsigned int id_j;

    for (unsigned int j = 0; j < n+1; j++) // devo arrivare fino a n iterazioni perchè sennò non trova l'ultimo lato
    {
        bool lato_valido = true;

        if (j==i){continue;}

        else if(num_iterazioni == n-1){ // ultima iterazione
            id_j = estremi[0];
        }

        else{
            id_j = estremi[j];

            Vector3d coord_i = sottopoligono.Cell0DCoordinates[id_i];
            Vector3d coord_j = sottopoligono.Cell0DCoordinates[id_j];
            Vector3d vec1 = coord_j - coord_i; // vettore direzione che congiunge i candidati vertici consecutivi

            // per ogni candidato lato devo verificare il prodotto vettoriale con il vettore congiungente un suo estremo con tutti gli altri punti che sono in totale n-2
            unsigned iter = 0;
            unsigned int m = n-2;
            vector<double> prodscalare;
            prodscalare.reserve(m);

            for (unsigned int k = 0; k < n; k++){
                if (k == i || k == j){
                    continue;
                }

                unsigned int id_k = estremi[k];
                Vector3d coord_k = sottopoligono.Cell0DCoordinates[id_k];
                Vector3d vec2 = coord_k - coord_i;
                Vector3d prodVett = {};
                prodVett = vec1.cross(vec2);
                prodscalare[iter] = prodVett.dot(vett_normale_frattura);


                if(iter > 0){
                    if ((prodscalare[iter] > 0 && prodscalare[iter -1] < 0) || (prodscalare[iter] < 0 && prodscalare[iter -1] > 0)){ // se è diverso dal precedente vuol dire che il lato non va bene
                        lato_valido = false;
                        continue;
                    }
                }
                iter += 1;
            }
        }

        // se sono arrivata qui => ho trovato un lato
        if(lato_valido){
            Vector2i l(id_i, id_j);
            id_estremi_lato.push_back(l);
            num_iterazioni += 1;

            // verifico se il lato è già presente in Cell1D
            auto it = find(sottopoligono.Cell1DVertices.begin(), sottopoligono.Cell1DVertices.end(), l);
            Vector2i l_inverso(id_j, id_i);
            auto it1 = find(sottopoligono.Cell1DVertices.begin(), sottopoligono.Cell1DVertices.end(), l_inverso);

            // se ho già trovato il lato per un nuovo sottopoligono => NON devo aggiornare nè Cell1DVertices nè Cell1DId MA devo fare il push_back a id_lati per aggiornare successivamente Cell2DEdges e Cell2DVertices
            if(it != sottopoligono.Cell1DVertices.end()){
                unsigned int posizione = distance(sottopoligono.Cell1DVertices.begin(), it); // rappresenta l'indice a cui it si riferisce nel vettore
                unsigned int id_ = sottopoligono.Cell1DId[posizione];
                id_lati.push_back(id_);
            }
            else if(it1 != sottopoligono.Cell1DVertices.end()){
                unsigned int posizione = distance(sottopoligono.Cell1DVertices.begin(), it1);
                unsigned int id_ = sottopoligono.Cell1DId[posizione];
                id_lati.push_back(id_);
            }
            else if(it == sottopoligono.Cell1DVertices.end() || it1 == sottopoligono.Cell1DVertices.end()){

                unsigned int id;

                if(sottopoligono.Cell1DId.empty())
                {
                    id = 0;
                }
                else if(!sottopoligono.Cell1DId.empty()){
                    id = sottopoligono.Cell1DId.back() + 1;
                }
                sottopoligono.Cell1DId.push_back(id);
                sottopoligono.Cell1DVertices.push_back(l);
                id_lati.push_back(id); // inizio a creare il vettore da inserire in Cell2DEdges (se li sto ordinando in senso orario piuttosto li inverto dopo)
            }

        }

        if(num_iterazioni <= n){
            i = j; // permette di trovare i lati in ordine

            id_i = estremi[i];// serve per andare avanti con i lati, altrimenti fa sempre riferimento al primo
        }

    }

    sottopoligono.NumberCell1D = sottopoligono.Cell1DId.size();

    // calcolo il baricentro del sottopoligono
    for (unsigned int i = 0; i < n; i++){
        unsigned int id_i = sottopoligono.Cell0DId[i];
        Vector3d coord_i = sottopoligono.Cell0DCoordinates[id_i];
        vertices.col(i) = coord_i;
    }
    array <double,3> bar = barycenter(vertices, n);
    Vector3d bar_vec(bar[0], bar[1], bar[2]);
    
    // i lati sono in ordine (devo solo verificare che siano in ordine antiorario e NON orario)
    unsigned int id_0 = id_estremi_lato[0][0];
    unsigned int id_1 = id_estremi_lato[0][1];
    Vector3d coord_0 = sottopoligono.Cell0DCoordinates[id_0];
    Vector3d coord_1 = sottopoligono.Cell0DCoordinates[id_1];

    // vettori che congiungono gli estremi del primo lato al baricentro
    Vector3d v1 = coord_0 - bar_vec;
    Vector3d v2 = coord_1 - bar_vec;
    Vector3d v1xv2 = vec_product(v1, v2);
    double prod_scal = v1xv2.dot(vett_normale_frattura);

    // se il prodotto scalare è negativo => devo prendere l'altro senso
    if(prod_scal < 0){
        reverse(id_estremi_lato.begin(), id_estremi_lato.end());
        // invertire l'ordine all'interno di ogni coppia in id_estremi_lato
        for (auto& coppia : id_estremi_lato) {
            swap(coppia[0], coppia[1]);
        }
        reverse(id_lati.begin(), id_lati.end());
    }

    // aggiorno Cell2D
    sottopoligono.Cell2DId.push_back(num_sottopoligono);
    sottopoligono.NumberVertices.push_back(n);
    sottopoligono.NumberEdges.push_back(n);

    // Cell2DVertices trasformo la lista delle coppie di estremi identificativi del lato in una sequenza di punti consecutivi
    unordered_set<int> id_estremi_set;
    vector<unsigned int> id_lati_vec;

    for (auto it = id_estremi_lato.begin(); it != id_estremi_lato.end(); ++it) {
        Vector2i vec = *it;
        // Inserisci il primo elemento se non presente nel set
        if (id_estremi_set.find(vec[0]) == id_estremi_set.end()) {
            id_estremi_set.insert(vec[0]);
            id_lati_vec.push_back(vec[0]);
        }
        // Inserisci il secondo elemento se non presente nel set
        if (id_estremi_set.find(vec[1]) == id_estremi_set.end()) {
            id_estremi_set.insert(vec[1]);
            id_lati_vec.push_back(vec[1]);
        }
    }

    sottopoligono.Cell2DVertices[num_sottopoligono] = id_lati_vec;
    sottopoligono.Cell2DEdges[num_sottopoligono] = id_lati;

}

}
