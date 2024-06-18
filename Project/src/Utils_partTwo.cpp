#include "Utils.hpp"
#include "FracturesTraces.hpp"
#include "inline.hpp"
#include <vector>
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <Eigen/Dense>


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

        // visto che push_back alloca memoria e inserisce dovrebbe funzionare comunque (reserve non sapevo come farlo)
        sottoPoligono.Cell0DId.push_back(NumPuntiFinora); //d'ora in poi userò come id dei vertici il numero di punti salvati
        sottoPoligono.Cell0DCoordinates.push_back(vertice);
        MatrixXd M(0,0);
        sottoPoligono.SequenzeXpunto.push_back(M);
        markerDiz[0].push_back(NumPuntiFinora); // marker con chiave 0 se vertici

        cout << "Vertice " << sottoPoligono.Cell0DId[i] << ": " << sottoPoligono.Cell0DCoordinates[i].transpose();
        cout <<" con Marker: ";

        bool foundidpunto = false;

        // Itera attraverso la mappa
        for (const auto& pair : markerDiz) {
            const unsigned int key = pair.first;
            const list<unsigned int>& values = pair.second;

            // Cerca l'ID del punto nella lista dei valori
            if (find(values.begin(), values.end(), sottoPoligono.Cell0DId[i]) != values.end()) {
                cout << key << endl;
                foundidpunto = true;
                break; // Esci dal ciclo una volta trovato l'ID
            }
        }

        if (!foundidpunto) {
            cout << "ID del punto " << sottoPoligono.Cell0DId[i] << " non trovato nella mappa." << endl;
        }


    }


    // salvo puntiestremi di tracce passanti traccia
    /// ciclo su fratture (cambiare TraceIdsPassxFracture con quello ordinato)
    unsigned int numTraccePassantiInZ = trace.TraceIdsPassxFracture[z].size();
    for (unsigned int i = 0; i<numTraccePassantiInZ ; i++ ){

        //Vector3d vecNormaleAfratt = fracture.vettoreNormalePiano[z]; //mi serve dopo

        unsigned int idTraccia = trace.TraceIdsPassxFracture[z][i];

        for (int t= 0; t<2; t++){
        Vector3d EstremoTraccia = trace.CoordinatesEstremiTraces[idTraccia].col(t); //estraggo estremo 1 della traccia i

        vector<Vector3d> VettoreCoordinateIn0D = sottoPoligono.Cell0DCoordinates;
        // CheckInserimento returna true se l'inserimento non è ancora avvenuto e false se è già avvenuto (inline fun)
        if (!checkInserimento(EstremoTraccia, VettoreCoordinateIn0D)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
        {
            continue;
        }

        cout<<"aggiungo estremi di traccia "<<idTraccia<<endl;
        addAndPrintPoint(sottoPoligono, markerDiz, EstremoTraccia, 1);

        }//fine for per i due estremi della traccia

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
            Vector3d t = Estremo2Traccia - Estremo1Traccia; //serve dopo
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

            cout<<"c'è intersezione tra traccia "<<idTraccia<< " e "<<idTraccia2<<" interna alla frattura. "<<endl;

            addAndPrintPoint(sottoPoligono, markerDiz, Punto02, 2);

        }


    }//fine for tracce passanti
    cout<<"\n\nmarker aggiornati alla fine del controllo sulle passanti: "<<endl;
    for (const auto& pair : markerDiz) {
        cout << "Marker " << pair.first << ": ";
        for (const auto& id : pair.second) {
            cout << id << " ";
        }
        cout << endl;
    }
    cout <<endl;


    /*
    // salvo punti che derivano da intersezione lato-traccia
    /// ciclo su fratture (cambiare TraceIdsPassxFracture con quello ordinato)
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
                Vertice2 = insiemeVerticiFrattZ.col(0);
            }
            else
            {
                Vertice2 = insiemeVerticiFrattZ.col(j+1);
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
            cout << "Punto0: " << Punto0<<endl;
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

            cout <<"++ Punto "<<sottoPoligono.Cell0DId[NumPuntiFinora]<<": "<<sottoPoligono.Cell0DCoordinates[NumPuntiFinora].transpose();

            // non so se sia necessario questo però sto allocando memoria
            MatrixXd M(0,0); //matrice vuota di dimensione
            sottoPoligono.SequenzeXpunto.push_back(M);
            markerDiz[1].push_back(NumPuntiFinora); //marker con chiave 1 per punti sui lati

            cout<<"\n\nmarker aggiornati: "<<endl;
            for (const auto& pair : markerDiz) {
                cout << "Marker " << pair.first << ": ";
                for (const auto& id : pair.second) {
                    cout << id << " ";
                }
                cout << endl;
            }
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
    }*/

}

void MemorizzaVerticiNonPassanti_Cell0Ds (const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z){
    // 1) calcolo tutti i punti "papabili" estremi
    MatrixXd insiemeVerticiFrattZ = fracture.CoordinatesVertice[z];// qui aggiungo queste due variabili anche se
    unsigned int numVerticiFrattZ = fracture.numVertices[z]; // sono comuni alle altre funzioni come MemorizzaVertici_Cello0Ds
    unsigned int numTracceNoPassantiInZ = trace.TraceIdsNoPassxFracture[z].size(); // sostituire con struttura ordinata
    unsigned int numTraccePassantiInZtrace = trace.TraceIdsPassxFracture[z].size();
    vector<Matrix<double, 3, 4>> NuoviEstremi = {};   //mi salvo i nuovi estremi delle tracce non passanti e
        // mi salvo i vettori direttori delle rette passanti per i nuovi estremi delle tracce
    NuoviEstremi.reserve(numTracceNoPassantiInZ);
    for (unsigned int i = 0; i<numTracceNoPassantiInZ ; i++ )      /// ciclo per le tracce non passanti
    {
        cout << "Traccia non passante: " << i << std::endl;
        vector<Vector3d> puntiIntersPapabili = {}; //salvo i punti di intersezione tra traccia non passante i e lati/tracce
        puntiIntersPapabili.reserve(numTraccePassantiInZtrace + numVerticiFrattZ);
        vector<Vector3d> vettoriDirettoriRette = {}; // salvo vettori direttori rette per ogni punto in puntiIntersPapabili
        vettoriDirettoriRette.reserve(numTraccePassantiInZtrace + numVerticiFrattZ);
        // in puntiIntersPapabili e vettoriDirettoriRette sfrutto la posizione vettore direttore di puntiIntersPapabili[j] = vettoriDirettoriRette[j]
        unsigned int idTraccia = trace.TraceIdsNoPassxFracture[z][i];
        Vector3d Estremo1Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(0); //estraggo estremo 1 della traccia i
        Vector3d Estremo2Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(1); //estraggo estremo 2 della traccia i
        Vector3d t = Estremo2Traccia - Estremo1Traccia;
        cout << "Estremi della traccia: (" << Estremo1Traccia.transpose() << "), (" << Estremo2Traccia.transpose() << ")" << std::endl;
        for (unsigned int j = 0; j < numVerticiFrattZ ; j++) //// ciclo su vertici e determino intersezione lati-traccia non passante
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
            if (!soluzione_sistema3x2(t, Vertice1,Vertice2,Estremo1Traccia, Punto0))
            {
                continue;
            }
            // controllo se il punto è compreso tra gli estremi della traccia in modo da
            // individuare solo punti all'interno della frattura
            // inline per combinazione convessa
            if (!combinazione_convessa(Vertice1, Vertice2, Punto0))
            {
                continue;
            }
            puntiIntersPapabili.push_back(Punto0); // superati i controlli inserisco il punto di intersezione tra i papabili
            Vector3d vettDir = Vertice2-Vertice1;
            vettoriDirettoriRette.push_back(vettDir);
        }
        cout << "Numero di punti di intersezione dopo i lati della frattura: " << puntiIntersPapabili.size() << std::endl;
        // ciclo su tracce passanti
        for (unsigned int j = 0; j<  numTraccePassantiInZtrace ; j++ ) //// ciclo su tracce passanti per determinare intersezione con traccia non passante
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
            if (!combinazione_convessa(Estremo1Traccia2, Estremo2Traccia2, Punto02)) //controllo che il punto non sia fuori dalla traccia passante in questione (altrimenti sarei fuori dalla frattura)
            {
                continue;
            }
            puntiIntersPapabili.push_back(Punto02); // superati i controlli inserisco il punto di intersezione tra i papabili
            Vector3d vettDir = Estremo1Traccia2-Estremo2Traccia2;
            vettoriDirettoriRette.push_back(vettDir);
        }
        cout << "Numero di punti di intersezione dopo le tracce passanti: " << puntiIntersPapabili.size() << std::endl;
        // ciclo su tracce non passanti già iterate
        if (NuoviEstremi.empty()) // non ci sono tracce non passanti già iterate
        {
            continue;
        }
        for (unsigned int j = 0; j< NuoviEstremi.size() ; j++ ) //// ciclo su tracce non passanti già iterate (quindi contenute in NuoviEstremi)
        {
            // vector<Matrix<double, 3, 4>> NuoviEstremi = {};
            Vector3d Estremo1Traccia3 = NuoviEstremi[j].col(0);
            Vector3d Estremo2Traccia3 = NuoviEstremi[j].col(2);
            // calcolo intersezione retta su cui giace il lato e retta su cui giace la traccia
            Vector3d Punto03 = {0, 0, 0};
            if (!soluzione_sistema3x2(t, Estremo1Traccia3, Estremo2Traccia3,Estremo1Traccia, Punto03))
            {
                continue;
            }
            if (!combinazione_convessa(Estremo2Traccia3, Estremo1Traccia3, Punto03)) // il punto di intersezione deve appartenere agli estremi della traccia non passante che sto iterando
            {
                continue;
            }
            puntiIntersPapabili.push_back(Punto03); // superati i controlli inserisco il punto di intersezione tra i papabili
            Vector3d vettDir = Estremo1Traccia3-Estremo2Traccia3;
            vettoriDirettoriRette.push_back(vettDir);
        }
        cout << "Numero di punti di intersezione dopo le tracce non passanti: " << puntiIntersPapabili.size() << std::endl;
        // 2) valutazione punti papabili
        // ora valuto tra i punti papabili quali sono i punti che mi interessano:
        // ci saranno due punti che saranno deputati all'essere estremi della traccia non passante (quello più grande tra il minimo e il più piccolo tra il massimo se guardo i parametri liberi)
        // e abbiamo quindi allungato la traccia non passante.
        // ci saranno altri punti (interni al segmento della traccia allungata che non saranno estremi ma andranno a inseriti nella mesh (in Cell0D) in quanto sono visibili.
        vector<array<double,2>> PerEstremoSinistro = {};
        vector<array<double,2>> PerEstremoDestro = {};
        PerEstremoDestro.reserve(puntiIntersPapabili.size());
        PerEstremoSinistro.reserve(puntiIntersPapabili.size());
        // osserviamo che per come defininiamo in seguito la combinazione convessa convenzionalmente avremo:
        double alphaEstTr1 = 1;
        double alphaEstTr2 = 0;
        for (unsigned int k = 0; k <  puntiIntersPapabili.size() ; k++ ) //ciclo sui punti papabili e calcolo il parametro libero considerando come retta quella della traccia non passante
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
            if (alpha <= alphaEstTr2) //controllo se alpha è minore dell'alpha minore tra i due alpha dei due estremi
            {
                array<double,2> arraino = {double(k),alpha}; //trasformo k in double per essere conforme a bubbleSort
                PerEstremoSinistro.push_back(arraino);
            }
            else if(alpha >= alphaEstTr1) //controllo se alpha è maggiore dell'alpha minore tra i due estremi
            {
                array<double,2> arraino = {double(k),alpha}; //trasformo k in double per essere conforme a bubbleSort
                PerEstremoDestro.push_back(arraino);
            }
            else
            {
                if (!checkInserimento(punto, sottoPoligono.Cell0DCoordinates)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
                {
                    continue;
                }
                // è un punto interno, lo salvo in Cell0D (senza fare altro)
                unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
                sottoPoligono.NumberCell0D = NumPuntiFinora  + 1;
                /// visto che push_back alloca memoria e inserisce dovrebbe funzionare comunque (reserve non sapevo come farlo)
                sottoPoligono.Cell0DId.push_back(NumPuntiFinora); // es: se ho già 2 punti questi hanno identificativo 0,1. Quando aggiungo il terzo questo avrà id = 2
                sottoPoligono.Cell0DCoordinates.push_back(punto);
                // non so se sia necessario questo però sto allocando memoria
                MatrixXd M(0,0); //matrice vuota di dimensione
                sottoPoligono.SequenzeXpunto.push_back(M); //pensare a un resize eventuale
                //markerDiz[2].push_back(NumPuntiFinora); //marker con chiave 2 per punti interni !marker qui non funziona, non so distinguere in maniera semplice se il punto è in una posizione particolare
            }
        }
        // esco dal ciclo sui punti
        BubbleSort_mod(PerEstremoDestro); // sto ordinando in base al parametro
        BubbleSort_mod(PerEstremoSinistro); // sto ordinando in base al parametro
        // in caso gestire il caso in cui non diamo nulla a bubbleSort_mod
        array<double, 2> UltimoArray = PerEstremoSinistro.back();         // devo prendere il massimo tra i parametri liberi minori = al parametro libero dell'estremo della traccia non passante
        array<double, 2> PrimoArray = PerEstremoDestro.front(); // devo prendere il minimo tra i parametri liberi maggiori = al parametro libero dell'estremo della traccia non passante
        // questo lo facciamo perchè gli estremi "nuovi" della traccia non passante devono essere i punti di intersezione più vicini agli estremi originali
        double pos1Double = UltimoArray[0]; // estraggo k in modo da riuscire a estrarre nuovamente i punti "buoni" da
        double pos2Double = PrimoArray[0];
        // consigliato da chat gpt per trasformare da double a int
        unsigned int pos1 = static_cast<unsigned int>(pos1Double);
        unsigned int pos2 = static_cast<unsigned int>(pos2Double);
        // inizio salvataggio in Cell0D
        Vector3d nuovoEstremo1 = puntiIntersPapabili[pos1];   // possibile problema di double - int
        if (checkInserimento(nuovoEstremo1, sottoPoligono.Cell0DCoordinates)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
        {
            unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
            sottoPoligono.NumberCell0D = NumPuntiFinora  + 1;
            /// visto che push_back alloca memoria e inserisce dovrebbe funzionare comunque (reserve non sapevo come farlo)
            sottoPoligono.Cell0DId.push_back(NumPuntiFinora); // es: se ho già 2 punti questi hanno identificativo 0,1. Quando aggiungo il terzo questo avrà id = 2
            sottoPoligono.Cell0DCoordinates.push_back(nuovoEstremo1);
            // non so se sia necessario questo però sto allocando memoria
            MatrixXd M(0,0); //matrice vuota di dimensione
            sottoPoligono.SequenzeXpunto.push_back(M);
            // markerDiz[1].push_back(NumPuntiFinora); qui per i marker è più complicato e visto che non è richiesto non mi complicherei la vita
        }
        Vector3d nuovoEstremo2 = puntiIntersPapabili[pos2];
        if (checkInserimento(nuovoEstremo2, sottoPoligono.Cell0DCoordinates)) // se esiste già lo stesso punto in Cell0D non aggiungerlo
        {
            unsigned int NumPuntiFinora = sottoPoligono.NumberCell0D;
            sottoPoligono.NumberCell0D = NumPuntiFinora  + 1;
            /// visto che push_back alloca memoria e inserisce dovrebbe funzionare comunque (reserve non sapevo come farlo)
            sottoPoligono.Cell0DId.push_back(NumPuntiFinora); // es: se ho già 2 punti questi hanno identificativo 0,1. Quando aggiungo il terzo questo avrà id = 2
            sottoPoligono.Cell0DCoordinates.push_back(nuovoEstremo2);
            // non so se sia necessario questo però sto allocando memoria
            MatrixXd M(0,0); //matrice vuota di dimensione
            sottoPoligono.SequenzeXpunto.push_back(M);
            // markerDiz[1].push_back(NumPuntiFinora); qui per i marker è più complicato e visto che non è richiesto non mi complicherei la vita
        }
        // fine salvataggio in Cell0D
        // inizio salvataggio in "NuoviEstremi"
        Matrix<double, 3, 4> matrice;
        Vector3d primaColonna = nuovoEstremo1; // sarebbe puntiIntersPapabili[pos1];
        Vector3d secondaColonna = vettoriDirettoriRette[pos1];
        Vector3d terzaColonna = nuovoEstremo2;
        Vector3d quartaColonna = vettoriDirettoriRette[pos2];
        matrice.col(0) = primaColonna;
        matrice.col(1) = secondaColonna;
        matrice.col(2) = terzaColonna;
        matrice.col(3) = quartaColonna; //si può compattare ma per capire cosa stiamo facendo va bene
        // tutto ciò dovrebbe funzionare perchè il pushback in vettoriDirettoriRette e in puntiIntersPapabili è contemporaneo, quindi le posizioni dovrebbero combaciare
        NuoviEstremi.push_back(matrice);
        // fine salvataggio in "Nuovi Estremi"
        /// parentesi fine traccia non passante salvataggio punti intersezione
    }

}//fine memo non passanti

void Creazioni_Sequenze_Passanti(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z)
{
    /// std::vector<Matrix<float, Dynamic, Dynamic>> SequenzeXpunto = {};
    /// numCol = numSequenze (=num sottopoligoni)
    /// numRow = numTraces
    Vector3d vecNormaleAfratt = fracture.vettoreNormalePiano[z]; //mi serve dopo
    for (unsigned int i = 0; i<trace.TraceIdsPassxFracture[z].size() ; i++ ) //ciclo di nuovo sulle tracce passanti
    {
        unsigned int idTraccia = trace.TraceIdsPassxFracture[z][i];
        Vector3d Estremo1Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(0); //estraggo estremo 1 della traccia i
        Vector3d Estremo2Traccia = trace.CoordinatesEstremiTraces[idTraccia].col(1); //estraggo estremo 2 della traccia i
        // controllo per ogni punto se è a destra o sinistra
        for(unsigned int j = 0; j < sottoPoligono.NumberCell0D; j++) //ciclo su tutti i punti salvati in precedenza in Cell0D
        {
            MatrixXd& M = sottoPoligono.SequenzeXpunto[j]; //estraggo la matrice riferita al punto con id = j dal vettore per riferimento: modificando M modifico quella nel vettore
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
                    numCols = M.cols(); //conto le colonne
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
    }
}

// quando la chiamo ho già fatto il controllo e ho trovato tutti i vertici con la stessa sequenza, ho incrementato il numero di Cell2D
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
    unsigned int num_iterazioni = 0; // all'ultima iterazione id_j deve essere uguale al primo id_i
    unsigned int i = 0;
    unsigned int id_i = estremi[i];
    unsigned int id_j;

    for (unsigned int j = 0; j < n+1; j++) // devo arrivare fino a n iterazioni perchè sennò non trovo l'ultimo lato
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

            for (unsigned int k = 0; k < n; k++){ // devo confrontarli con tutti gli altri punti, se avessi fatto con k = j+1 mi sarei persa il confronto con gli i e j precedenti
                if (k == i || k == j){ // se trovo k uguale a uno dei vertici che sto considerando come estremi => incremento k
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
                //sottopoligono.Cell1DVertices.push_back(l);
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
            i = j; // mi permette di trovare i lati in ordine
            id_i = estremi[i];//mi serve per andare avanti con i lati, altrimenti fa sempre riferimento al primo
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

    // i lati che trovo sono in ordine (devo solo verificare che siano in ordine antiorario e NON orario)
    // devo prendere id_estremi_lato[0]
    unsigned int id_0 = id_estremi_lato[0][0];
    unsigned int id_1 = id_estremi_lato[0][1];
    Vector3d coord_0 = sottopoligono.Cell0DCoordinates[id_0];
    Vector3d coord_1 = sottopoligono.Cell0DCoordinates[id_1];

    // trovo i vettori che congiungono gli estremi del primo lato al baricentro
    Vector3d v1 = coord_0 - bar_vec;
    Vector3d v2 = coord_1 - bar_vec;
    Vector3d v1xv2 = vec_product(v1, v2);
    double prod_scal = v1xv2.dot(vett_normale_frattura);

    // se il prodotto scalare è negativo => devo prendere l'altro senso
    if(prod_scal < 0){
        reverse(id_estremi_lato.begin(), id_estremi_lato.end());
        reverse(id_lati.begin(), id_lati.end());
    }

    // aggiorno Cell2D
    sottopoligono.Cell2DId.push_back(num_sottopoligono); // aggiorno l'id del sottopoligono
    sottopoligono.NumberVertices.push_back(n); // aggiorno num vertici
    sottopoligono.NumberEdges.push_back(n); // aggiorno num lati

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





}//del namespace
