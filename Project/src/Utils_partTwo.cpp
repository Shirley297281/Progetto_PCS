#include "Utils.hpp"
#include "FracturesTraces.hpp"
#include "inline.hpp"
#include <vector>
#include "Eigen/Eigen"
#include <iostream>
#include <vector>
#include <list>
#include <set>


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

void MemorizzaVertici_Cell0Ds(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z){

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


void Creazioni_Sequenze(const Fractures& fracture, const Traces& trace, Polygons& sottoPoligono, unsigned int z)
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
    unsigned int n = estremi.size();
    vector<Vector2i> id_estremi_lato; // lati identificati dagli id degli estremi
    id_estremi_lato.reserve(n); // n vertici => avrò n lati
    MatrixXd vertices(3, n); // per baricentro
    Vector3d vett_normale_frattura = fracture.vettoreNormalePiano[num_fracture];
    vector<unsigned int> id_lati;
    id_lati.reserve(n);


    // con i primi due for seleziono due vertici e verifico che possano essere consecutivi (faccio il controllo con il terzo for)
    unsigned int num_iterazioni = 0;
    unsigned int i = 0;
    //for(unsigned int i = 0; i < n; i++ ){
    unsigned int id_i = estremi[i];
    unsigned int id_i_primo = id_i;
    unsigned int id_j = 0;
    Vector3d coord_i = sottopoligono.Cell0DCoordinates[id_i]; // seleziono il primo estremo

    for (unsigned int j = 0; j < n+1; j++)
    {
        bool lato_valido = true;

        if (j==i){continue;} // MANCA UN LATO
        if(num_iterazioni == n-1){
            id_i = id_estremi_lato[n-2][1];
            id_j = id_i_primo;
        }
        // aggiorna!!!!!
        else{
            id_j = estremi[j];
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
        if(lato_valido && num_iterazioni < 4){
            Vector2i l(id_i, id_j);
            id_estremi_lato.push_back(l);

            cout << "id_estremi_lato di " << i << ": "<<id_estremi_lato[i].transpose()<<endl;

            // verifico se il lato è già presente in Cell1D
            auto it = find(sottopoligono.Cell1DVertices.begin(), sottopoligono.Cell1DVertices.end(), l);
            Vector2i l_inverso(id_j, id_i);
            auto it1 = find(sottopoligono.Cell1DVertices.begin(), sottopoligono.Cell1DVertices.end(), l_inverso);

            if(it != sottopoligono.Cell1DVertices.end() || it1 != sottopoligono.Cell1DVertices.end()){ // l'ho già trovato
                num_iterazioni++;
                continue;
            }
            else if(it == sottopoligono.Cell1DVertices.end() || it1 == sottopoligono.Cell1DVertices.end()){
                sottopoligono.Cell1DVertices.push_back(l);
                unsigned int id;

                if(sottopoligono.Cell1DId.empty())
                {
                    id =0;

                }
                else if(!sottopoligono.Cell1DId.empty()){
                    id = sottopoligono.Cell1DId.back() + 1;
                }
                sottopoligono.Cell1DId.push_back(id);
                id_lati.push_back(id); // inizio a creare il vettore da inserire in Cell2DEdges (se li sto ordinando in senso orario piuttosto li inverto dopo)
                num_iterazioni += 1;
            }



        }
        i = j; // mi permette di trovare i lati in ordine
        id_i = estremi[i];//mi serve per andare avanti con i lati, altrimenti fa sempre riferimento al primo

        //break;

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
    set<int> id_estremi;
    for(auto it = id_estremi_lato.begin(); it!= id_estremi_lato.end(); ++it){
        Vector2i vec = *it;
        id_estremi.insert(vec[0]);
        id_estremi.insert(vec[1]);
    }

    vector<unsigned int> id_lati_vec(id_estremi.begin(), id_estremi.end());
    sottopoligono.Cell2DVertices[num_sottopoligono] = id_lati_vec;

    // Cell2DEdges: ogni sottopoligono identificato da vettore contenente l'id dei lati

    sottopoligono.Cell2DEdges[num_sottopoligono] = id_lati; // NON LO RIEMPIE




}





}//del namespace
