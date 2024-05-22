#include "Utils.hpp"
#include "FracturesTraces.hpp"
#include <vector>
#include "Eigen/Eigen"
#include <cmath> // per sqrt

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono> //for counting time
#include <stdlib.h>

double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices)
{


    vector<double> distance_vector(num_vertices, 0.0);
    // Itero sulle colonne
    for(unsigned int i = 0; i < num_vertices - 1; ++i){
        Vector3d point_a = m.col(i);
        Vector3d point_b = m.col(i + 1);
        distance_vector[i] = euclidean_distance(point_a, point_b);
    }

    // Calcolo la distanza tra l'ultimo punto e il primo
    Vector3d last_point = m.col(num_vertices - 1);
    Vector3d first_point = m.col(0);
    distance_vector[num_vertices - 1] = euclidean_distance(last_point, first_point);

    // Trovo la distanza massima
    auto it_max_distance = max_element(distance_vector.begin(), distance_vector.end());
    double max_distance = *it_max_distance;

    return max_distance;

}


//baricentro
array <double,3> barycenter(const MatrixXd& m, unsigned int num_vertices)
{

    array<double,3> barycenter_coord ={};

    for(unsigned int i = 0; i < 3; ++i){
        for(unsigned int j = 0; j < num_vertices ; ++j){
            barycenter_coord[i] += m(i,j);
        }
        barycenter_coord[i] /= num_vertices;
    }
    return barycenter_coord;
}

//funzione sfere
bool check_sphere(const array<double,3> bar1, const array<double,3> bar2, const double l1, const double l2)
{
    //controlliamo la distanza tra i due baricentri
    double distance_bar = 0.0;
    distance_bar = sqrt( (bar1[0] - bar2[0])*(bar1[0] - bar2[0]) +
                         (bar1[1] - bar2[1])*(bar1[1] - bar2[1]) +
                         (bar1[2] - bar2[2])*(bar1[2] - bar2[2]) );

    double max_distance = 0.0;
    max_distance = (l1 + l2) ;///h

    if (distance_bar > max_distance){
        return false; // le fratture non si intersecano
    }

    return true;
}

//normale al piano
Vector3d normal_vector(const MatrixXd& m)
{

    Vector3d v1 = {}; // vettore1
    Vector3d v2 = {}; // vettore2
    Vector3d v = {};

    v1[0] = m(0,1) - m(0,0);
    v1[1] = m(1,1) - m(1,0);
    v1[2] = m(2,1) - m(2,0);

    v2[0] = m(0,1) - m(0,2);
    v2[1] = m(1,1) - m(1,2);
    v2[2] = m(2,1) - m(2,2);

    // prodotto vettoriale
    v = vec_product(v1, v2);

    return v;

}


namespace GeometryLibrary
{


bool ImportFR(const string &filename,
                   Fractures& fracture)
{

    ifstream file;
    file.open(filename);

    if (file.fail())
    {
        cout << "Errore durante l'apertura del file." << endl;
        return false;
    }

    cout << "File "<<filename<<" aperto correttamente." << endl;

    string line;
    getline(file,line); //ignora prima riga
    getline(file,line);
    int num_fractures;
    istringstream converter(line);
    string token;

    getline(converter, token); // Utilizza converter invece di file per ottenere il token dalla stringa
    istringstream(token) >> num_fractures;

    cout << "num fractrues = "<< num_fractures<<endl;

    fracture.NumFractures = num_fractures;

    fracture.IdFractures.reserve(num_fractures);
    fracture.numVertices.reserve(num_fractures);
    fracture.lenghtMaxEdges.reserve(num_fractures);
    fracture.baricentro.reserve(num_fractures);
    //fracture.CoordinatesVertice.resize(num_fractures);

    //introduco un count per far arrestare la lettura
    int count = 0;

    while (count < num_fractures) {

        getline(file,line);//ignora la riga "# FractureId; NumVertices"
        getline(file,line);
        //memorizza i due valori id e num_vertices
        int id_fracture;
        istringstream ss(line);
        string token;
        getline(ss, token, ';'); // Utilizza ss invece di file per ottenere il token dalla stringa
        istringstream(token) >> id_fracture; // Usa istringstream per convertire il token in un intero
        fracture.IdFractures.push_back(id_fracture);

        ss.ignore(1); // Ignora lo spazio dopo il separatore
        token.clear();

        int num_vertices;
        getline(ss, token);
        istringstream(token) >> num_vertices;
        fracture.numVertices.push_back(num_vertices);

        //matrice dei vertici
        MatrixXd vertices(3, num_vertices);


        //ignora la riga "# Vertices"
        getline(file, line);
        //for sulle tre linee che segono per le 3d
        for(int i = 0; i<3; i++){
            //ignora la riga "# Vertices"
            getline(file,line);
            double coord;
            istringstream cc(line);
            string coordinate;

            for (int j=0; j<num_vertices;j++){

                getline(cc, coordinate, ';'); //legge il valore prima del ;
                istringstream(coordinate) >> coord; // Usa istringstream per convertire il coordinate in un double

                vertices(i,j)=coord;

                cc.ignore(1); // Ignora lo spazio dopo il ;
                coordinate.clear();
            }

        }

        fracture.CoordinatesVertice.push_back(vertices);

        double length;
        length = max_euclidean_distance(vertices, num_vertices);
        fracture.lenghtMaxEdges.push_back(length);

        array <double,3> bary;
        bary = barycenter(vertices, num_vertices);
        fracture.baricentro.push_back(bary);

        Vector3d v_normal;
        v_normal = normal_vector(vertices);
        fracture.vettoreNormalePiano.push_back(v_normal);


        count ++;//frattura letta √

    }//eof

    file.close();


    //controllo di quel che si ha memorizzato

    // /!\ FORSE CI CONVIENE CONTROLLARE CON L'OUTPUT QUALCOSA
    //     DI SIGNIFICATIVO SULLO STILE DEI MERKER DELL'ESERCITAZIONE 5
    cout << "Dati delle fratture memorizzati correttamente." << endl;
    /*for (size_t i = 0; i < fracture.IdFractures.size(); ++i) {
        cout << "--- Frattura " << fracture.IdFractures[i] << " ---" << endl;
        cout << "Numero vertici: " << fracture.numVertices[i] << endl;
        cout << "Coordinate vertici:" << endl;
        cout << fracture.CoordinatesVertice[i] << endl;
        cout << endl;
    }

    // display the vector elements using a for loop
    for (int i = 0; i < num_fractures; i++) {
        cout << "\nmax length of fracture[" << i << "] = " <<scientific<<setprecision(6)<< fracture.lenghtMaxEdges[i] << endl;
        cout << "baricentro: (" << fracture.baricentro[i][0] << ", " << fracture.baricentro[i][1] << ", " << fracture.baricentro[i][2] << ")" << std::endl;
        cout << "vettore normale: (" << fracture.vettoreNormalePiano[i][0] << ", " << fracture.vettoreNormalePiano[i][1] << ", " << fracture.vettoreNormalePiano[i][2] << ")" << endl;
    }*/



    return true;

}



void CalcoloTracce(Fractures& fracture, Traces& trace)
{
    int escluse = 0;

    /*
    // sovrastima ragionevole basata sul numero massimo di combinazioni di fratture = binomiale
    const unsigned int estimatedNumTraces = (fracture.NumFractures * (fracture.NumFractures - 1)) * 0.5;
    // Stima del numero di tracce: uso la sovrastima perchè riservare troppo spazio in anticipo può sprecare memoria,
    //  ma le operazioni di inserimento saranno efficienti. se SOTTOSTIMASSI dovrei gestire riallocazioni multiple (tempo++)

    // Riserva di spazio per migliorare l'efficienza delle allocazioni di memoria
    trace.IdTraces.reserve(estimatedNumTraces);
    trace.CoordinatesEstremiTraces.reserve(estimatedNumTraces);
    trace.lengthTraces.reserve(estimatedNumTraces);
    trace.vectorTips.reserve(estimatedNumTraces);
    */

    for (unsigned int i = 0; i< fracture.NumFractures - 1 ; i++ ){
        for (unsigned int j=i+1; j < fracture.NumFractures; j++ ){

            cout << "\nAnalisi coppia di fratture: " << i << " e " << j << endl;


            // richiamo funzione sfera
            if ( !check_sphere( fracture.baricentro[i], fracture.baricentro[j], fracture.lenghtMaxEdges[i], fracture.lenghtMaxEdges[j]) )
            {
                escluse++;
                continue;
            }

            // t è la tangente che inidividua la direzione della retta di intersezione
            Vector3d t = {};
            t = vec_product(fracture.vettoreNormalePiano[i], fracture.vettoreNormalePiano[j]);

            //il nostro P0 è il baricentro (COMPUTATIONAL GEOMETRY 2, PROBLEMA 4)
            //Point è il punto che appartiene alla retta(intersezione tra piani) che ha parametro libero x e direzione t (tangente)
            Vector3d Point = {};
            bool ris = system_solution(fracture.vettoreNormalePiano[i], fracture.vettoreNormalePiano[j],
                                    fracture.baricentro[i], fracture.baricentro[j],
                                    t, Point);
            //tutto questo ci serve per trovare la retta di intersezione r(x) = x*t_ + Point;



            ///X SHY
    // in calcolo tracce ma allinizio bisogna:
    // 2) inizializzare numTraces, CoordinatesEstremiTraces, lenghtTraces, vectorTips che vengono usate nel codice successivo
    // il problema è che mi sembra complicato stimare il numero di tracce per fare un primo resize quindi potremmo pensare di partire da un numero
    // non altissimo e di aggiustare manmano con l'operazione di raddoppio spiegata da Berrone nella prima lezione (se non ricordo male c'è proprio una
    // funzione raddoppio nelle slide) giuro che qua la pianto con i commenti, ecco inizia il codice da inserire a riga 261
            if (!ris)
            {
                cout<<" || ! ||  non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
                continue;
            }


    // TRACCE SI o NO?
    /// controllo se possano esistere tracce per la frattura i
    ///
            vector<Vector3d> vecI = {};  // contiene i punti di intersezione calcolati per la frattura i --> vertici della traccia
            unsigned int iterI = Calcolo_par(t, Point, i, vecI, fracture);

            // cout << "Frattura " << i << " ha " << iterI << " punti di intersezione." << endl;

            if (iterI != 2) // la frattura non può avere traccia con l'altra frattura, mi fermo
            {
                vecI.clear();
                cout<<"non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
                continue;
            }


            vector<Vector3d> vecJ = {};  // contiene i punti di intersezione calcolati per la frattura i --> vertici della traccia
            unsigned int iterJ = Calcolo_par(t, Point, j, vecJ, fracture);

            // cout << "Frattura " << i << " ha " << iterI << " punti di intersezione." << endl;

            if (iterJ != 2) // la frattura non può avere traccia con l'altra frattura, mi fermo
            {
                vecJ.clear();
                cout<<"non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
                continue;

            }


            // Memorizzo i dati nella struttura Traces
            trace.numTraces++;
            trace.IdTraces.push_back((trace.numTraces)-1);
            //è sbagliato qui, da spostare dopo
            /*Matrix<double, 3, 2> estremiTraccia;
            estremiTraccia.col(0) = vecI[0];
            estremiTraccia.col(1) = vecI[1];
            trace.CoordinatesEstremiTraces.push_back(estremiTraccia);*/
            ///da modificare
            /* double lunghezzaTraccia = euclidean_distance(vecI[0], vecI[1]);
            trace.lengthTraces.push_back(lunghezzaTraccia);*/

            // Visualizza i dati della traccia
            cout << "\nTraccia " << trace.numTraces << ":" << endl;
            cout << " - Frattura 1: " << i << endl;
            cout << " - Frattura 2: " << j << endl;
            //da spostare dopo
            //cout << " - Estremi traccia: (" << vecI[0].transpose() << "), (" << vecI[1].transpose() << ")" << endl;
            //cout << " - Lunghezza traccia: " << lunghezzaTraccia << endl;


            int result = Controllo_tracce2(fracture,trace,vecI,vecJ,Point,t,i,j);
            if (result == 1)
            {
                cout<<"non ho intersezione perchè sono nel caso 1 e 2"<<endl;
                continue;
            }
            else if (result == 3)
            {
                cout<<"non ho considerato questa casistica: "<< i << " e "<< j<< endl;
                continue;
            }




        }//end for j



    }//end for i


    cout << "\nescluse in principio "<< escluse<< " possibili intersezioni! SBAM."<<endl;

//magari ancora qua dentro facciamo il primo file output


}//Calcolo tracce


unsigned int Calcolo_par(Vector3d& t, Vector3d& Point, int i, vector<Vector3d>& vec, Fractures& fracture){

    MatrixXd matrixVertices = fracture.CoordinatesVertice[i];
    unsigned int iter = 0; //valore da resituire = numero di intersezioni trovate
    int numColonne = matrixVertices.cols();
    vec.reserve(numColonne);//vado a salvarmi le intersezioni tra retta passante per due vertici e retta che individua l'intersezioni tra piani. In generale
    Vector3d Punto0 = {}; // è il punto di intersezione tra la retta tangente e la retta tra due vertici V1 e V2

    for (int z = 0; z < numColonne; z++) {
        Vector3d V1 = matrixVertices.col(z);
        Vector3d V2;
        if (z == numColonne - 1) { // Accoppio l'ultimo vertice con il primo
            V2 = matrixVertices.col(0);
        } else {
            V2 = matrixVertices.col(z + 1);
        }

        bool a;
        a = soluzione_sistema3x2(t,V1,V2,Point,Punto0);

        if (!a) //se il sistema non ha soluzione cambio lato
        {
            //cout << "\t sono parallele!!!"<<endl;
            continue;
        }else{
            //cout <<"punto intersezione tra retta int piani e retta tra due vertici: " << Punto0<<endl;
        }

        double freeParP0 = 0.0; // parametro libero che corrisponde al Punto0
        for (int i=0;i<3;i++) //presupponiamo che le fratture non siano degeneri
        {
            if (abs(V2[i]-V1[i])>1e-14)
            {
                freeParP0 = (Punto0[i]-V1[i])/(V2[i]-V1[i]);
                break;
            }
        }

        /// valuto Punto0, V1 e V2 nella retta passante per i due vertici
        /// (P=s*(V2-V1)+V1 ad esempio) e calcolo il valore che assume il parametro libero per ognuno

        double freeParV1 = 0;
        double freeParV2 = 1;
        /// controllo se il Punto0 può essere scritto come combinazione convessa dei due vertici ovvero freeParP0 appartiene a 0 o 1
        if (freeParP0 >= freeParV1 - 1e-15 && freeParP0 <= freeParV2 + 1e-15) {

            //controllo se c'è già il punto
            if (find(vec.begin(), vec.end(), Punto0) == vec.end()) {
                // Se il valore non è presente, aggiungilo al vettore
                vec.push_back(Punto0);
                iter++;
            }

        }


    }//fine for sui vertici della frattura

    //cout << " il punto di intersezione tra il lato della frattura "<<i<<" e "<<j<<" è : "<< Punto0.transpose() <<endl;
    //cout<<"\t\tfinito il controllo tra "<<i<<" e "<<j<<endl;

    return iter;


}

int Controllo_tracce2(Fractures fracture, Traces trace, const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ,
                      const Vector3d Point, Vector3d t, const unsigned int i, const unsigned int j)
{   // i e j solo per il cout
    double freeParP1;
    double freeParP2;
    double freeParP3;
    double freeParP4;
    for (unsigned int s=0; s<3; s++)
    {
        if (abs(t[s])>1e-15){ //mi baso sul fatto che esista almeno una coordinata di t diversa da zero altrimenti non saremmo arrivati qua
            freeParP1=(vecI[0][s]-Point[s])/(t[s]);
            freeParP2=(vecI[1][s]-Point[s])/(t[s]);
            freeParP3=(vecJ[0][s]-Point[s])/(t[s]);
            freeParP4=(vecJ[1][s]-Point[s])/(t[s]);
        }
    }
    // creo un dizionario per poter riottenere le informazioni sul punto sulla retta a partire dal parametro libero
    map<double, Vector3d> dizfreeParToVec;
    // popolo il dizionario (non mi interessa se ho una chiave che si ripete (quindi questa verrebbe sovrascritta) tanto il dizionario viene usato negli ultimi elseif)
    dizfreeParToVec[freeParP1] = vecI[0];
    dizfreeParToVec[freeParP2] = vecI[1];
    dizfreeParToVec[freeParP3] = vecJ[0];
    dizfreeParToVec[freeParP4] = vecJ[1];
    /// in calcolo tracce (nel seguito non so se bisogna usare tolleranze o no.
    // ordino i parametri liberi di ciascun vec
    Vector2d Par1 = {freeParP1, freeParP2};
    Vector2d Par2 = {freeParP3, freeParP4};
    sort(Par1.begin(), Par2.end());
    sort(Par1.begin(), Par2.end());  //finchè non ho bubble sort uso sort della libreria STL



    ///SI INTROMETTE QUELLA FICCANASO DI SHY
    ///
    ///

    Tips_Shy(fracture, trace, vecI, vecJ, i, j, dizfreeParToVec, freeParP1, freeParP2, freeParP3, freeParP4);

    ///
    ///
    /// fine





    /// i casi visionati di seguito coprono le varie possibilità e sono in ordine rispetto a quanto mostrato su file inviato su whatsapp scritto su tablet
    if (Par1[1]<=Par2[0] || Par1[0]>=Par2[1])     // caso 1 e 2
    {
        // devo scartare le tracce
        return 1;
    }
    else if (Par1[0]==Par2[0] && Par2[1] == Par1[1]) // per entrambe le fratture la traccia è passante, caso 3
    {
        cout<<"traccia passante in "<<i<<" e "<<j<<endl;
        //euclidean distance tra vecI[0] e vecI[1] per calcolare lunghezza traccia e inserirla in lenghtTraces
        //inserire vecI[0] e vecI[1] in CoordinatesEstremiTraces
        //numTraces ++
        //inserire Id traccia (che coindice con il numtraces) in TraceIdsPassxFracture sia per i che per j (usare metodo insert visto in esercitazione 5
        //in modo da gestire l'inserimento se la chiave = IdFrattura è già inserita)
    }
    else if (Par1[0]>= Par2[0] && Par1[1] <= Par2[1]) //solo per la frattura i la traccia è passante, caso 4 e 6
    {
        cout<<"traccia passante solo per "<<i<<" ma non per "<<j<<endl;
        //euclidean distance tra vecI[0] e vecI[1]  per calcolare lunghezza traccia e inserirla in lenghtTraces
        //inserire vecI[0] e vecI[1] in CoordinatesEstremiTraces
        //numTraces ++
        //inserire l'Id traccia (che coindice con il numtraces) in TraceIdsPassxFracture con chiave i
        //inserire l'Id traccia (che coindice con il numtraces) in TraceIdsNoPassxFracture con chiave j
    }
    else if (Par2[0]>= Par1[0] && Par2[1] <= Par1[1]) //solo per la fratture j la traccia è passante, caso 4 e 6
    {
        cout<<"traccia passante solo per "<<j<<" ma non per "<<i<<endl;
        //euclidean distance tra vecJ[0] e vecJ[1]  per calcolare lunghezza traccia e inserirla in lenghtTraces
        //inserire vecJ[0] e vecJ[1] in CoordinatesEstremiTraces
        //numTraces ++
        //inserire l'Id traccia (che coindice con il numtraces) in TraceIdsPassxFracture con chiave j
        //inserire l'Id traccia (che coindice con il numtraces) in TraceIdsNoPassxFracture con chiave i
    }
    else if (Par2[1] > Par1[0])  //non passante per entrambe, caso 5
    {
        cout<<"traccia non passante sia per "<<j<<" che per "<<i<<endl;
        //euclidean distance tra dizfreeParToVecI[Par2[1]] e dizfreeParToVecI[Par1[1]]  per calcolare lunghezza traccia e inserirla in lenghtTraces
        //inserire dizfreeParToVecI[Par2[1]]  e dizfreeParToVecI[Par1[1]] in CoordinatesEstremiTraces
        //numTraces ++
        //inserire Id traccia (che coindice con il numtraces) in TraceIdsNoPassxFracture sia per i che per j
    }
    else if ( Par1[1] > Par2[0])   //non passante per entrambe, caso 5
    {
        //euclidean distance tra dizfreeParToVecI[Par1[1]] e dizfreeParToVecI[Par2[0]]  per calcolare lunghezza traccia e inserirla in lenghtTraces
        //inserire dizfreeParToVecI[Par1[1]]  e dizfreeParToVecI[Par2[0]] in CoordinatesEstremiTraces
        //numTraces ++
        //inserire Id traccia (che coindice con il numtraces) in TraceIdsNoPassxFracture sia per i che per j
    }
    else
    {
        return 3; //controllo se ci sono casistiche non considerate
    }
    return 0;
}


bool Tips_Shy(Fractures fracture, Traces trace,const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ, const int i, const int j ,map<double, Vector3d>& dizfreeParToVec,
              double freeParP1,
              double freeParP2,
              double freeParP3,
              double freeParP4) {

    //creo un vettore di 4 elementi da ordinare
    vector<double> intersezioni= {freeParP1, freeParP2, freeParP3, freeParP4};

    BubbleSort(intersezioni);

    //creo la matrcie 3 righe 2 colonne da inserire nel vettore delle matrici degli estremi.
    Matrix<double, 3, 2> Estremi;

    // Controllo che le chiavi esistano nel dizionario prima di accedervi
    if (dizfreeParToVec.find(intersezioni[1]) != dizfreeParToVec.end() &&
        dizfreeParToVec.find(intersezioni[2]) != dizfreeParToVec.end()) {

        // Riempio la matrice con i punti dai valori del dizionario
        Estremi.col(0) = dizfreeParToVec[intersezioni[1]];
        Estremi.col(1) = dizfreeParToVec[intersezioni[2]];
    } else {
        cerr << "Chiavi non trovate nel dizionario!" << endl;
        return false; // O gestisci l'errore in modo appropriato
    }

    trace.CoordinatesEstremiTraces.push_back(Estremi);

    cout << " - Estremi traccia: (" << Estremi.col(0).transpose() << "), (" << Estremi.col(1).transpose() << ")" << endl;

    trace.lengthTraces.push_back(euclidean_distance(Estremi.col(0), Estremi.col(1)));


    int touchCount = 0;
    bool passa;

    MatrixXd vertice = fracture.CoordinatesVertice[i];
        for (int k = 0; k < vertice.cols(); ++k) {
            Vector3d v1 = vertice.col(k);
            Vector3d v2 = vertice.col((k + 1) % vertice.cols());;
            double edgeLength = euclidean_distance(v1, v2);

            for (int e = 0; e < 2; ++e) {
                Vector3d extremity = Estremi.col(e);
                if (abs(euclidean_distance(extremity, v1) + euclidean_distance(extremity, v2) - edgeLength) < 1e-14) {
                    touchCount++;
                    break;
                }
            }
        }

    if (touchCount == 2){
        //la traccia è passante
        passa = true;

    }else{
        passa = false;

    }

    cout << " - Passante per la frattura " <<i<< " : " << (passa ? "Sì" : "No") << endl;


    touchCount = 0;
    vertice = fracture.CoordinatesVertice[j];
    for (int k = 0; k < vertice.cols(); ++k) {
        Vector3d v1 = vertice.col(k);
        Vector3d v2 = vertice.col((k + 1)% vertice.cols());
        double edgeLength = euclidean_distance(v1, v2);

        for (int e = 0; e < 2; ++e) {
            Vector3d extremity = Estremi.col(e);
            if (abs(euclidean_distance(extremity, v1) + euclidean_distance(extremity, v2) - edgeLength) < 1e-14) {
                touchCount++;
                break;
            }
        }
    }

    if (touchCount == 2){
        //la traccia è passante
        passa = true;

    }else{
        passa = false;

    }

    cout << " - Passante per la frattura " <<j<< " : " << (passa ? "Sì" : "No") << endl;
    cout << endl;


}

}
