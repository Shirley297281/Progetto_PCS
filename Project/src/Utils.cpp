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


    vector<double> distance_vector(num_vertices-1, 0.0);
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

    fracture.NumFractures = num_fractures;

    fracture.IdFractures.reserve(num_fractures);
    fracture.numVertices.reserve(num_fractures);
    fracture.lenghtMaxEdges.reserve(num_fractures);
    fracture.baricentro.reserve(num_fractures);

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

            cout << "Analisi coppia di fratture: " << i << " e " << j << endl;


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
                cout<<"non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
                continue;
            }
    // TRACCE SI o NO?
    /// controllo se possano esistere tracce per la frattura i

            MatrixXd matrixVerticesI = fracture.CoordinatesVertice[i];
            unsigned int iterI = 0; //
            int numColonneI = matrixVerticesI.cols();
            vector<Vector3d> vecI; // contiene i punti di intersezione calcolati per la frattura i --> vertici della traccia
            vecI.reserve(numColonneI);//vado a salvarmi le intersezioni tra retta passante per due vertici e retta che individua l'intersezioni tra piani. In generale
            Vector3d Punto0 = {};
            VectorXd freePar;

            for (int z = 0; z < numColonneI; z++) {
                Vector3d V1 = matrixVerticesI.col(z);
                Vector3d V2;
                if (z == numColonneI - 1) { // Accoppio l'ultimo vertice con il primo
                    V2 = matrixVerticesI.col(0);
                } else {
                    V2 = matrixVerticesI.col(z + 1);
                }
             /// calcolo retta (GeometryComputational 1) dove giacciono i due vertici



                bool a;
                a = soluzione_sistema3x2(t,V1,V2,Point,Punto0);

                cout <<"punto intersezione tra retta int piani e retta tra due vertici: " << Punto0<<endl;

                if (!a) //se il sistema non ha soluzione cambio lato
                {
                    cout << "\t sono parallele!!!"<<endl;
                    continue;
                }

                double freeParP0;
                for (int i=0;i<3;i++) //presupponiamo che le fratture non siano degeneri
                {
                    if (abs(V2[i]-V1[i])>1e-14)
                    {
                        freeParP0 = (Punto0[i]-V1[i])/(V2[i]-V1[i]);
                        break;
                    }
                }
                cout<<"    "<<endl;
                cout <<"parametro libero di P0: " <<freeParP0<<endl;
                cout<<"    "<<endl;



                /// valuto Punto0, V1 e V2 nella retta passante per i due vertici
                /// (P=s*(V2-V1)+V1 ad esempio) e calcolo il valore che assume il parametro libero per ognuno

                double freeParV1 = 0;
                double freeParV2 = 1;
                /// controllo se il Punto0 può essere scritto come combinazione convessa dei due vertici ovvero freeParP0 appartiene a 0 o 1
                if (freeParP0 >= freeParV1 - 1e-15 && freeParP0 <= freeParV2 + 1e-15) {
                    iterI++;
                    vecI.push_back(Punto0);
                }
            }//fine for sui punti

            //cout << " il punto di intersezione tra il lato della frattura "<<i<<" e "<<j<<" è : "<< Punto0.transpose() <<endl;
            //cout<<"\t\tfinito il controllo tra "<<i<<" e "<<j<<endl;
            if (iterI != 2) // la frattura non può avere traccia con l'altra frattura, mi fermo
            {
                vecI.clear();
                continue;
            }



            /// controllo se possano esistere tracce per la frattura j (stesso procedimento cambio solo i nomi)
            MatrixXd matrixVerticesJ = fracture.CoordinatesVertice[j];
            unsigned int iterJ = 0;
            int numVerticiJ = matrixVerticesJ.cols();
            vector<Vector3d> vecJ; //vertici traccia per la frattura j
            vecJ.reserve(numVerticiJ);
            for (int z = 0; z < numVerticiJ; z++)
            {
                Vector3d V3 = matrixVerticesJ.col(z);
                Vector3d V4;
                if (z == numVerticiJ - 1) {
                    V4 = matrixVerticesJ.col(0);
                } else {
                    V4 = matrixVerticesJ.col(z + 1);
                }
                bool b;
                b = soluzione_sistema3x2(t,V3,V4,Point,Punto0);
                if (!b) //se il sistema non ha soluzione cambio lato
                {
                    continue;
                }

                double freeParP0;
                for (int i=0;i<3;i++) //presupponiamo che le fratture non siano degeneri
                {
                    if (abs(V4[i]-V3[i])>1e-15)
                    {
                        freeParP0 = (Punto0[i]-V3[i])/(V4[i]-V3[i]);
                        break;
                    }
                }
                cout<<"    "<<endl;
                cout <<"parametro libero di P0: " <<freeParP0<<endl;
                cout<<"    "<<endl;



                double freeParV3 = 0;
                double freeParV4 = 1;
                if (freeParP0 >= freeParV3 - 1e-15 && freeParP0 <= freeParV4 + 1e-15) {
                    iterJ++;
                    vecJ.push_back(Punto0);
                }
            }

            // Stampa di debug per iterI e iterJ
            cout << "Frattura " << i << " ha " << iterI << " punti di intersezione." << endl;
            cout << "Frattura " << j << " ha " << iterJ << " punti di intersezione." << endl;

            if (iterJ != 2 ) // la frattura non può avere traccia con l'altra frattura e passo ad altre due fratture
            {
                // non so se sia necessario, chat dice che uscendo di scope si distrugge automaticamente
                vecJ.clear();
                cout << "non c'è traccia tra "<< i << " e "<< j << endl;

                continue;
            }



            // Memorizzo i dati nella struttura Traces
            trace.numTraces++;
            trace.IdTraces.push_back(trace.numTraces);
            Matrix<double, 3, 2> estremiTraccia;
            estremiTraccia.col(0) = vecI[0];
            estremiTraccia.col(1) = vecI[1];
            trace.CoordinatesEstremiTraces.push_back(estremiTraccia);
            ///da modificare
            double lunghezzaTraccia = euclidean_distance(vecI[0], vecI[1]);
            trace.lengthTraces.push_back(lunghezzaTraccia);

            // Visualizza i dati della traccia
            cout << "Traccia " << trace.numTraces << ":" << endl;
            cout << " - Frattura 1: " << i << endl;
            cout << " - Frattura 2: " << j << endl;
            cout << " - Estremi traccia: (" << vecI[0].transpose() << "), (" << vecI[1].transpose() << ")" << endl;
            cout << " - Lunghezza traccia: " << lunghezzaTraccia << endl;


        // PASSANTI O NO? Metodo di Cri
           /// se siamo arrivati a questo punto è perchè posso identificare la traccia in passante e non passante
            //riasseganre i "punteggi"
            double freeParP1=(vecI[0][0]-Point[0])/(t[0]);
            double freeParP2=(vecI[1][0]-Point[0])/(t[0]);
            double freeParP3 = (vecJ[0][0] - Point[0]) / t[0];
            double freeParP4 = (vecJ[1][0] - Point[0]) / t[0];
            Vector2d Par1 = {freeParP1, freeParP2 };
            Vector2d Par2 = {freeParP3, freeParP4};
            sort(Par1.begin(), Par2.end());
            sort(Par1.begin(), Par2.end());  //finchè non ho bubble sort uso sort della libreria STL



            /// i casi visionati di seguito coprono le varie possibilità e sono in ordine rispetto a quanto mostrato su file latex
            // creo un dizionario per poter riottenere le informazioni sul punto sulla retta a partire dal parametro libero
            map<double, Vector3d> dizfreeParToVecI;
            dizfreeParToVecI[freeParP1] = vecI[0];//primo vertice della traccia per la frattura i
            dizfreeParToVecI[freeParP2] = vecI[1];//secondo vertice della traccia per la frattura j
            dizfreeParToVecI[freeParP3] = vecJ[0]; //primo vertice della traccia per la j
            dizfreeParToVecI[freeParP4] = vecJ[1]; //secondo vertice della traccia per la j
            // popolo il dizionario (non mi interessa se ho una chiave che si ripete (quindi questa verrebbe sovrascritta) tanto il dizionario viene usato negli ultimi elseif)

            /// in calcolo tracce (nel seguito non so se bisogna usare tolleranze o no)
            if (Par1[1] <= Par2[0] || Par1[0] <= Par2[1]) {
                continue;
            } else {
                Vector3d point1, point2;
                bool isPassante;

                if (Par1[0] == Par2[0] && Par2[1] == Par1[1]) {
                    point1 = vecI[0];
                    point2 = vecI[1];
                    isPassante = true;
                } else if (Par1[0] >= Par2[0] && Par1[1] <= Par2[1]) {
                    point1 = vecI[0];
                    point2 = vecI[1];
                    isPassante = true;
                } else if (Par2[0] >= Par1[0] && Par2[1] <= Par1[1]) {
                    point1 = vecJ[0];
                    point2 = vecJ[1];
                    isPassante = true;
                } else if (Par2[1] > Par1[0]) {
                    point1 = dizfreeParToVecI[Par2[1]];
                    point2 = dizfreeParToVecI[Par1[1]];
                    isPassante = false;
                } else if (Par1[1] > Par2[0]) {
                    point1 = dizfreeParToVecI[Par1[1]];
                    point2 = dizfreeParToVecI[Par2[0]];
                    isPassante = false;
                } else {
                    cout << "Non ho considerato questo caso: fratture " << i << " e " << j << endl;
                    continue;
                }

                trace.numTraces++;
                trace.IdTraces.push_back(trace.numTraces);

                Matrix<double, 3, 2> estremi;
                estremi.col(0) = point1;
                estremi.col(1) = point2;
                trace.CoordinatesEstremiTraces.push_back(estremi);

                double length = (point1 - point2).norm();
                trace.lengthTraces.push_back(length);

                trace.vectorTips.push_back(isPassante);

                cout << "Numero di tracce totali: " << trace.numTraces << endl;
                for (unsigned int i = 0; i < trace.numTraces; ++i) {

                    cout << "Traccia " << trace.IdTraces[i] << ":" << endl;
                    cout << " - Frattura 1: " << i << endl;
                    cout << " - Frattura 2: " << j << endl;
                    cout << "Coordinate Estremi:" << endl;
                    cout << "  Estremo 1: ("
                         << trace.CoordinatesEstremiTraces[i](0, 0) << ", "
                         << trace.CoordinatesEstremiTraces[i](1, 0) << ", "
                         << trace.CoordinatesEstremiTraces[i](2, 0) << ")" << endl;
                    cout << "  Estremo 2: ("
                         << trace.CoordinatesEstremiTraces[i](0, 1) << ", "
                         << trace.CoordinatesEstremiTraces[i](1, 1) << ", "
                         << trace.CoordinatesEstremiTraces[i](2, 1) << ")" << endl;

                    cout << "Lunghezza della traccia: " << trace.lengthTraces[i] << endl;
                    cout << "Traccia passante: " << trace.vectorTips[i]  << endl;


                }



            }

            /*bool s;
            s=Tips_Shy(vecI, vecJ);


            cout << " - Passante: " << (s ? "Sì" : "No") << endl;
            cout << endl;*/


        }//end for j



    }//end for i
    cout << "\nescluse in principio "<< escluse<< " possibili intersezioni! SBAM."<<endl;

//magari ancora qua dentro facciamo il primo file output


}//Calcolo tracce

bool Tips_Shy(const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ) {
    // Verifica se almeno un estremo di vecI tocca un lato di vecJ
    for (const auto& pointI : vecI) {
        for (size_t k = 0; k < vecJ.size(); ++k) {
            size_t next = (k + 1) % vecJ.size();
            double dist1 = euclidean_distance(pointI, vecJ[k]);
            double dist2 = euclidean_distance(pointI, vecJ[next]);
            double segmentLength = euclidean_distance(vecJ[k], vecJ[next]);
            if (abs(dist1 + dist2 - segmentLength) < 1e-6) {
                // L'estremo di vecI tocca il lato di vecJ
                // Continua con la verifica dell'altro estremo di vecI
                // e verifica se anche tocca un lato di vecJ
                for (const auto& pointJ : vecJ) {
                    dist1 = euclidean_distance(pointJ, vecJ[k]);
                    dist2 = euclidean_distance(pointJ, vecJ[next]);
                    segmentLength = euclidean_distance(vecJ[k], vecJ[next]);
                    if (abs(dist1 + dist2 - segmentLength) < 1e-6) {
                        // Entrambi gli estremi di vecI e vecJ toccano i lati l'uno dell'altro
                        return true; // Traccia passante
                    }
                }
            }
        }
    }
    // Nessun estremo di vecI tocca un lato di vecJ
    return false; // Traccia non passante
}

}
