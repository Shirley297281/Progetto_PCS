#include "namespace.hpp" //contiene gli header di tutte le funzione definite come GeometryLibrary (anche Tol)
#include "Utils.hpp"
#include "inline.hpp"


#include <vector>
#include "Eigen/Eigen"
#include <cmath> // per sqrt

#include <iostream>
#include <vector>
#include <chrono> //for counting time
#include <stdlib.h>

using namespace std;
using namespace Eigen;


namespace GeometryLibrary{


void CalcoloTracce(Fractures& fracture, Traces& trace)
{
    int escluse = 0;


    // sovrastima ragionevole basata sul numero massimo di combinazioni di fratture = binomiale
    const unsigned int estimatedNumTraces = (fracture.NumFractures * (fracture.NumFractures - 1)) * 0.5;
    // Stima del numero di tracce: uso la sovrastima perchè riservare troppo spazio in anticipo può sprecare memoria,
    //  ma le operazioni di inserimento saranno efficienti. se SOTTOSTIMASSI dovrei gestire riallocazioni multiple (tempo++)

    // Riserva di spazio per migliorare l'efficienza delle allocazioni di memoria
    trace.IdTraces.reserve(estimatedNumTraces);
    trace.CoordinatesEstremiTraces.reserve(estimatedNumTraces);
    trace.lengthTraces.reserve(estimatedNumTraces);
    //aggiunta da shy come prova
    trace.TraceIdsPassxFracture.resize(fracture.NumFractures);
    trace.TraceIdsNoPassxFracture.resize(fracture.NumFractures);


    for (unsigned int i = 0; i< fracture.NumFractures - 1 ; i++ ){
        //aggiunta da shy come prova
        trace.TraceIdsPassxFracture[i].reserve(estimatedNumTraces);
        trace.TraceIdsNoPassxFracture[i].reserve(estimatedNumTraces);

        for (unsigned int j=i+1; j < fracture.NumFractures; j++ ){

            ///cout << "\nAnalisi coppia di fratture: " << i << " e " << j << endl;


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
                ///cout<<" || ! ||  non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
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
                ///cout<<"non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
                continue;
            }

            /// controllo se possano esistere tracce per la frattura j
            ///
            vector<Vector3d> vecJ = {};  // contiene i punti di intersezione calcolati per la frattura i --> vertici della traccia
            unsigned int iterJ = Calcolo_par(t, Point, j, vecJ, fracture);

            // cout << "Frattura " << i << " ha " << iterI << " punti di intersezione." << endl;

            if (iterJ != 2) // la frattura non può avere traccia con l'altra frattura, mi fermo
            {
                vecJ.clear();
                ///cout<<"non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
                continue;

            }

            //Riservo spazio stimato

            const unsigned int estimatedSize = (fracture.NumFractures * (fracture.NumFractures - 1)) / 2;

            // Reserving space for unordered_maps if applicable
            trace.TraceIdsPassxFracture.reserve(estimatedSize);
            trace.TraceIdsNoPassxFracture.reserve(estimatedSize);


            int result = Controllo_tracce2(fracture, trace,vecI,vecJ,Point,t,i,j);
            if (result == 1)
            {
                ///cout<<"non ho intersezione perchè sono nel caso 1 e 2"<<endl;
                continue;
            }
            else if (result == 2)
            {
                ///cout<<"non ho considerato questa casistica: "<< i << " e "<< j<< endl;
                continue;
            }




        }//end for j



    } //end for i

    /*
    cout << "TraceIdsPassxFracture contenuto:" << endl;
    for (unsigned int i = 0; i < trace.TraceIdsPassxFracture.size(); ++i) {
        cout << "Frattura " << i << ":";
        for (unsigned int j = 0; j < trace.TraceIdsPassxFracture[i].size(); ++j) {
            cout << " " << trace.TraceIdsPassxFracture[i][j];
        }
        cout << endl;
    }

    cout << "TraceIdsNoPassxFracture contenuto:" << endl;
    for (unsigned int i = 0; i < trace.TraceIdsNoPassxFracture.size(); ++i) {
        cout << "Frattura " << i << ":";
        for (unsigned int j = 0; j < trace.TraceIdsNoPassxFracture[i].size(); ++j) {
            cout << " " << trace.TraceIdsNoPassxFracture[i][j];
        }
        cout << endl;
    }


    cout << "\nescluse in principio "<< escluse<< " possibili intersezioni! SBAM."<<endl;*/

}//Calcolo tracce




int Controllo_tracce2(Fractures& fracture, Traces& trace, const vector<Vector3d>& vecI, const vector<Vector3d>& vecJ,
                      const Vector3d& Point, Vector3d& t, unsigned int i, unsigned int j)
{   // i e j solo per il cout
    double freeParP1 = 0.0;
    double freeParP2 = 0.0;
    double freeParP3 = 0.0;
    double freeParP4 = 0.0;

    vector<array<double,2>> idpar = {}; //è un vector di array da 2 che contengono l'id della frattura e il parametro del punto
    idpar.reserve(4);

    for (unsigned int s=0; s<3; s++)
    {
        if (abs(t[s])>tolDefault){ //mi baso sul fatto che esista almeno una coordinata di t diversa da zero altrimenti non saremmo arrivati qua
            freeParP1=(vecI[0][s]-Point[s])/(t[s]);
            freeParP2=(vecI[1][s]-Point[s])/(t[s]);
            freeParP3=(vecJ[0][s]-Point[s])/(t[s]);
            freeParP4=(vecJ[1][s]-Point[s])/(t[s]);
            break;
        }
    }

    // Riempimento del vector di array
    idpar.push_back({double(i), freeParP1});
    idpar.push_back({double(i), freeParP2});
    idpar.push_back({double(j), freeParP3});
    idpar.push_back({double(j), freeParP4});

    //ora lo devo ordinare rispetto a freeParP_
    BubbleSort_mod(idpar);

    // creo un dizionario per poter riottenere le informazioni sul punto sulla retta a partire dal parametro libero
    map<double, Vector3d> dizfreeParToVec;
    // popolo il dizionario (non mi interessa se ho una chiave che si ripete (quindi questa verrebbe sovrascritta) tanto il dizionario viene usato negli ultimi elseif)
    dizfreeParToVec[freeParP1] = vecI[0];
    dizfreeParToVec[freeParP2] = vecI[1];
    dizfreeParToVec[freeParP3] = vecJ[0];
    dizfreeParToVec[freeParP4] = vecJ[1];

    //estraiamo gli estremi delle tracce
    double id1 = idpar[0][0];
    double id2 = idpar[1][0];

    if (id1==id2){
        return 1;
    }
    else if (id1 != id2){

        // Memorizzo i dati nella struttura Traces
        trace.numTraces++;
        trace.IdsFractures.resize(trace.numTraces);
        trace.IdTraces.push_back(trace.numTraces-1);


        // Visualizza i dati della traccia
        /*cout << "\nTraccia " << trace.numTraces-1 << ":" << endl;
        cout << " - Frattura 1: " << i << endl;
        cout << " - Frattura 2: " << j << endl;*/

        array<unsigned int,2> vector_id_fractures = {};
        vector_id_fractures[0] = i;
        vector_id_fractures[1] = j;
        trace.IdsFractures[trace.numTraces - 1]=vector_id_fractures;

        //determino gli estremi

        //creo la matrice 3 righe 2 colonne da inserire nel vettore delle matrici degli estremi.
        Matrix<double, 3, 2> Estremi;

        // Controllo che le chiavi esistano nel dizionario prima di accedervi
        if (dizfreeParToVec.find(idpar[1][1]) != dizfreeParToVec.end() &&
            dizfreeParToVec.find(idpar[2][1]) != dizfreeParToVec.end()) {

            // Riempio la matrice con i punti dai valori del dizionario
            Estremi.col(0) = dizfreeParToVec[idpar[1][1]];
            Estremi.col(1) = dizfreeParToVec[idpar[2][1]];
        } else {
            cerr << "Chiavi non trovate nel dizionario!" << endl;
            return false; // O gestisci l'errore in modo appropriato
        }

        trace.CoordinatesEstremiTraces.push_back(Estremi);
        ///cout << " - Estremi traccia: (" << Estremi.col(0).transpose() << "), (" << Estremi.col(1).transpose() << ")" << endl;
        trace.lengthTraces.push_back(euclidean_distance(Estremi.col(0), Estremi.col(1)));
        ///cout<<"trace length: "<<euclidean_distance(Estremi.col(0), Estremi.col(1))<<endl;
    }else
    {
        return 2; //controllo se ci sono casistiche non considerate
    }


    //controllo PASSANTi O NO?

    // cout << Tips_Shy(fracture, trace, i, j)<<endl;

    std::chrono::steady_clock::time_point t_begin= chrono::steady_clock::now();

    int pass = 0;

    //evitiamo cancellazione numerica con la sottrazione
    if (abs(idpar[0][1]- idpar[1][1]) < tolDefault && abs(idpar[2][1]- idpar[3][1]) < tolDefault){ //passante per entrambe le fratture

        pass = 0;
        inserimento_map(pass,idpar[0][0], trace);
        inserimento_map(pass,idpar[1][0], trace);

        ///cout << "   Passante per entrambi le fratture " <<idpar[0][0]<<" "<<idpar[1][0]<< " ."<<endl;

        //lascio solo il primo per ricordarci che abbiamo fatto un controllo sul tempo di esecuzione e avere un'idea sul costo computazionale
        // poi l'abbiamo fatto per ogni condizione
        /*std::chrono::steady_clock::time_point t_end= chrono::steady_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();
        cout<<"\t ------- FREEPAR elapsed time: "<<duration<<" microseconds\n" <<endl;*/


        return 3;


    }else if (  (idpar[0][0] == double(j) &&  idpar[3][0] == double(j))
                || (idpar[0][0]==double(i) && idpar[2][0]==double(i) && abs(idpar[0][1]-idpar[1][1])<tolDefault)
                || (idpar[0][0]==double(j) && idpar[2][0]==double(j) && abs(idpar[2][1]-idpar[3][1])<tolDefault) ) // passante solo per i
    {

        pass = 0;
        inserimento_map(pass,i, trace);
        ///cout << "   Passante per la frattura " <<i<< " ."<<endl;
        pass = 1;
        inserimento_map(pass,j, trace);
        ///cout << "   NON Passante per la frattura " <<j<< " ."<<endl;

        return 4;


    }
    else if((idpar[0][0] == double(i) && idpar[3][0] == double(i))
              ||(idpar[0][0]== double(j) && idpar[2][0] == double(j) && abs(idpar[0][1]-idpar[1][1])<tolDefault)
               || (idpar[0][0]== double(i) && idpar[2][0] == double(i) && abs(idpar[2][1]-idpar[3][1])<tolDefault) ) // passante solo per j
    {
        pass = 0;
        inserimento_map(pass,idpar[1][0], trace);
        ///cout << "   Passante per la frattura " <<j<< " ."<<endl;
        pass = 1;
        inserimento_map(pass,idpar[0][0], trace);
        ///cout << "   NON Passante per la frattura " <<i<< " ."<<endl;

        return 4;
    }
    else{ //non passante per entrambe
        pass = 1;
        inserimento_map(pass,idpar[1][0], trace);
        inserimento_map(pass,idpar[0][0], trace);
        ///cout << "   NON Passante per entrambi le fratture " <<idpar[0][0]<<" "<<idpar[1][0]<< " ."<<endl;

        return 6;
    }

    return 0;
}



bool Tips_Shy(Fractures& fracture, Traces& trace, const int i, const int j)
{
    auto t_begin = chrono::steady_clock::now();
    const double tolerance = tolDefault;  // Tolleranza unificata per i confronti

    auto is_passing = [&](const MatrixXd& vertice, const Matrix<double, 3, 2>& trace_coords) {
        int touchCount = 0;
        double dist1 = 0;
        double dist2 = 0;

        for (int k = 0; k < vertice.cols(); ++k) {
            Vector3d v1 = vertice.col(k);
            Vector3d v2 = vertice.col((k + 1) % vertice.cols());
            double edgeLength = (v1 - v2).norm();

            for (int e = 0; e < 2; ++e) {
                Vector3d extremity = trace_coords.col(e);
                dist1 = (extremity- v1).norm();
                dist2 = (extremity- v2).norm();
                if (abs(dist1 + dist2 - edgeLength) <= tolerance) {
                    touchCount++;
                    break;
                }
            }
        }
        return (touchCount == 2);
    };

    bool passa_i = is_passing(fracture.CoordinatesVertice[i], trace.CoordinatesEstremiTraces[trace.numTraces-1]);
    //cout << " - Passante per la frattura " << i << " : " << (passa_i ? "Sì" : "No") << endl;

    bool passa_j = is_passing(fracture.CoordinatesVertice[j], trace.CoordinatesEstremiTraces[trace.numTraces-1]);
    //cout << " - Passante per la frattura " << j << " : " << (passa_j ? "Sì" : "No") << endl;

    auto t_end = chrono::steady_clock::now();
    double duration = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();
    //cout << "\t ------- Tips_sHy elapsed time: " << duration << " microseconds" << endl;

    return passa_i && passa_j;
}


unsigned int Calcolo_par(Vector3d& t, Vector3d& Point, int i, vector<Vector3d>& vec, Fractures& fracture)
{

    auto t_begin = chrono::steady_clock::now();

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
            if (abs(V2[i]-V1[i])>tolDefault)
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
        if (freeParP0 >= freeParV1 - tolDefault && freeParP0 <= freeParV2 + tolDefault) {

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
}
