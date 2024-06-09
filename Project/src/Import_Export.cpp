#include "Utils.hpp"
#include "FracturesTraces.hpp"
#include "namespace.hpp"//contiene gli header di tutte le funzione definite come GeometryLibrary
#include <vector>
#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include "inline.hpp"

//facciamo anche l'export??

//using namespace Eigen;
//using namespace std;

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

bool exportFR1(const string &filename, const Traces& trace)
{
    // #include <fstream> in main
    ofstream outFile("1TRACES_"+filename);
    if (!outFile) {
        cerr << "Errore nell'apertura del file per la scrittura." << endl;
        return 1;
    }
    outFile << "#Number of traces" << endl;
    outFile << trace.numTraces << endl;
    outFile << "#TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;
    for (unsigned int i = 0; i < trace.numTraces; i++)
    {
        outFile  << trace.IdTraces[i] <<"; ";
        outFile  << trace.IdsFractures[i][0] <<"; ";
        outFile  << trace.IdsFractures[i][1];

        Matrix<double, 3, 2> vertices = trace.CoordinatesEstremiTraces[i];
        Vector3d V1 = vertices.col(0);
        Vector3d V2 = vertices.col(1);
        for (unsigned int k = 0; k < 3; k++)
        {
            outFile << "; " << V1[k];
        }
        for (unsigned int k = 0; k < 3; k++)
        {
            outFile << "; " << V2[k];
        }
        outFile<< endl;
    }
    // Chiudere il file
    outFile.close();
    return 8;
}


bool secondoOutput(const string &filename, const Fractures& fracture, Traces& trace)
{
    ofstream outFile("2TRACESforeachFRECTURE_"+filename);
    if (!outFile) {
        cerr << "Errore nell'apertura del file per la scrittura." << endl;
        return false;
    }


    for (unsigned int i=0; i<fracture.NumFractures; i++)
    {


        vector<unsigned int> VecIdTracesPassI = {};
        vector<unsigned int> VecIdTracesNoPassI ={};
        unsigned int dimPass = 0;
        unsigned int dimNoPass = 0;



        if(!trace.TraceIdsPassxFracture[i].empty()){

           VecIdTracesPassI = trace.TraceIdsPassxFracture[i]; //estraggo il vettore di tracce passanti della frattura i
           dimPass =  VecIdTracesPassI.size();

        }
        if (!trace.TraceIdsNoPassxFracture[i].empty()){
            VecIdTracesNoPassI = trace.TraceIdsNoPassxFracture[i];  //estraggo il vettore di tracce non passanti della frattura i
            dimNoPass =  VecIdTracesNoPassI.size();
        }

        if(dimPass == 0 && dimNoPass == 0){
            continue;
        }



        unsigned int numTracesXi =  dimPass + dimNoPass;
        outFile << "#FractureId; NumTraces" << endl;
        outFile << i << "; " << numTracesXi << endl;
        outFile << "#TraceId; Tips; Length" <<endl;

        // Id traccia in ordine di lunghezza
        vector<array<double,2>> idLengthsPass = {};  // costruisco vettori di array [idtraccia, lunghezzatraccia]

        vector<array<double,2>> idLengthsNoPass = {};


        if (dimPass != 0) //altrimenti non faccio nulla
        {
            idLengthsPass.reserve(dimPass);
            bool tips = false; // false se passante
            for (unsigned int j=0; j<dimPass; j++) //riempio idLengthsPass
            {
                unsigned int tracciaIdPass = VecIdTracesPassI[j];
                double lunghezza = trace.lengthTraces[tracciaIdPass]; //estraggo la sua lunghezza tramite la posizione
                array<double,2> ArrayDiSupporto = {double(tracciaIdPass), lunghezza}; // {id, lunghezza corrispondente a quell'id}
                idLengthsPass.push_back(ArrayDiSupporto);
            }
            // uso bubble sort già creato apposta per ordinare secondo il secondo elemento (lunghezza)
            BubbleSort_mod(idLengthsPass);

            for (unsigned int j = 0; j < dimPass; j++)
            {
                // printo da idLenghsPass a sto punto ordinato
                outFile << idLengthsPass[dimPass-j-1][0] << "; "<< tips << "; "<< idLengthsPass[dimPass-j-1][1]<<endl; // CREARE VARIABILE TIPS, è specificatamente richiesto e sostituire passante con "fracture.Tips[idLengthsPass[i][0]]"
                trace.TraceIdsPassxFracture[i][j] = idLengthsPass[dimPass-j-1][0];
            }
        }
        if (dimNoPass != 0)
        {
            idLengthsNoPass.reserve(dimNoPass);
            bool tips = true; // false se non passante

            for (unsigned int j=0; j<dimNoPass; j++)
            {
                unsigned int tracciaIdNoPass = VecIdTracesNoPassI[j] ;
                double lunghezza = trace.lengthTraces[tracciaIdNoPass];
                array<double,2> ArrayDiSupporto = {double(tracciaIdNoPass), lunghezza}; // {id, lunghezza corrispondente a quell'id}
                idLengthsNoPass.push_back(ArrayDiSupporto);
            }
            // uso bubble sort già creato apposta per ordinare secondo il secondo elemento
            BubbleSort_mod(idLengthsNoPass);

            for (unsigned int j = 0; j < dimNoPass; j++)
            {
                outFile << idLengthsNoPass[dimNoPass - j-1][0] << "; " << tips << "; " << idLengthsNoPass[dimNoPass - j-1][1] << endl; // sostituire con fracture.Tips[idLengthsPass[i][0]]
                trace.TraceIdsNoPassxFracture[i][j] = idLengthsNoPass[dimNoPass-j-1][0];
            }

        }
    }

    // Funzione per stampare il contenuto di una mappa
    cout <<"\n\nper il check post riordinamento della TestId%PassxFracture\n";
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

    return true;
}


}
