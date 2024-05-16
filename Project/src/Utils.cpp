#include "Utils.hpp"
#include "FracturesTraces.hpp"
#include <vector>
#include "Eigen/Eigen"
#include <cmath> // per sqrt
#include <algorithm> // qui si trova std::max_element
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono> //for counting time

//distanza massima
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices)
{

    vector<double> distance_vector(num_vertices-1, 0.0);
    // itero sulle colonne
    for(unsigned int i = 0; i < num_vertices - 1; ++i){
        distance_vector[i] = sqrt( (m(0,i) - m(0,i+1))*(m(0,i) - m(0,i+1)) +
                                   (m(1,i) - m(1,i+1))*(m(1,i) - m(1,i+1)) +
                                   (m(2,i) - m(2,i+1))*(m(2,i) - m(2,i+1)) );

        cout <<"*: "<< distance_vector[i]<<endl;
        distance_vector[i] = sqrt(pow((m(0,i) - m(0,i+1)),2) + pow((m(1,i) - m(1,i+1)),2) + pow((m(2,i) - m(2,i+1)),2));
        cout <<"^: "<< distance_vector[i]<<endl;

    }

    distance_vector[num_vertices-1]=sqrt(pow((m(0,num_vertices-1) - m(0,0)),2) +
                                             pow((m(1,num_vertices-1) - m(1,0)),2) +
                                             pow((m(2,num_vertices-1) - m(2,0)),2));

    auto it_max_distance = max_element(distance_vector.begin(), distance_vector.end());
    double max_distance = *it_max_distance;
    it_max_distance = distance_vector.begin(); //per inizializzare correttamente ad ogni richiamo della funzione

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
    distance_bar = sqrt(pow( bar1[0] - bar2[0],2) +
                        pow( bar1[1] - bar2[1],2) +
                        pow( bar1[2] - bar2[2],2));

    double max_distance = 0.0;
    max_distance = (l1 + l2) * 0.5;

    if (distance_bar > max_distance){
        return false; // le fratture non si intersecano
    }

    return true;
}

//normale al piano
array<double, 3> normal_vector(const MatrixXd& m)
{

    array <double, 3> v1 = {}; // vettore1
    array <double, 3> v2 = {}; // vettore2
    array <double, 3> v = {};

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

        array <double, 3> v_normal;
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
    int countinter = 0;

    cout << "Ci sono "<<(fracture.NumFractures * (fracture.NumFractures - 1)) * 0.5 <<" possibili intersezioni"<<endl;

    for (unsigned int i = 0; i< fracture.NumFractures - 1 ; i++ ){
        for (unsigned int j=i+1; j < fracture.NumFractures; j++ ){
            // richiamo funzione sfera
            if ( !check_sphere( fracture.baricentro[i], fracture.baricentro[j], fracture.lenghtMaxEdges[i], fracture.lenghtMaxEdges[j]) )
            {
                escluse++;
                continue;
            }

            // t è la tangente che inidividua la direzione della retta di intersezione
            array <double,3> t = {};
            t = vec_product(fracture.vettoreNormalePiano[i], fracture.vettoreNormalePiano[j]);

            //il nostro P0 è il baricentro (COMPUTATIONAL GEOMETRY 2, PROBLEMA 4)
            Vector3d Point = {};
            bool ris = system_solution(fracture.vettoreNormalePiano[i], fracture.vettoreNormalePiano[j],
                                    fracture.baricentro[i], fracture.baricentro[j],
                                    t, Point);
            //tutto questo ci serve per trovare la retta di intersezione r(x) = x*t_ + Point;

            if (ris == 0){
                cout<<"non c'è intersezione tra frattura "<<i<<" e "<<j<<endl;
            }else{
                //cout <<"possibile intersezione tra "<<i<<" e "<<j<<endl;
            }




            MatrixXd matrixVerticesI = fracture.CoordinatesVertice[i];
            unsigned int iterI = 0;
            int numColonneI = matrixVerticesI.cols();
            vector<double> vecI;
            vecI.reserve(numColonneI);// RICORDARSI DI RESETTARLO ALLA FINE QUANDO NON SERVE PIu con .empty()
            for (int z = 0; z < numColonneI; z++) {
                Vector3d V1 = matrixVerticesI.col(z);
                Vector3d V2;
                if (z == numColonneI - 1) { // Accoppio l'ultimo vertice con il primo
                    V2 = matrixVerticesI.col(0);
                } else {
                    V2 = matrixVerticesI.col(z + 1);
                }
            /// calcolo retta (GeometryComputational 1) dove giacciono i due vertici

                Vector3d Punto0 = {};

                bool a;
                a = soluzione_sistema3x2(t,V1,V2,Point,Punto0);


                /// calcolo intersezione tra retta di intersezione piani e retta dei due vertici
                /// ... (= Vector3d Punto0 = X_punto0,Y_punto0,Z_punto0)
                /// valuto Punto0, V1 e V2 nella retta dei due vertici e calcolo il valore che assume il parametro libero
                double freeParP0 = Punto0[0];  // finire formula inversa usando
                double freeParV1 = V1[0];   //finire formula inversa
                double freeParV2 = V2[0];   //finire formula inversa
                if ((freeParP0<=freeParV1 && freeParP0>=freeParV2)||(freeParP0<=freeParV2 && freeParP0>=freeParV1) )
                {
                    iterI = iterI + 1 ;
                    //vecI.push_back(Punto0);
                }
            }
            //cout<<"\t\tfinito il controllo tra "<<i<<" e "<<j<<endl;
            if (iterI != 2) // la frattura non può avere traccia con l'altra frattura, mi fermo
            {
                continue;

            }





        }//end for j

    }//end for i
    cout << "\nescluse in principio "<< escluse<< " possibili intersezioni! SBAM."<<endl;




}

}
