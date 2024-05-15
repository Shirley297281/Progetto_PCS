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

//distanza massima
double max_euclidean_distance(const MatrixXd& m, unsigned int num_vertices)
{

    vector<double> distance_vector(num_vertices-1, 0.0);
    // itero sulle colonne
    for(unsigned int i = 0; i < num_vertices - 1; ++i){
        distance_vector[i] = sqrt(pow((m(0,i) - m(0,i+1)),2) + pow((m(1,i) - m(1,i+1)),2) + pow((m(2,i) - m(2,i+1)),2));
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
vector <double> barycenter(const MatrixXd& m, unsigned int num_vertices)
{

    vector<double> barycenter_coord(3, 0.0);

    for(unsigned int i = 0; i < 3; ++i){
        for(unsigned int j = 0; j < num_vertices ; ++j){
            barycenter_coord[i] += m(i,j);
        }
        barycenter_coord[i] /= num_vertices;
    }
    return barycenter_coord;
}

//funzione sfere
bool check_sphere(const vector<double> bar1, const vector<double> bar2, const double l1, const double l2)
{
    //controlliamo la distanza tra i due baricentri
    double distance_bar = 0.0;
    distance_bar = sqrt(pow( bar1[0] - bar2[0],2) +
                        pow( bar1[1] - bar2[1],2) +
                        pow( bar1[2] - bar2[2],2));

    double max_distance = 0.0;
    max_distance = (l1 + l2) / 2.0;

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

    if (file.fail()) {
        cout << "Errore durante l'apertura del file." << endl;
        return false;
    }

    cout << "File aperto correttamente." << endl;

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

        vector <double> bary;
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
    cout << "Dati delle fratture memorizzati:" << endl;
    for (size_t i = 0; i < fracture.IdFractures.size(); ++i) {
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
    }



    return true;

}



void CalcoloTracce(Fractures& fracture, Traces& trace){
    int escluse = 0;

    for (unsigned int i = 0; i< fracture.NumFractures - 1 ; i++ ){
        for (unsigned int j=i+1; j < fracture.NumFractures; j++ ){
            // richiamo funzione sfera
            if ( !check_sphere( fracture.baricentro[i], fracture.baricentro[j], fracture.lenghtMaxEdges[i], fracture.lenghtMaxEdges[j]) )
            {
                escluse++;
                continue;
            }
            else{
                //cout<<"\n\n Forse c'è intersezione tra la frattura  "<<i<< " e la frattura "<<j<<endl;
                //lavorare con le tracce
            }

        }
    }

    cout << "\n\nescluse in principio "<< escluse<< " possibili intersezioni! SBAM."<<endl;
}

}
