#include "utils.hpp"
#include "FracturesTracesPolygons.hpp"
#include "inline.hpp"
#include <vector>
#include "Eigen/Eigen"
#include <cmath> // per sqrt
#include <vector>

double tolDefault = 100 * std::numeric_limits<double>::epsilon();

using namespace std;

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
    //controllo la distanza tra i due baricentri
    double distance_bar = 0.0;
    distance_bar = sqrt( (bar1[0] - bar2[0])*(bar1[0] - bar2[0]) +
                         (bar1[1] - bar2[1])*(bar1[1] - bar2[1]) +
                         (bar1[2] - bar2[2])*(bar1[2] - bar2[2]) );

    double max_distance = 0.0;
    max_distance = (l1 + l2) ;

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

// funzione implementata per risolvere un sistema 3x3 per trovare l’intersezione fra piani
bool intersezione_piani (Vector3d& n1,
                        Vector3d& n2,
                        array <double, 3>& b1,
                        array <double, 3>& b2,
                        Vector3d& t,
                        Vector3d& Point){
    double d1 = 0.0; //termine noto piano 1
    double d2 = 0.0; //termine noto piano 2
    // termine noto piano d'intersezione = 0 come prodotto scalare tra vettore normale
    // al piano e punto generico della frattura (scelto il baricentro)
    d1 = n1[0] * b1[0] + n1[1] * b1[1]  + n1[2] * b1[2];
    d2 = n2[0] * b2[0] + n2[1] * b2[1]  + n2[2] * b2[2];

    Vector3d vec = n1.cross(n2);
    // parametrizzo la z = t
    // t individua la direzione della retta di intersezione tra i due piani. Si trova un Point qualsiasi sulla retta,
    // quindi un suo parametro può essere scelto a piacere. Si sceglie la z del Point come parametro libero, ma si può fare se la t[2]≠0; dunque selezione diversa se t[2]=0
    // if (t[2]>1e-14)//se t[2] != 0
    //{
    // Crea la matrice A e la inizializza con i vettori  n1, n2, t
    Matrix3d A;
    A << n1[0], n1[1],n1[2],
        n2[0], n2[1],n2[2],
        t[0], t[1], t[2];
    //A.row(2) = t;

    // risoluzione sistema
    // det(A) != 0 : per verificare che il sistema abbia soluzione (i due piani si intersecano)
    // controllo che il modulo di t sia diverso da 0 per evitare prodotto vettoriale != 0

    // Verifica della possibilità di intersezione
    if (vec.dot(t) > tolDefault) {
        // Risoluzione del sistema lineare per trovare il punto di intersezione
        Vector3d B;
        B << d1, d2,0;
        Point = A.colPivHouseholderQr().solve(B); //(x,y,t)
        return true;

    } else {
        // I piani sono paralleli o coincidenti
        return false;
    }

    /*}else{//se t[2] == 0
        // scelta di un altro parametro libero
        Matrix2d A;
        A << n1[0], n1[2],
            n2[0], n2[2];

        //A.row(2) = t;

        // risoluzione sistema
        // det(A) != 0 :

        // Verifica della possibilità di intersezione
        if (vec.dot(t) > 1e-14) {
            // Risoluzione del sistema lineare per trovare il punto di intersezione
            Vector2d B;
            B << d1, d2;
            Vector2d duecompPoint = A.colPivHouseholderQr().solve(B); //(x,y,t)
            Point[0] = duecompPoint[0];
            Point[2] = duecompPoint[1];
            Point[1] = 0;
            return true;

        } else {
            // I piani sono paralleli o coincidenti
            return false;
        }
    }*/
}


//troviamo il punto di intersezione tra la retta (Point, t)  e la retta di prolungamento del segmento V1V2
bool intersezione_rette (Vector3d& t,
                        Vector3d& V1,
                        Vector3d& V2,
                        Vector3d& Point,
                        Vector3d& Punto0)
{

    Vector3d vettoreDirezioneI = V1 - V2;
    VectorXd termine_noto = V1 - Point;
    /*if (abs(t[2])< 1e-14){

        MatrixXd A_(2,2);
        A_.col(0) << t[0], t[1];
        A_.col(1) << vettoreDirezioneI[0], vettoreDirezioneI[1];

        MatrixXd Completa_(2,3); //è A con ultima colonna termine noto
        Completa_.col(0) << A_.col(0);
        Completa_.col(1) << A_.col(1);
        Completa_.col(2) << termine_noto[0],termine_noto[1];

        // Calcola il rango utilizzando la decomposizione LU
        int rankA = A_.fullPivLu().rank();

        int rankC = Completa_.fullPivLu().rank();

        if (rankA == rankC && rankA == 2) {
            VectorXd sol_;
            if (A_.rows() >= A_.cols()) {

            std::chrono::steady_clock::time_point t_begin= chrono::steady_clock::now();
            sol = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(termine_noto);
            std::chrono::steady_clock::time_point t_end= chrono::steady_clock::now();
            double duration = chrono::duration_cast<chrono::microseconds>(t_end - t_begin).count();
            cout<<"\tjacobisvd.solve time: "<<duration<<" microseconds\n" <<endl;


                //chrono::steady_clock::time_point t_begin2= chrono::steady_clock::now();
                //sol_ = (A_.transpose() * A_).ldlt().solve(A_.transpose() * termine_noto);
                //chrono::steady_clock::time_point t_end2= chrono::steady_clock::now();
                //double duration2 = chrono::duration_cast<chrono::microseconds>(t_end2 - t_begin2).count();
                //cout<<"\tother time: "<<duration2<<" microseconds\n" <<endl;

                //ora abbiamo trovato i coeffieiente s e aplha

                // Calcola il punto di intersezione Punto0
                //Punto0 = Point + sol_[0] * t ;
                //return true;
            //}
        //}
    //}
    //else{ //terza componente di t != 0*/
    // si risolve un problema 2x2 escludendo la componente minima di t
    int pos = 0;
    if(abs(t[1]) < abs(t[0]) + tolDefault){
        pos = 1;
    }
    if(abs(t[2]) < (t[pos]) + tolDefault){
        pos = 2;
    }

    Matrix2d A;
    Vector2d b;
    if (pos == 0) {
        A << t[1], vettoreDirezioneI[1],
            t[2], vettoreDirezioneI[2];
        b << termine_noto[1], termine_noto[2];
    } else if (pos == 1) {
        A << t[0], vettoreDirezioneI[0],
            t[2], vettoreDirezioneI[2];
        b << termine_noto[0], termine_noto[2];
    } else {
        A << t[0], vettoreDirezioneI[0],
            t[1], vettoreDirezioneI[1];
        b << termine_noto[0], termine_noto[1];
    }

    Vector3d vettopa = vec_product(t,vettoreDirezioneI);
    double norm = vettopa.norm();

    if (norm>tolDefault) {
        Vector2d sol;
        sol = A.colPivHouseholderQr().solve(b);

        // Calcola il punto di intersezione Punto0
        Punto0 = Point + sol[0] * t ;
        return true;

    }

    return false;

}

