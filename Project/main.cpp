#include "src/Utils.hpp"
#include "src/FracturesTraces.hpp"
#include <sstream>

using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;

int main()
{
    Fractures fracture;
    string filename = "FR3_data.txt";

    if( !ImportFR(filename, fracture) )
        return 1;

    // richiamo funzione calcola tracce
    Traces trace;
    CalcoloTracce(fracture, trace);

    return 0;
}
