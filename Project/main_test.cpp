#include "gtest/gtest.h"
#include "src/Utils.hpp"
#include "src/FracturesTracesPolygons.hpp"
#include <sstream>
#include "TestFractures.hpp"

using namespace std;
using namespace Eigen;
using namespace GeometryLibrary;
int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
