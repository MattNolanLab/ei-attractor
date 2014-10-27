#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "gtest/gtest.h"
#include "nmda.h"

using namespace std;
using ::testing::InitGoogleTest;

namespace {

bool almost_equal(double x, double y, double rtol=1e-5) {
    return abs((x - y)) / y < rtol;
}

TEST(iaf_gridcells, test_nmda_multiplier) {
    string line;
    stringstream convert;
    double Vm = 0.;
    double C_Mg = 0.;
    double correct_s = 0.;
    const double rtol = 1e-10;

    ifstream fin("nmda_test_vectors.txt");
    ASSERT_TRUE(fin.is_open());

    // The first line is a header
    getline(fin, line);
    while (fin >> Vm >> C_Mg >> correct_s) {
        if (fin.good()) {
            //cout << Vm << " " << C_Mg << " " << correct_s << endl;
            double calc_s = nmda_mg_multiplier(Vm, C_Mg);
            ASSERT_TRUE(almost_equal(correct_s, calc_s, rtol)) <<
                "Values " << correct_s << " and " << calc_s <<
                " are not equal with rtol = " << rtol << endl;
        }
    }
}

}  // namespace


int main(int argc, char **argv) {
  InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
