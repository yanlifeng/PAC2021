#include <stdio.h>
#include "util.h"
#include "mrc.h"
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char **argv) {
    if (argc != 3) {
        cout << "Invalid arguments!" << endl;
        cout << "Usage: ./validate mrc1_path mrc_2path" << endl;
        return 0;
    }
    string mrc1_path = argv[1];
    string mrc2_path = argv[2];
    cout << "MRC1 Path:" << mrc1_path << endl;
    cout << "MRC2 Path:" << mrc2_path << endl;

    MRC mrc1(mrc1_path.c_str(), "rb");
    MRC mrc2(mrc2_path.c_str(), "rb");

    if (mrc1.getNx() != mrc2.getNx()
        || mrc1.getNy() != mrc2.getNy()
        || mrc1.getNz() != mrc2.getNz()
            ) {
        cout << "MRCs' size don't match!" << endl;
    }

    float *mrc1_data = new float[mrc1.getNx() * mrc1.getNy()];
    float *mrc2_data = new float[mrc2.getNx() * mrc2.getNy()];

    float mrc12error = 0;
    double mrc12errord = 0;
    float mrc1sum = 0;
    float mrc2sum = 0;

    float mrc1mean = mrc1.getMean();
    float mrc2mean = mrc2.getMean();
    cout << "mrc1mean: " << mrc1mean << endl;
    cout << "mrc1min: " << mrc1.getMin() << endl;
    cout << "mrc1max: " << mrc1.getMax() << endl;
    cout << "mrc2mean: " << mrc2mean << endl;
    cout << "mrc2min: " << mrc2.getMin() << endl;
    cout << "mrc2max: " << mrc2.getMax() << endl;


    int out_cnt = 0;

    for (int i = 0; i < mrc1.getNz(); i++) {
        mrc1.read2DIm_32bit(mrc1_data, i);
        mrc2.read2DIm_32bit(mrc2_data, i);

        for (int j = 0; j < mrc1.getNx() * mrc1.getNy(); j++) {
            mrc1sum += mrc1_data[j];
            mrc2sum += mrc2_data[j];
//            if (fabs(mrc1_data[j] - mrc2_data[j]) > 1e-7 && out_cnt < 100) {
//                printf("%.6f %.6f\n", mrc1_data[j], mrc2_data[j]);
//                out_cnt++;
//            }
            mrc12error += fabs(mrc1_data[j] - mrc2_data[j]);
            mrc12errord += fabs(mrc1_data[j] - mrc2_data[j]);
        }

    }

    float absMeanError = mrc12error / (mrc1.getNx() * mrc1.getNy() * mrc1.getNz());
    float rltMeanError = mrc12error / (mrc1.getNx() * mrc1.getNy() * mrc1.getNz()) / fabs(mrc1.getMean());

    double absMeanErrord = mrc12errord / (mrc1.getNx() * mrc1.getNy() * mrc1.getNz());
    double rltMeanErrord = mrc12errord / (mrc1.getNx() * mrc1.getNy() * mrc1.getNz()) / fabs(mrc1.getMean());


    cout << "range " << mrc1.getNx() * mrc1.getNy() * mrc1.getNz() << endl;

    cout << "The error accumulation of the two volume is: " << mrc12error << endl;
    cout << "The error accumulation of the two volume is(double): " << mrc12errord << endl;
    cout << "The summation error of the two volume is: " << fabs(mrc1sum - mrc2sum) << endl;
    cout << "The mean error of the two volume is: " << fabs(mrc1mean - mrc2mean) << endl;
    cout << "The absolute mean error on a single voxel is: " << absMeanError << endl;
    cout << "The relative mean error on a single voxel is: " << rltMeanError << endl;
    cout << "The absolute mean error on a single voxel is(double): " << absMeanErrord << endl;
    cout << "The relative mean error on a single voxel is(double): " << rltMeanErrord << endl;

    if (fabs(rltMeanError) < 1e-3) {
        cout << "Validation Passed!" << endl;
    } else {
        cout << "Validation Failed!" << endl;
    }

    delete[] mrc1_data;
    delete[] mrc2_data;

    return 0;
}
