/*******************************************************************
 *       Filename:  ReconstructionAlgo.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  07/07/2020 05:40:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ruan Huabin
 *          Email:  ruanhuabin@tsinghua.edu.cn
 *        Company:  Dep. of CS, Tsinghua Unversity
 *
 *******************************************************************/
#include "ReconstructionAlgo_WBP_RAM.h"
#include "mrc.h"
#include "CTF.h"
#include "math.h"
#include "fftw3.h"
#include "omp.h"
#include "util.h"
#include <sys/time.h>
#include <cassert>

#include "omp.h"

#include <immintrin.h>

#define mycos(x) (1-x*x/2+x*x*x*x/24)

const int threadNumber = 1;

const float eps = 1e-7;


double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

static void buf2fft(float *buf, float *fft, int nx, int ny) {
    int nxb = nx + 2 - nx % 2;
    int i;
    for (i = 0; i < (nx + 2 - nx % 2) * ny; i++) {
        fft[i] = 0.0;
    }
    for (i = 0; i < ny; i++) {
        memcpy(fft + i * nxb, buf + i * nx, sizeof(float) * nx);
    }
}

static void fft2buf(float *buf, float *fft, int nx, int ny) {
    int nxb = nx + 2 - nx % 2;
    int i;
    for (i = 0; i < nx * ny; i++) {
        buf[i] = 0.0;
    }
    for (i = 0; i < ny; i++) {
        memcpy(buf + i * nx, fft + i * nxb, sizeof(float) * nx);
    }
}

static void buf2fft_padding_1D(float *buf, float *fft, int nx_orig, int nx_final, int ny_orig, int ny_final) {
    int nxb = nx_final + 2 - nx_final % 2;
    int nxp = nx_final - nx_orig + 1;
    int i, j;
#pragma omp parallel for num_threads(threadNumber)
    for (i = 0; i < (nx_final + 2 - nx_final % 2) * ny_final; i++) {
        fft[i] = 0.0;
    }
#pragma omp parallel for num_threads(threadNumber)
    for (i = 0; i < ny_orig; i++) {
        memcpy(fft + i * nxb, buf + i * nx_orig, sizeof(float) * nx_orig);
        for (j = nx_orig; j < nx_final; j++)    // padding for continuity (FFT本质是周期延拓后做DFT，将首尾连接以保证连续性)
        {
            fft[i * nxb + j] = (float(j - nx_orig + 1) / float(nxp)) * buf[i * nx_orig] +
                               (float(nx_final - j) / float(nxp)) * buf[(i + 1) * nx_orig - 1];
        }
    }
}

static void fft2buf_padding_1D(float *buf, float *fft, int nx_orig, int nx_final, int ny_orig, int ny_final) {
    int nxb = nx_final + 2 - nx_final % 2;
    int i;
#pragma omp parallel for num_threads(threadNumber)
    for (i = 0; i < nx_orig * ny_orig; i++) {
        buf[i] = 0.0;
    }
#pragma omp parallel for num_threads(threadNumber)
    for (i = 0; i < ny_orig; i++) {
        memcpy(buf + i * nx_orig, fft + i * nxb, sizeof(float) * nx_orig);
    }
}

static void filter_weighting(float *data, int Nx, float radial, float sigma) {
    int Nx_padding = int(Nx / 10);
    int Nx_final = Nx + Nx_padding;
    fftwf_plan plan_fft, plan_ifft;
    float *bufc = new float[Nx_final + 2 - Nx_final % 2];
    plan_fft = fftwf_plan_dft_r2c_1d(Nx_final, (float *) bufc, reinterpret_cast<fftwf_complex *>(bufc), FFTW_ESTIMATE);
    plan_ifft = fftwf_plan_dft_c2r_1d(Nx_final, reinterpret_cast<fftwf_complex *>(bufc), (float *) bufc, FFTW_ESTIMATE);
    buf2fft_padding_1D(data, bufc, Nx, Nx_final, 1, 1);
    fftwf_execute(plan_fft);

    int radial_Nx = int(floor((Nx_final) * radial));
    float sigma_Nx = float(Nx_final) * sigma;

    // loop: Nx_final+2-Nx_final%2 (all Fourier components)
    for (int i = 0; i < Nx_final + 2 - Nx_final % 2; i += 2)   // radial filtering
    {
        if (i == 0)    // DC
        {
            bufc[i] *= 0.2;
            bufc[i + 1] *= 0.2;
        } else if (i / 2 <= radial_Nx)  // radial
        {
            bufc[i] *= (i / 2);
            bufc[i + 1] *= (i / 2);
        } else    // Gaussian falloff
        {
            bufc[i] = bufc[i] * float(radial_Nx) *
                      exp(-float((i / 2 - radial_Nx) * (i / 2 - radial_Nx)) / (sigma_Nx * sigma_Nx));
            bufc[i + 1] = bufc[i + 1] * float(radial_Nx) *
                          exp(-float((i / 2 - radial_Nx) * (i / 2 - radial_Nx)) / (sigma_Nx * sigma_Nx));
        }
    }

    fftwf_execute(plan_ifft);
    fft2buf_padding_1D(data, bufc, Nx, Nx_final, 1, 1);
    for (int i = 0; i < Nx; i++)   // normalization
    {
        data[i] /= Nx;
    }
    fftwf_destroy_plan(plan_fft);
    fftwf_destroy_plan(plan_ifft);
    delete[] bufc;
}

static void filter_weighting_1D_many(float *data, int Nx, int Ny, float radial, float sigma) {
    int Nx_padding = int(Nx / 10);
    int Nx_final = Nx + Nx_padding;
    fftwf_plan plan_fft, plan_ifft;
    float *bufc = new float[(Nx_final + 2 - Nx_final % 2) * Ny];

    plan_fft = fftwf_plan_many_dft_r2c(1, &Nx_final, Ny, (float *) bufc, NULL, 1, (Nx_final + 2 - Nx_final % 2),
                                       reinterpret_cast<fftwf_complex *>(bufc), NULL, 1,
                                       (Nx_final + 2 - Nx_final % 2) / 2, FFTW_ESTIMATE);
    plan_ifft = fftwf_plan_many_dft_c2r(1, &Nx_final, Ny, reinterpret_cast<fftwf_complex *>(bufc), NULL, 1,
                                        (Nx_final + 2 - Nx_final % 2) / 2, (float *) bufc, NULL, 1,
                                        (Nx_final + 2 - Nx_final % 2), FFTW_ESTIMATE);

    buf2fft_padding_1D(data, bufc, Nx, Nx_final, Ny, Ny);

    fftwf_execute(plan_fft);

    int radial_Nx = int(floor((Nx_final) * radial));
    float sigma_Nx = float(Nx_final) * sigma;
    // loop: Ny (all Fourier components for y-axis)
#pragma omp parallel for num_threads(threadNumber)
    for (int j = 0; j < Ny; j++) {
        // loop: Nx_final+2-Nx_final%2 (all Fourier components for x-axis)
        for (int i = 0; i < Nx_final + 2 - Nx_final % 2; i += 2) {
            if (i == 0)    // DC
            {
                bufc[i + j * (Nx_final + 2 - Nx_final % 2)] *= 0.2;
                bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] *= 0.2;
            } else if (i / 2 <= radial_Nx)  // radial
            {
                bufc[i + j * (Nx_final + 2 - Nx_final % 2)] *= (i / 2);
                bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] *= (i / 2);
            } else    // Gaussian falloff
            {
                bufc[i + j * (Nx_final + 2 - Nx_final % 2)] =
                        bufc[i + j * (Nx_final + 2 - Nx_final % 2)] * float(radial_Nx) *
                        exp(-float((i / 2 - radial_Nx) * (i / 2 - radial_Nx)) / (sigma_Nx * sigma_Nx));
                bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] =
                        bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] * float(radial_Nx) *
                        exp(-float((i / 2 - radial_Nx) * (i / 2 - radial_Nx)) / (sigma_Nx * sigma_Nx));
            }
        }
    }

    fftwf_execute(plan_ifft);
    fft2buf_padding_1D(data, bufc, Nx, Nx_final, Ny, Ny);
#pragma omp parallel for num_threads(threadNumber)
    for (int i = 0; i < Nx * Ny; i++)   // normalization
    {
        data[i] /= Nx;
    }

    fftwf_destroy_plan(plan_fft);
    fftwf_destroy_plan(plan_ifft);
    delete[] bufc;
}

static void filter_weighting_2D(float *image, int Nx, int Ny, float radial, float sigma, float psi_rad) {
    int Nx_padding = 0;
    int Nx_final = Nx_padding + Nx;
    fftwf_plan plan_fft, plan_ifft;
    float *bufc = new float[(Nx_final + 2 - Nx_final % 2) * Ny];
    plan_fft = fftwf_plan_dft_r2c_2d(Ny, Nx_final, (float *) bufc, reinterpret_cast<fftwf_complex *>(bufc),
                                     FFTW_ESTIMATE);
    plan_ifft = fftwf_plan_dft_c2r_2d(Ny, Nx_final, reinterpret_cast<fftwf_complex *>(bufc), (float *) bufc,
                                      FFTW_ESTIMATE);
    buf2fft(image, bufc, Nx, Ny);
    fftwf_execute(plan_fft);

    float radial_pi = M_PI * radial, sigma_pi = M_PI * sigma;
    float fft_x[Nx_final], fft_y[Ny];
    // loop: Nx_final (all Fourier components for x-axis)
    for (int i = 0; i < Nx_final; i++) {
        fft_x[i] = 2 * M_PI / (Nx_final) * i;
        if (fft_x[i] > M_PI) {
            fft_x[i] = fft_x[i] - 2 * M_PI;
        }
    }
    // loop: Ny (all Fourier components for y-axis)
    for (int j = 0; j < Ny; j++) {
        fft_y[j] = 2 * M_PI / Ny * j;
        if (fft_y[j] > M_PI) {
            fft_y[j] = fft_y[j] - 2 * M_PI;
        }
    }

    // loop: Ny (all Fourier components for y-axis)
    for (int j = 0; j < Ny; j++) {
        // loop: Nx_final+2-Nx_final%2 (all Fourier components for x-axis)
        for (int i = 0; i < Nx_final + 2 - Nx_final % 2; i += 2)   // 2D radial filtering
        {
            float x_rot = fft_x[i / 2] * cos(-psi_rad) - fft_y[j] * sin(-psi_rad), y_rot =
                    fft_x[i / 2] * sin(-psi_rad) + fft_y[j] * cos(-psi_rad);
            if (fabs(x_rot) < 1e-6)    // DC
            {
                bufc[i + j * (Nx_final + 2 - Nx_final % 2)] *= (0.2 * fft_x[1] * cos(-psi_rad));
                bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] *= (0.2 * fft_x[1] * cos(-psi_rad));
            } else if (fabs(x_rot) <= radial_pi)   // radial
            {
                bufc[i + j * (Nx_final + 2 - Nx_final % 2)] *= x_rot;
                bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] *= x_rot;
            } else    // Gaussian falloff
            {
                bufc[i + j * (Nx_final + 2 - Nx_final % 2)] = bufc[i + j * (Nx_final + 2 - Nx_final % 2)] * radial_pi *
                                                              exp(-((x_rot - radial_pi) * (x_rot - radial_pi)) /
                                                                  (sigma_pi * sigma_pi));
                bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] =
                        bufc[i + 1 + j * (Nx_final + 2 - Nx_final % 2)] * radial_pi *
                        exp(-((x_rot - radial_pi) * (x_rot - radial_pi)) / (sigma_pi * sigma_pi));
            }
        }
    }

    fftwf_execute(plan_ifft);
    fft2buf(image, bufc, Nx, Ny);
    for (int i = 0; i < Nx * Ny; i++)   // normalization
    {
        image[i] = image[i] / (Nx * Ny);
    }
    fftwf_destroy_plan(plan_fft);
    fftwf_destroy_plan(plan_ifft);
    delete[] bufc;
}

static void
ctf_correction(float *image, int Nx, int Ny, CTF ctf, bool flip_contrast, float z_offset)   // z_offset in pixels
{
    fftwf_plan plan_fft, plan_ifft;

    float *bufc = new float[(Nx + 2 - Nx % 2) * Ny];
#pragma omp critical
    {
        plan_fft = fftwf_plan_dft_r2c_2d(Ny, Nx, (float *) bufc, reinterpret_cast<fftwf_complex *>(bufc),
                                         FFTW_ESTIMATE);
        plan_ifft = fftwf_plan_dft_c2r_2d(Ny, Nx, reinterpret_cast<fftwf_complex *>(bufc), (float *) bufc,
                                          FFTW_ESTIMATE);
    }
    buf2fft(image, bufc, Nx, Ny);
    fftwf_execute(plan_fft);
    // loop: Ny (all Fourier components for y-axis)
    for (int j = 0; j < Ny; j++) {
        // loop: Nx+2-Nx%2 (all Fourier components for x-axis)
        for (int i = 0; i < (Nx + 2 - Nx % 2); i += 2) {
            float ctf_now = ctf.computeCTF2D(i / 2, j, Nx, Ny, true, flip_contrast, z_offset);
            bufc[i + j * (Nx + 2 - Nx % 2)] *= ctf_now;
            bufc[(i + 1) + j * (Nx + 2 - Nx % 2)] *= ctf_now;
        }
    }

    fftwf_execute(plan_ifft);
    fft2buf(image, bufc, Nx, Ny);
    for (int i = 0; i < Nx * Ny; i++)   // normalization
    {
        image[i] = image[i] / (Nx * Ny);
    }
#pragma omp critical
    {
        fftwf_destroy_plan(plan_fft);
        fftwf_destroy_plan(plan_ifft);
    }

    delete[] bufc;
}


ReconstructionAlgo_WBP_RAM::~ReconstructionAlgo_WBP_RAM() {

}

void ReconstructionAlgo_WBP_RAM::doReconstruction(map<string, string> &inputPara, map<string, string> &outputPara) {
    cout << "Run doReconstruction() in ReconstructionAlgo_WBP_RAM" << endl;

    // Input
    map<string, string>::iterator it = inputPara.find("path");
    string path;
    if (it != inputPara.end()) {
        path = it->second;
        cout << "File path: " << path << endl;
    } else {
        cout << "No specifit file path, set default: ./" << endl;
        path = "./";
    }

    it = inputPara.find("input_mrc");
    string input_mrc;
    if (it != inputPara.end()) {
        input_mrc = path + "/" + it->second;
        cout << "Input file name: " << input_mrc << endl;
    } else {
        cerr << "No input file name!" << endl;
        abort();
    }
    MRC stack_orig(input_mrc.c_str(), "rb");
    if (!stack_orig.hasFile()) {
        cerr << "Cannot open input mrc stack!" << endl;
        abort();
    }

    it = inputPara.find("output_mrc");
    string output_mrc;
    if (it != inputPara.end()) {
        output_mrc = path + "/" + it->second;
        cout << "Output file name: " << output_mrc << endl;
    } else {
        cout << "No output file name, set default: tomo.rec" << endl;
        output_mrc = "tomo.rec";
    }

    it = inputPara.find("prfx");
    string prfx;
    if (it != inputPara.end()) {
        prfx = path + "/" + it->second;
        cout << "Prefix: " << prfx << endl;
    } else {
        cout << "No prfx, set default: tomo" << endl;
        prfx = path + "/" + "tomo";
    }

    it = inputPara.find("h");
    int h;
    if (it != inputPara.end()) {
        h = atoi(it->second.c_str());
        cout << "Reconstruction height: " << h << endl;
    } else {
        h = int(stack_orig.getNx() / 4);
        cout << "No height for reconstruction, set default (Nx/4): " << h << endl;
    }

    bool skip_ctfcorrection, skip_weighting, skip_3dctf;
    it = inputPara.find("skip_ctfcorrection");
    if (it != inputPara.end()) {
        skip_ctfcorrection = atoi(it->second.c_str());
        cout << "Skip CTF correction: " << skip_ctfcorrection << endl;
    } else {
        cout << "No skip_ctfcorrection, set default: 0" << endl;
        skip_ctfcorrection = 0;
    }
    it = inputPara.find("skip_3dctf");
    if (it != inputPara.end()) {
        skip_3dctf = atoi(it->second.c_str());
        cout << "Skip 3D-CTF: " << skip_3dctf << endl;
    } else {
        cout << "No skip_3dctf, set default: 0 (Perform 3D-CTF correction)" << endl;
        skip_3dctf = 0;
    }
    it = inputPara.find("skip_weighting");
    if (it != inputPara.end()) {
        skip_weighting = atoi(it->second.c_str());
        cout << "Skip weighting: " << skip_weighting << endl;
    } else {
        cout << "No skip_weighting, set default: 0 (Perform WBP)" << endl;
        skip_weighting = 0;
    }

    float weighting_radial = 0.05, weighting_sigma = 0.5;
    if (!skip_weighting) {
        cout << "Weighting parameters:" << endl;
        it = inputPara.find("weighting_radial");
        if (it != inputPara.end()) {
            weighting_radial = atof(it->second.c_str());
            cout << "Radial: " << weighting_radial << endl;
        } else {
            cout << "No weighting_radial, set default: 0.05" << endl;
            weighting_radial = 0.05;
        }
        it = inputPara.find("weighting_sigma");
        if (it != inputPara.end()) {
            weighting_sigma = atof(it->second.c_str());
            cout << "Sigma: " << weighting_sigma << endl;
        } else {
            cout << "No weighting_sigma, set default: 0.5" << endl;
            weighting_sigma = 0.5;
        }
    }

    it = inputPara.find("input_tlt");
    string input_tlt;
    if (it != inputPara.end()) {
        input_tlt = path + "/" + it->second;
        cout << "Input tlt file name: " << input_tlt << endl;
    } else {
        cerr << "No input tlt file name!" << endl;
        abort();
    }
    FILE *ftlt = fopen(input_tlt.c_str(), "r");
    if (ftlt == NULL) {
        cerr << "Cannot open tlt file!" << endl;
        abort();
    }
    float theta[stack_orig.getNz()];
    float theta_max = 0.0;
    for (int n = 0; n < stack_orig.getNz(); n++) {
        fscanf(ftlt, "%f", &theta[n]);
        if (fabs(theta[n]) > theta_max) {
            theta_max = fabs(theta[n]);
        }
    }
    fflush(ftlt);
    fclose(ftlt);

    bool unrotated_stack = false;
    it = inputPara.find("unrotated_stack");
    if (it != inputPara.end()) {
        unrotated_stack = atoi(it->second.c_str());
        cout << "Input unrotated stack: " << unrotated_stack << endl;
    } else {
        cout << "No unrotated stack, input rotated stack" << endl;
        unrotated_stack = 0;
    }

    int h_tilt_max = int(ceil(
            float(stack_orig.getNx()) * sin(theta_max / 180 * M_PI) + float(h) * cos(theta_max / 180 * M_PI))) +
                     1;    // the maximum height after tilt


    string path_psi = path + "/" + "psi.txt";
    FILE *fpsi = fopen(path_psi.c_str(), "r");
    if (!fpsi) {
        cerr << "No psi found!" << endl;
        abort();
    }
    float psi_deg, psi_rad;
    fscanf(fpsi, "%f", &psi_deg);
    psi_rad = psi_deg * M_PI / 180;
    fflush(fpsi);
    fclose(fpsi);
    cout << "psi: " << psi_deg << endl;

    CTF ctf_para[stack_orig.getNz()];
    float Cs, pix, volt, w_cos;
    bool flip_contrast = false;
    int defocus_step = 1;
    it = inputPara.find("pixel_size");
    if (it != inputPara.end()) {
        pix = atof(it->second.c_str());
        cout << "Pixel size (A): " << pix << endl;
    } else {
        cerr << "No pixel size!" << endl;
        abort();
    }
    if (!skip_ctfcorrection) // read in defocus file for CTF correction
    {
        it = inputPara.find("Cs");
        if (it != inputPara.end()) {
            Cs = atof(it->second.c_str());
            cout << "Cs (mm): " << Cs << endl;
        } else {
            cerr << "No Cs!" << endl;
            abort();
        }
        it = inputPara.find("voltage");
        if (it != inputPara.end()) {
            volt = atof(it->second.c_str());
            cout << "Accelerating voltage (kV): " << volt << endl;
        } else {
            cerr << "No accelerating voltage!" << endl;
            abort();
        }
        it = inputPara.find("w");
        if (it != inputPara.end()) {
            w_cos = atof(it->second.c_str());
            cout << "Amplitude contrast: " << w_cos << endl;
        } else {
            cerr << "No amplitude contrast!" << endl;
            abort();
        }

        it = inputPara.find("flip_contrast");
        if (it != inputPara.end()) {
            flip_contrast = atoi(it->second.c_str());
            cout << "Flip contrast: " << flip_contrast << endl;
        } else {
            cout << "No flip contrast, set default: 0" << endl;
            flip_contrast = false;
        }

        it = inputPara.find("defocus_file");
        string defocus_file;
        if (it != inputPara.end()) {
            defocus_file = path + "/" + it->second;
            cout << "Defocus file name: " << defocus_file << endl;
        } else {
            cout << "No defocus file name, set default: defocus_file.txt" << endl;
            defocus_file = "defocus_file.txt";
        }
        FILE *fdefocus = fopen(defocus_file.c_str(), "r");
        if (!fdefocus) {
            cerr << "Cannot open defocus file!" << endl;
            abort();
        }

        for (int n = 0; n < stack_orig.getNz(); n++) {
            ctf_para[n].setN(n);
            ctf_para[n].setAllImagePara(pix, volt, Cs);
            float defocus_tmp[7];
            for (int i = 0; i < 7; i++)    // CTFFIND4 style
            {
                fscanf(fdefocus, "%f", &defocus_tmp[i]);
            }
            if (unrotated_stack) {
                ctf_para[n].setAllCTFPara(defocus_tmp[1], defocus_tmp[2], defocus_tmp[3], defocus_tmp[4], w_cos);
            } else {
                ctf_para[n].setAllCTFPara(defocus_tmp[1], defocus_tmp[2], defocus_tmp[3] - psi_deg, defocus_tmp[4],
                                          w_cos);   // 特别注意：目前的CTF估计结果取自CTFFIND4，是用原图（即未经旋转的图）估计的，因此对于重构旋转后的图，像散角（astig）也要旋转对应的角度！！！
            }
        }
        fflush(fdefocus);
        fclose(fdefocus);

        if (!skip_3dctf) {
            it = inputPara.find("defocus_step");
            if (it != inputPara.end()) {
                defocus_step = atoi(it->second.c_str());
                cout << "Defocus step (pixels): " << defocus_step << endl;
            } else {
                cout << "No defocus step, set default: 1" << endl;
                defocus_step = 1;
            }
        }
    }



    // Reconstruction
    cout << endl << "Reconstruction with (W)BP in RAM:" << endl << endl;
    if (!unrotated_stack)    // input rotated stack (y-axis as tilt axis)
    {
        cout << "Using rotated stack" << endl;

        float *stack_recon[stack_orig.getNy()]; // (x,z,y)
        for (int j = 0; j < stack_orig.getNy(); j++) {
            stack_recon[j] = new float[stack_orig.getNx() * h];
            for (int i = 0; i < stack_orig.getNx() * h; i++) {
                stack_recon[j][i] = 0.0;
            }
        }

        cout << "Start reconstruction:" << endl;
        float x_orig_offset = float(stack_orig.getNx()) / 2.0, z_orig_offset = float(h) / 2.0;
        // loop: Nz (number of images)

        double cost = 0;
        double cost0 = 0;
        double cost1 = 0;
        double cost2 = 0;
        double cost3 = 0;
        double t0, t1, t2, t3;

        t0 = GetTime();

        cout << "range " << stack_orig.getNx() << " " << stack_orig.getNy() << " " << stack_orig.getNz() << endl;

        fftwf_init_threads();


        for (int n = 0; n < stack_orig.getNz(); n++)   // loop for every micrograph
        {
            t1 = GetTime();
            t2 = GetTime();
            cout << "Image " << n << ":" << endl;
            float theta_rad = theta[n] / 180 * M_PI;
            double theta_rad_cos = cos(theta_rad);
            double theta_rad_sin = sin(theta_rad);
            float *image_now = new float[stack_orig.getNx() * stack_orig.getNy()];
            stack_orig.read2DIm_32bit(image_now, n);
            float *image_now_backup = new float[stack_orig.getNx() * stack_orig.getNy()];
            memcpy(image_now_backup, image_now, sizeof(float) * stack_orig.getNx() * stack_orig.getNy());

            t3 = GetTime();
            cost += t3 - t2;

            if (skip_ctfcorrection)  // no correction, simple (W)BP
            {
                cout << "\tSkip CTF correction, use raw micrograph!" << endl;

                // weighting
                if (!skip_weighting) {
                    cout << "\tStart weighting: " << endl;
                    filter_weighting_1D_many(image_now, stack_orig.getNx(), stack_orig.getNy(), weighting_radial,
                                             weighting_sigma);
                    cout << "\tDone!" << endl;
                }

                // reconstruction
                cout << "\tStart reconstruction, loop over xz-plane:" << endl;
                float *strip_now;
                float *recon_now;

                strip_now = new float[stack_orig.getNx()];
                recon_now = new float[stack_orig.getNx() * h];  // 第一维x，第二维z


                // loop: Ny (number of xz-slices)
                for (int j = 0; j < stack_orig.getNy(); j++) {
                    memcpy(strip_now, image_now + j * stack_orig.getNx(), sizeof(float) * stack_orig.getNx());

                    // BP
                    // loop: Nx*h (whole xz-slice)
                    for (int i = 0; i < stack_orig.getNx() * h; i++) {
                        recon_now[i] = 0.0;
                    }
                    // loop: Nx*h (whole xz-slice)
                    for (int k = 0; k < h; k++) {
                        for (int i = 0; i < stack_orig.getNx(); i++)   // loop for the xz-plane to perform BP
                        {
                            float x_orig = (float(i) - x_orig_offset) * cos(theta_rad) -
                                           (float(k) - z_orig_offset) * sin(theta_rad) + x_orig_offset;
                            float z_orig = (float(i) - x_orig_offset) * sin(theta_rad) +
                                           (float(k) - z_orig_offset) * cos(theta_rad) + z_orig_offset;
                            float coeff = x_orig - floor(x_orig);
                            if (floor(x_orig) >= 0 && ceil(x_orig) < stack_orig.getNx()) {
                                recon_now[i + k * stack_orig.getNx()] = (1 - coeff) * strip_now[int(floor(x_orig))] +
                                                                        (coeff) * strip_now[int(ceil(x_orig))];
                            } else {
                                recon_now[i + k * stack_orig.getNx()] = 0.0;
                            }
                        }
                    }
                    // loop: Nx*h (whole xz-slice)
                    for (int i = 0; i < stack_orig.getNx() * h; i++) {
                        stack_recon[j][i] += recon_now[i];
                    }
                }

                delete[] strip_now;
                delete[] recon_now;


                cout << "\tDone" << endl;
            } else    // perform CTF correction
            {
                if (skip_3dctf)  // perform simple CTF correction (no consideration of height)
                {
                    cout << "\tPerform uniform CTF correction!" << endl;

                    // correction
                    cout << "\tStart correction: " << endl;
                    ctf_correction(image_now, stack_orig.getNx(), stack_orig.getNy(), ctf_para[n], flip_contrast, 0.0);
                    cout << "\tDone!" << endl;

                    // weighting
                    if (!skip_weighting) {
                        cout << "\tStart weighting: " << endl;
                        filter_weighting_1D_many(image_now, stack_orig.getNx(), stack_orig.getNy(), weighting_radial,
                                                 weighting_sigma);
                        cout << "\tDone!" << endl;
                    }

                    // recontruction
                    cout << "\tStart reconstuction, loop over xz-plane:" << endl;

                    float *strip_now;
                    float *recon_now;

                    strip_now = new float[stack_orig.getNx()];
                    recon_now = new float[stack_orig.getNx() * h];  // 第一维x，第二维z


                    // loop: Ny (number of xz-slices)
                    for (int j = 0; j < stack_orig.getNy(); j++) {
                        memcpy(strip_now, image_now + j * stack_orig.getNx(), sizeof(float) * stack_orig.getNx());

                        // BP
                        // loop: Nx*h (whole xz-slice)
                        for (int i = 0; i < stack_orig.getNx() * h; i++) {
                            recon_now[i] = 0.0;
                        }
                        // loop: Nx*h (whole xz-slice)
                        for (int k = 0; k < h; k++) {
                            for (int i = 0; i < stack_orig.getNx(); i++)   // loop for the xz-plane to perform BP
                            {
                                float x_orig = (float(i) - x_orig_offset) * cos(theta_rad) -
                                               (float(k) - z_orig_offset) * sin(theta_rad) + x_orig_offset;
                                float z_orig = (float(i) - x_orig_offset) * sin(theta_rad) +
                                               (float(k) - z_orig_offset) * cos(theta_rad) + z_orig_offset;
                                float coeff = x_orig - floor(x_orig);
                                if (floor(x_orig) >= 0 && ceil(x_orig) < stack_orig.getNx()) {
                                    recon_now[i + k * stack_orig.getNx()] =
                                            (1 - coeff) * strip_now[int(floor(x_orig))] +
                                            (coeff) * strip_now[int(ceil(x_orig))];
                                } else {
                                    recon_now[i + k * stack_orig.getNx()] = 0.0;
                                }
                            }
                        }
                        // loop: Nx*h (whole xz-slice)
                        for (int i = 0; i < stack_orig.getNx() * h; i++) {
                            stack_recon[j][i] += recon_now[i];
                        }
                    }

                    delete[] strip_now;
                    delete[] recon_now;


                    cout << "\tDone" << endl;
                } else    // perform 3D-CTF correction
                {
                    cout << "\tPerform 3D-CTF correction!" << endl;

                    // write all corrected and weighted images in one mrc stack
                    cout << "\tPerform 3D correction & save corrected stack:" << endl;

                    t2 = GetTime();

                    int Nx = stack_orig.getNx();
                    int Ny = stack_orig.getNy();
                    int Nx2 = Nx + 2 - Nx % 2;
                    int Nxh = Nx / 2 + 1;
                    int Nyh = Ny / 2 + 1;
                    int Nyh2 = (Ny + 1) / 2;
//                    float *stack_corrected[int(h_tilt_max / defocus_step) + 1]; // 第一维遍历不同高度，第二维x，第三维y


                    float *stack_corrected = new float[(int(h_tilt_max / defocus_step) + 1) * Nx2 * Ny];


                    int n_zz = 0;
                    fftwf_plan_with_nthreads(4);


                    // weighting
                    if (!skip_weighting) {
                        cout << "\tStart weighting..." << endl;
//                        printf("now fftw use %d\n", fftw_planner_nthreads());

                        filter_weighting_1D_many(image_now, stack_orig.getNx(), stack_orig.getNy(), weighting_radial,
                                                 weighting_sigma);
                        cout << "\tDone" << endl;
                    }
                    t3 = GetTime();
                    cost0 += t3 - t2;

                    // 3D-CTF correction
                    //hotspot 1
                    cout << "\tStart 3D-CTF correction..." << endl;
                    t2 = GetTime();

                    //把第一波fft操作拿出来，预处理好放到bufc_pre中
                    float *image_pre;
                    image_pre = new float[Nx * Ny];
                    memcpy(image_pre, image_now, sizeof(float) * Nx * Ny);
                    float *bufc_pre = new float[(Nx + 2 - Nx % 2) * Ny];

                    fftwf_plan plan_fft_pre = fftwf_plan_dft_r2c_2d(Ny, Nx, (float *) bufc_pre,
                                                                    reinterpret_cast<fftwf_complex *>(bufc_pre),
                                                                    FFTW_ESTIMATE);
                    buf2fft(image_pre, bufc_pre, Nx, Ny);
                    fftwf_execute(plan_fft_pre);
                    fftwf_destroy_plan(plan_fft_pre);
                    fftwf_plan_with_nthreads(1);


                    int n_zz_thread[threadNumber];
                    for (int i = 0; i < threadNumber; i++) {
                        n_zz_thread[i] = 0;
                    }
                    int zl = -int(h_tilt_max / 2);
                    int zr = int(h_tilt_max / 2);

                    double a_w_cos = acos(w_cos);

                    float *atan_xy = new float[Ny * Nx];
                    double *sin_atan_xy = new double[Ny * Nx];
                    double *cos_atan_xy = new double[Ny * Nx];
                    for (int j = 1; j < Ny; j++) {
                        for (int i = 1; i < Nx; i++) {
                            float x_norm = i;
                            float y_norm = (j >= Nyh) ? (j - Ny) : (j);
                            float x_real = float(x_norm) / float(Nx) * (1 / pix);
                            float y_real = float(y_norm) / float(Ny) * (1 / pix);
                            float alpha = atan(y_real / x_real);
                            atan_xy[j * Nx + i] = alpha;
//                            sin_atan_xy[j * Nx + i] = sin(2 * alpha);
//                            cos_atan_xy[j * Nx + i] = cos(2 * alpha);

                        }
                    }


                    //loop: number of blocks (for 3D-CTF correction)

                    /*

init cost 0.00113
for cost 0.03127
fftw cost 0.00842
fft2buf cost 0.00098
last cost 0.00014
                     */

#pragma omp parallel for num_threads(threadNumber)
                    for (int zz = zl; zz < zr; zz += defocus_step) {
                        int thread_id = omp_get_thread_num();

                        int n_z = (zz + int(h_tilt_max / 2)) / defocus_step;
//                        stack_corrected[n_z] = new float[Nx2 * Ny];
                        float *image = stack_corrected + n_z * Nx2 * Ny;
                        CTF ctf = ctf_para[n];
                        float z_offset = float(zz) + float(defocus_step - 1) / 2;
                        memcpy(image, bufc_pre, sizeof(float) * (Nx2 * Ny));
                        float defocus1 = ctf.defocus1;
                        float defocus2 = ctf.defocus2;
                        float astig = ctf.astig;
                        float lambda = ctf.lambda;
                        float phase_shift = ctf.phase_shift;
                        float w_sin = ctf.w_sin;
                        float w_cos = ctf.w_cos;
                        float pix = ctf.pix;
                        float Cs = ctf.Cs;
                        float A = (defocus1 + defocus2 - 2 * z_offset * pix);
                        double sin2ast = sin(2 * astig);
                        double cos2ast = cos(2 * astig);
                        for (int j = 0; j < min(1, Ny); j++) {
                            for (int i = 0; i < Nxh; i++) {
                                float x = i, y = j;
                                float x_norm = x;
                                float y_norm = (y >= Nyh) ? (y - Ny) : (y);
                                float x_real = float(x_norm) / float(Nx) * (1 / pix);
                                float y_real = float(y_norm) / float(Ny) * (1 / pix);
                                float alpha;
                                if (x_norm == 0) {
                                    if (y_norm > 0) {
                                        alpha = M_PI_2;
                                    } else if (y_norm < 0) {
                                        alpha = -M_PI_2;
                                    } else {
                                        alpha = 0.0;
                                    }
                                } else {
                                    alpha = atan(y_real / x_real);
                                }
                                float freq2 = x_real * x_real + y_real * y_real;
                                float df_now = (A + (defocus1 - defocus2) * cos(2 * (alpha - astig))) / 2.0;
                                float chi = M_PI * lambda * df_now * freq2 -
                                            M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2 + phase_shift;
                                float ctf_now_tmp = w_sin * sin(chi) + w_cos * cos(chi);
                                if (flip_contrast) {
                                    ctf_now_tmp = -ctf_now_tmp;
                                }
                                if (ctf_now_tmp < 0) {
                                    image[i * 2 + j * Nx2] *= -1;
                                    image[(i * 2 + 1) + j * Nx2] *= -1;
                                }
                            }
                        }
                        for (int j = 0; j < Ny; j++) {
                            for (int i = 0; i < min(1, Nxh); i++) {
                                float x = i, y = j;
                                float x_norm = x;
                                float y_norm = (y >= Nyh) ? (y - Ny) : (y);
                                float x_real = float(x_norm) / float(Nx) * (1 / pix);
                                float y_real = float(y_norm) / float(Ny) * (1 / pix);
                                float alpha;
                                if (x_norm == 0) {
                                    if (y_norm > 0) {
                                        alpha = M_PI_2;
                                    } else if (y_norm < 0) {
                                        alpha = -M_PI_2;
                                    } else {
                                        alpha = 0.0;
                                    }
                                } else {
                                    alpha = atan(y_real / x_real);
                                }
                                float freq2 = x_real * x_real + y_real * y_real;
                                float df_now = (A + (defocus1 - defocus2) * cos(2 * (alpha - astig))) / 2.0;
                                float chi = M_PI * lambda * df_now * freq2 -
                                            M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2 + phase_shift;
                                float ctf_now_tmp = w_sin * sin(chi) + w_cos * cos(chi);
                                if (flip_contrast) {
                                    ctf_now_tmp = -ctf_now_tmp;
                                }
                                if (ctf_now_tmp < 0) {
                                    image[i * 2 + j * Nx2] *= -1;
                                    image[(i * 2 + 1) + j * Nx2] *= -1;
                                }
                            }
                        }


                        for (int j = 1; j < Ny; j++) {
//#pragma simd
                            for (int i = 1; i < Nxh; i++) {
                                float x_norm = i;
                                float y_norm = (j >= Nyh) ? (j - Ny) : (j);
                                float x_real = float(x_norm) / float(Nx) * (1 / pix);
                                float y_real = float(y_norm) / float(Ny) * (1 / pix);
                                float alpha = atan_xy[j * Nx + i];

                                float freq2 = x_real * x_real + y_real * y_real;
//                                float df_now = (A + (defocus1 - defocus2) * (cos_atan_xy[j * Nx + i] * cos2ast +
//                                                                             sin_atan_xy[j * Nx + i] * sin2ast)) / 2.0;
                                float df_now = (A + (defocus1 - defocus2) * cos(2 * (alpha - astig))) / 2.0;

//                                float df_now = (A + (defocus1 - defocus2) * mycos(2 * (alpha - astig))) / 2.0;
                                //TODO defocus1 - defocus2 is too small that remove it can even pass check
                                float chi = M_PI * lambda * df_now * freq2 -
                                            M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2
                                            + phase_shift;
//                                float ctf_now_tmp = w_sin * sin(chi) + w_cos * cos(chi);
                                double rrr = chi - a_w_cos;
                                int mul_base = 1;
                                //TODO check if rrr will bigger than 2*pi
                                if (fabs(rrr) > M_PI_2)mul_base = -1;
                                if (flip_contrast) {
                                    mul_base = -mul_base;
                                }
                                if (mul_base < 0) {
                                    image[i * 2 + j * Nx2] *= -1;
                                    image[(i * 2 + 1) + j * Nx2] *= -1;
                                }
                            }
                        }



//                        printf("now fftw use %d\n", fftw_planner_nthreads());

                        fftwf_plan plan_ifft;
#pragma omp critical
                        {
                            plan_ifft = fftwf_plan_dft_c2r_2d(Ny, Nx,
                                                              reinterpret_cast<fftwf_complex *>(image),
                                                              (float *) image, FFTW_ESTIMATE);
                        }

                        fftwf_execute(plan_ifft);
#pragma omp critical
                        {
                            fftwf_destroy_plan(plan_ifft);
                        }

                        for (int i = 0; i < Nx2 * Ny; i++)image[i] = image[i] / (Nx * Ny);
                        n_zz_thread[thread_id]++;
                    }
                    for (int i = 0; i < threadNumber; i++) {
                        n_zz += n_zz_thread[i];
                    }
                    delete[] bufc_pre;


                    t3 = GetTime();
                    printf("3DF cost %.3f\n", t3 - t2);
                    cost1 += t3 - t2;
                    cout << "\tDone!" << endl;



                    // recontruction

                    cout << "\tPerform reconstruction:" << endl;
                    t2 = GetTime();

                    // loop: Ny (number of xz-slices)
#pragma omp parallel for num_threads(threadNumber)
                    for (int j = 0; j < Ny; j++) {

                        // loop: Nx*h (whole xz-slice)
                        for (int k = 0; k < h; k++) {
                            double l = 0;
                            double r = Nx - 1;

                            l = max(l, ((k - z_orig_offset) * theta_rad_sin - x_orig_offset)
                                       / theta_rad_cos + x_orig_offset);
                            r = min(r, ((k - z_orig_offset) * theta_rad_sin - x_orig_offset + Nx - 1)
                                       / theta_rad_cos + x_orig_offset);
                            if (theta_rad_cos < 0)swap(l, r);

                            double ll = 0;
                            double rr = Nx - 1;
                            ll = max(ll, ((z_orig_offset - k) * theta_rad_cos - int(h_tilt_max / 2))
                                         / theta_rad_sin + x_orig_offset);

                            rr = min(rr, ((z_orig_offset - k) * theta_rad_cos - int(h_tilt_max / 2) +
                                          n_zz * defocus_step) / theta_rad_sin + x_orig_offset);
                            if (theta_rad_sin < 0)swap(ll, rr);
                            l = max(l, ll);
                            r = min(r, rr);

                            int li = ceil(l);
                            int ri = floor(r);
                            double A = (k - z_orig_offset) * theta_rad_sin - x_orig_offset;
                            double B = (k - z_orig_offset) * theta_rad_cos + z_orig_offset;
                            float C = -z_orig_offset + int(h_tilt_max / 2);

                            int i = li;
                            __m512d a_con = _mm512_set1_pd(A);
                            __m512d b_con = _mm512_set1_pd(B);
                            __m512d cos_con = _mm512_set1_pd(theta_rad_cos);
                            __m512d sin_con = _mm512_set1_pd(theta_rad_sin);
                            __m512 offset_con = _mm512_set1_ps(x_orig_offset);
                            __m512 eps_con = _mm512_set1_ps(1 - eps);
                            __m512 c_con = _mm512_set1_ps(C);
                            __m512 d_step_con = _mm512_set1_ps(1.0 / defocus_step);
                            __m512i idx = _mm512_set_epi32(li + 15, li + 14, li + 13, li + 12, li + 11, li + 10, li + 9,
                                                           li + 8, li + 7, li + 6, li + 5, li + 4, li + 3, li + 2,
                                                           li + 1, li + 0);
                            __m512i nxy_con = _mm512_set1_epi32(Nx2 * Ny);
                            __m512i nx2_con = _mm512_set1_epi32(Nx2);
                            __m512i j_con = _mm512_set1_epi32(j);
                            __m512i si_con = _mm512_set1_epi32(16);

                            __m512 one_con = _mm512_set1_ps(1);

                            static int cnt = 100;
                            bool checcc = 0;
                            for (; i + 16 <= ri + 1; i += 16) {



                                //float x_orig = (i - x_orig_offset) * theta_rad_cos - A;
                                //float z_orig = (i - x_orig_offset) * theta_rad_sin + B;

                                if (cnt <= 10) {
                                    int tmp[16];
                                    _mm512_store_epi32(tmp, idx);
                                    cout << "idx cal:" << endl;
                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
                                    cout << endl;
                                    cout << "idx std:" << endl;
                                    for (int ii = 0; ii < 16; ii++)printf("%d ", i + ii);
                                    cout << endl;
                                }

//                                cnt++;
                                //i - x_orig_offset
                                __m512 sub_tmp = _mm512_sub_ps(_mm512_cvtepi32_ps(idx), offset_con);

                                if (cnt <= 10) {
                                    float tmp[16];
                                    _mm512_store_ps(tmp, sub_tmp);
                                    cout << "sub_tmp cal:" << endl;
                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
                                    cout << endl;
                                    cout << "sub_tmp std:" << endl;
                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", i + ii - x_orig_offset);
                                    cout << endl;
                                }

                                //select pre 8 float to f0, last 8 to f1
                                __m256 f0 = _mm512_extractf32x8_ps(sub_tmp, 0);
                                __m256 f1 = _mm512_extractf32x8_ps(sub_tmp, 1);
                                __m512d d0 = _mm512_cvtps_pd(f0);
                                __m512d d1 = _mm512_cvtps_pd(f1);
//                                if (cnt <= 10) {
//                                    double tmp[8];
//                                    _mm512_store_pd(tmp, d0);
//                                    cout << "d0 cal:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "d0 std:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)printf("%.3f ", i + ii - x_orig_offset);
//                                    cout << endl;
//
//                                    _mm512_store_pd(tmp, d1);
//                                    cout << "d1 cal:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "d1 std:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)printf("%.3f ", i + ii + 8 - x_orig_offset);
//                                    cout << endl;
//                                }
                                // sub_tmp * theta_rad_cos - A
                                __m512d x_tmp0 = _mm512_sub_pd(_mm512_mul_pd(d0, cos_con), a_con);
                                __m512d x_tmp1 = _mm512_sub_pd(_mm512_mul_pd(d1, cos_con), a_con);
//                                if (cnt <= 10) {
//                                    double tmp[8];
//                                    _mm512_store_pd(tmp, x_tmp0);
//                                    cout << "x_tmp0 cal:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "x_tmp0 std:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)
//                                        printf("%.3f ", (i + ii - x_orig_offset) * theta_rad_cos - A);
//                                    cout << endl;
//                                    _mm512_store_pd(tmp, x_tmp1);
//                                    cout << "x_tmp1 cal:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "x_tmp1 std:" << endl;
//                                    for (int ii = 0; ii < 8; ii++)
//                                        printf("%.3f ", (i + ii + 8 - x_orig_offset) * theta_rad_cos - A);
//                                    cout << endl;
//                                }
                                //merge 8+8 double to 16 float
                                __m512 x_ori;
                                x_ori = _mm512_insertf32x8(x_ori, _mm512_cvtpd_ps(x_tmp0), 0);
                                x_ori = _mm512_insertf32x8(x_ori, _mm512_cvtpd_ps(x_tmp1), 1);

//                                if (cnt <= 10) {
//                                    float tmp[16];
//                                    _mm512_store_ps(tmp, x_ori);
//                                    cout << "x_ori cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "x_ori std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%.3f ", (i + ii - x_orig_offset) * theta_rad_cos - A);
//                                    cout << endl;
//                                }

                                // sub_tmp * theta_rad_sin + B
                                __m512d z_tmp0 = _mm512_add_pd(_mm512_mul_pd(d0, sin_con), b_con);
                                __m512d z_tmp1 = _mm512_add_pd(_mm512_mul_pd(d1, sin_con), b_con);
                                //merge 8+8 double to 16 float
                                __m512 z_ori;
                                z_ori = _mm512_insertf32x8(z_ori, _mm512_cvtpd_ps(z_tmp0), 0);
                                z_ori = _mm512_insertf32x8(z_ori, _mm512_cvtpd_ps(z_tmp1), 1);
//                                if (cnt <= 10) {
//                                    float tmp[16];
//                                    _mm512_store_ps(tmp, z_ori);
//                                    cout << "z_ori cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "z_ori std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%.3f ", (i + ii - x_orig_offset) * theta_rad_sin + B);
//                                    cout << endl;
//                                }
                                if (checcc) {
                                    float tmp1[16];
                                    _mm512_store_ps(tmp1, x_ori);
                                    float tmp2[16];
                                    _mm512_store_ps(tmp2, z_ori);
                                    for (int ii = 0; ii < 16; ii++) {
                                        float x_orig = (ii + i - x_orig_offset) * theta_rad_cos - A;
                                        float z_orig = (ii + i - x_orig_offset) * theta_rad_sin + B;
                                        if (fabs(x_orig - tmp1[ii]) > 1e-5) {
                                            printf("GG\n");
                                            printf("x_orig on %d %d %d\n", j, k, i + ii);
                                            printf("%.6f %.6f\n", x_orig, tmp1[ii]);
                                            exit(0);
                                        }
                                        if (fabs(z_orig - tmp2[ii]) > 1e-5) {
                                            printf("GG\n");
                                            printf("z_orig on %d %d %d\n", j, k, i + ii);
                                            printf("%.6f %.6f\n", z_orig, tmp2[ii]);
                                            exit(0);
                                        }
                                    }
                                }

                                //int x1 = int(x_orig);
                                __m512i x1 = _mm512_cvt_roundps_epi32(x_ori, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
//                                if (cnt <= 10) {
//                                    int tmp[16];
//                                    _mm512_store_epi32(tmp, x1);
//                                    cout << "x1 cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "x1 std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%d ", int((i + ii - x_orig_offset) * theta_rad_cos - A));
//                                    cout << endl;
//                                }

                                //int x2 = int(x_orig + 1 - eps);
                                __m512i x2 = _mm512_cvt_roundps_epi32(_mm512_add_ps(x_ori, eps_con),
                                                                      _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

//                                if (cnt <= 10) {
//                                    int tmp[16];
//                                    _mm512_store_epi32(tmp, x2);
//                                    cout << "x2 cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "x2 std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%d ", int((i + ii - x_orig_offset) * theta_rad_cos - A + 1 - eps));
//                                    cout << endl;
//                                }

                                //float coeff = x_orig - x1;
                                __m512 coeff = _mm512_sub_ps(x_ori, _mm512_cvtepi32_ps(x1));

                                //int n_z = int((z_orig + C) / defocus_step);
                                //TODO change div to mul
                                __m512 nz_tmp = _mm512_mul_ps(_mm512_add_ps(z_ori, c_con), d_step_con);

//                                if (cnt <= 10) {
//                                    float tmp[16];
//                                    _mm512_store_ps(tmp, nz_tmp);
//                                    cout << "nz_tmp cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "nz_tmp std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%.3f ",
//                                               (((i + ii - x_orig_offset) * theta_rad_sin + B + C) / defocus_step));
//                                    cout << endl;
//                                }
                                __m512i n_z = _mm512_cvt_roundps_epi32(nz_tmp,
                                                                       _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);


//                                if (cnt <= 10) {
//                                    int tmp[16];
//                                    _mm512_store_epi32(tmp, n_z);
//                                    cout << "n_z cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "n_z std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%d ",
//                                               int(((i + ii - x_orig_offset) * theta_rad_sin + B + C) / defocus_step));
//                                    cout << endl;
//                                }

                                //cal  n_z * Nx2 * Ny + j * Nx2
                                __m512i tt0 = _mm512_mullo_epi32(n_z, nxy_con);
                                __m512i tt1 = _mm512_mullo_epi32(j_con, nx2_con);
//                                if (cnt <= 10) {
//                                    int tmp[16];
//                                    _mm512_store_epi32(tmp, n_z);
//                                    cout << "n_z cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//                                    _mm512_store_epi32(tmp, nxy_con);
//                                    cout << "nxy_con cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//
//                                    _mm512_store_epi32(tmp, tt0);
//                                    cout << "tt0 cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//                                }


                                __m512i ids_base = _mm512_add_epi32(tt0, tt1);

//                                if (cnt <= 10) {
//                                    int tmp[16];
//                                    _mm512_store_epi32(tmp, ids_base);
//                                    cout << "ids_base cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//                                    int nz_id[16];
//                                    _mm512_store_epi32(nz_id, n_z);
//
//                                    cout << "ids_base std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%d ", nz_id[ii] * Nx2 * Ny + j * Nx2);
//                                    cout << endl;
//                                }

                                //cal  n_z * Nx2 * Ny + j * Nx2 + x1
                                __m512i ids1 = _mm512_add_epi32(ids_base, x1);
                                //cal  n_z * Nx2 * Ny + j * Nx2 + x2
                                __m512i ids2 = _mm512_add_epi32(ids_base, x2);

//                                if (cnt <= 10) {
//                                    int tmp[16];
//                                    _mm512_store_epi32(tmp, ids1);
//                                    cout << "ids1 cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%d ", tmp[ii]);
//                                    cout << endl;
//                                    int base_id[16];
//                                    _mm512_store_epi32(base_id, ids_base);
//                                    int x1_id[16];
//                                    _mm512_store_epi32(x1_id, x1);
//                                    cout << "ids1 std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)
//                                        printf("%d ", base_id[ii] + x1_id[ii]);
//                                    cout << endl;
//                                }


                                __m512 cc0 = _mm512_i32gather_ps(ids1, stack_corrected, 4);

                                __m512 cc1 = _mm512_i32gather_ps(ids2, stack_corrected, 4);

//                                cout << "111" << endl;

                                //(1 - coeff) * stack_corrected[n_z][j * Nx2 + x1]
                                __m512 c0 = _mm512_mul_ps(_mm512_sub_ps(one_con, coeff), cc0);


                                //(coeff) * stack_corrected[n_z][j * Nx2 + x2]
                                __m512 c1 = _mm512_mul_ps(coeff, cc1);


                                if (checcc) {
                                    float tmp1[16];
                                    _mm512_store_ps(tmp1, c0);
                                    float tmp2[16];
                                    _mm512_store_ps(tmp2, c1);
                                    for (int ii = 0; ii < 16; ii++) {
                                        float x_orig = (i + ii - x_orig_offset) * theta_rad_cos - A;
                                        float z_orig = (i + ii - x_orig_offset) * theta_rad_sin + B;
                                        int x1 = int(x_orig);
                                        int x2 = int(x_orig + 1 - eps);
                                        float coeff = x_orig - x1;
                                        int n_z = int((z_orig + C) / defocus_step);
                                        // the num in the corrected stack for the current height

                                        float ttr0 = (1 - coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x1];
                                        float ttr1 = (coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x2];

                                        if (fabs(ttr0 - tmp1[ii]) > 1e-5) {
                                            printf("GG\n");
                                            printf("c0 on %d %d %d\n", j, k, i + ii);
                                            printf("%.6f %.6f\n", ttr0, tmp1[ii]);
                                            exit(0);
                                        }
                                        if (fabs(ttr1 - tmp2[ii]) > 1e-5) {
                                            printf("GG\n");
                                            printf("c1 on %d %d %d\n", j, k, i + ii);
                                            printf("%.6f %.6f\n", ttr1, tmp2[ii]);
                                            exit(0);
                                        }


                                    }
                                }

//                                if (cnt <= 10) {
//                                    float tmp[16];
//                                    _mm512_store_ps(tmp, c0);
//                                    cout << "c0 cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "c0 std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++) {
//                                        int idd = i + ii;
//                                        float x_orig = (idd - x_orig_offset) * theta_rad_cos - A;
//                                        float z_orig = (idd - x_orig_offset) * theta_rad_sin + B;
//                                        int x1 = int(x_orig);
//                                        int x2 = int(x_orig + 1 - eps);
//                                        float coeff = x_orig - x1;
//                                        int n_z = int((z_orig + C) / defocus_step);
//                                        float tttx = (1 - coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x1];
//                                        printf("%.3f ", tttx);
//                                    }
//                                    cout << endl;
//                                    _mm512_store_ps(tmp, c1);
//                                    cout << "c1 cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "c1 std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++) {
//                                        int idd = i + ii;
//                                        float x_orig = (idd - x_orig_offset) * theta_rad_cos - A;
//                                        float z_orig = (idd - x_orig_offset) * theta_rad_sin + B;
//                                        int x1 = int(x_orig);
//                                        int x2 = int(x_orig + 1 - eps);
//                                        float coeff = x_orig - x1;
//                                        int n_z = int((z_orig + C) / defocus_step);
//                                        float tttx = (coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x2];
//                                        printf("%.3f ", tttx);
//                                    }
//                                    cout << endl;
//                                }

//                                cout << "222" << endl;


                                float *p_now = stack_recon[j] + k * Nx + i;
                                __m512 tmplod = _mm512_load_ps(p_now);
//                                if (cnt <= 10) {
//                                    float tmp[16];
//                                    _mm512_store_ps(tmp, tmplod);
//                                    cout << "tmplod cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "tmplod std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++) {
//                                        int idd = i + ii;
//                                        printf("%.3f ", stack_recon[j][idd + k * Nx]);
//                                    }
//                                    cout << endl;
//                                }
                                //last baba
                                __m512 resss = _mm512_add_ps(_mm512_add_ps(c0, c1), tmplod);
                                _mm512_store_ps(p_now, resss);

//                                if (cnt <= 10) {
//                                    __m512 rechr = _mm512_load_ps(p_now);
//                                    float tmp[16];
//                                    _mm512_store_ps(tmp, rechr);
//                                    cout << "rechr cal:" << endl;
//                                    for (int ii = 0; ii < 16; ii++)printf("%.3f ", tmp[ii]);
//                                    cout << endl;
//                                    cout << "c0 std:" << endl;
//                                    for (int ii = 0; ii < 16; ii++) {
//                                        int idd = i + ii;
//                                        printf("%.3f ", 0);
//                                    }
//                                    cout << endl;
//                                }

                                idx = _mm512_add_epi32(idx, si_con);

//                                exit(0);
//                                cout << "333" << endl;

                                // the num in the corrected stack for the current height
//                                stack_recon[j][i + k * Nx] +=
//                                        (1 - coeff) * stack_corrected[n_z][j * Nx2 + x1] +
//                                        (coeff) * stack_corrected[n_z][j * Nx2 + x2];
                            }
                            for (; i <= ri; i++)   // loop for the xz-plane to perform BP
                            {
                                float x_orig = (i - x_orig_offset) * theta_rad_cos - A;
                                float z_orig = (i - x_orig_offset) * theta_rad_sin + B;
                                int x1 = int(x_orig);
                                int x2 = int(x_orig + 1 - eps);
                                float coeff = x_orig - x1;
                                int n_z = int((z_orig + C) / defocus_step);
                                // the num in the corrected stack for the current height
                                stack_recon[j][i + k * Nx] +=
                                        (1 - coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x1] +
                                        (coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x2];
                            }

                        }
                    }
                    t3 = GetTime();
                    printf("reconstruction cost %.3f\n", t3 - t2);
                    cost2 += t3 - t2;
                    cout << "\tDone" << endl;

                    t2 = GetTime();
                    delete stack_corrected;

                    t3 = GetTime();
                    cost3 += t3 - t2;
                }
            }
            delete[] image_now;
            delete[] image_now_backup;
            printf("image %d cost %.3f\n", n, GetTime() - t1);
        }

        printf("main for cost %.3f\n", GetTime() - t0);
        printf("init and read cost %.3f\n", cost);
        printf("weight cost %.3f\n", cost0);
        printf("3DF cost %.3f\n", cost1);
        printf("reconstruction cost %.3f\n", cost2);
        printf("free cost %.3f\n", cost3);

        // write out final result
        cout << "Wrtie out final reconstruction result:" << endl;
        MRC stack_final(output_mrc.c_str(), "wb");
        stack_final.createMRC_empty(stack_orig.getNx(), h, stack_orig.getNy(), 2); // (x,z,y)
        // loop: Ny (number of xz-slices)
        for (int j = 0; j < stack_orig.getNy(); j++) {
            stack_final.write2DIm(stack_recon[j], j);
        }

        // update MRC header
        int threads = 3;
        float min_thread[threads], max_thread[threads];
        double mean_thread[threads];
        for (int th = 0; th < threads; th++) {
            min_thread[th] = stack_recon[0][0];
            max_thread[th] = stack_recon[0][0];
            mean_thread[th] = 0.0;
        }

        for (int j = 0; j < stack_orig.getNy(); j++) {
            double mean_now = 0.0;
            for (int i = 0; i < stack_orig.getNx() * h; i++) {
                mean_now += stack_recon[j][i];
                if (min_thread[j % threads] > stack_recon[j][i]) {
                    min_thread[j % threads] = stack_recon[j][i];
                }
                if (max_thread[j % threads] < stack_recon[j][i]) {
                    max_thread[j % threads] = stack_recon[j][i];
                }
            }
            mean_thread[j % threads] += (mean_now / (stack_orig.getNx() * h));
        }
        float min_all = min_thread[0];
        float max_all = max_thread[0];
        double mean_all = 0;
        for (int th = 0; th < threads; th++) {
            mean_all += mean_thread[th];
            if (min_all > min_thread[th]) {
                min_all = min_thread[th];
            }
            if (max_all < max_thread[th]) {
                max_all = max_thread[th];
            }
        }
        mean_all /= stack_orig.getNy();
        stack_final.computeHeader(pix, false, min_all, max_all, float(mean_all));

        for (int j = 0; j < stack_orig.getNy(); j++) {
            delete[] stack_recon[j];
        }

        stack_final.close();
        cout << "Done" << endl;

    }

    stack_orig.close();

    cout << endl << "Finish reconstruction successfully!" << endl;
    cout << "All results save in: " << path << endl << endl;

}