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
#include <stdio.h> /* for perror() */
#include <unistd.h> /* for syscall() */
#include <sys/syscall.h> /* for __NR_* definitions */
#include <linux/aio_abi.h> /* for AIO types and constants */
#include <fcntl.h> /* O_RDWR */
#include <string.h> /* memset() */
#include <inttypes.h> /* uint64_t */
#include "omp.h"

#include <immintrin.h>


const int threadNumber = 64;

const float eps = 1e-7;


double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

inline int io_setup(unsigned nr, aio_context_t *ctxp) {
    return syscall(__NR_io_setup, nr, ctxp);
}

inline int io_destroy(aio_context_t ctx) {
    return syscall(__NR_io_destroy, ctx);
}

inline int io_submit(aio_context_t ctx, long nr, struct iocb **iocbpp) {
    return syscall(__NR_io_submit, ctx, nr, iocbpp);
}

inline int
io_getevents(aio_context_t ctx, long min_nr, long max_nr, struct io_event *events, struct timespec *timeout) {
    return syscall(__NR_io_getevents, ctx, min_nr, max_nr, events, timeout);
}

static void buf2fft(float *buf, float *fft, int nx, int ny) {
    int nxb = nx + 2 - nx % 2;
    int i;
#pragma omp parallel for num_threads(threadNumber)
    for (i = 0; i < (nx + 2 - nx % 2) * ny; i++) {
        fft[i] = 0.0;
    }
#pragma omp parallel for num_threads(threadNumber)
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

void ReconstructionAlgo_WBP_RAM::doReconstruction(map <string, string> &inputPara, map <string, string> &outputPara) {

    double t_start = GetTime();
    double tt_t = GetTime();

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

    printf("para cost %.3f\n", GetTime() - tt_t);


    // Reconstruction
    cout << endl << "Reconstruction with (W)BP in RAM:" << endl << endl;
    if (!unrotated_stack)    // input rotated stack (y-axis as tilt axis)
    {
        double cost0 = 0;
        double cost1 = 0;
        double cost2 = 0;
        double cost2_5 = 0;
        double cost3 = 0;
        double cost4 = 0;
        double cost5 = 0;
        double cost6 = 0;
        double cost7 = 0;
        double cost8 = 0;


        double t0, t1, t2, t3;
        t0 = GetTime();

        cout << "Using rotated stack" << endl;

//        float *stack_recon[stack_orig.getNy()]; // (x,z,y)
//#pragma omp parallel for num_threads(threadNumber)
//        for (int j = 0; j < stack_orig.getNy(); j++) {
//            stack_recon[j] = new float[stack_orig.getNx() * h];
//            for (int i = 0; i < stack_orig.getNx() * h; i++) {
//                stack_recon[j][i] = 0.0;
//            }
//        }
        float *stack_recon = new float[stack_orig.getNy() * stack_orig.getNx() * h]; // (x,z,y)
#pragma omp parallel for num_threads(threadNumber)
        for (int j = 0; j < stack_orig.getNy(); j++) {
            for (int i = 0; i < stack_orig.getNx() * h; i++) {
                stack_recon[j * stack_orig.getNx() * h + i] = 0.0;
            }
        }

        cout << "Start reconstruction:" << endl;
        float x_orig_offset = float(stack_orig.getNx()) / 2.0, z_orig_offset = float(h) / 2.0;
        // loop: Nz (number of images)




        cout << "range " << stack_orig.getNx() << " " << stack_orig.getNy() << " " << stack_orig.getNz() << endl;

        fftwf_init_threads();
        fftwf_plan_with_nthreads(1);
        int Nx = stack_orig.getNx();
        int Ny = stack_orig.getNy();
        int Nx2 = Nx + 2 - Nx % 2;
        int Nxh = Nx / 2 + 1;
        int Nyh = Ny / 2 + 1;
        int Nyh2 = (Ny + 1) / 2;
        int Nzz = (int(h_tilt_max / defocus_step) + 1);
        float *stack_corrected = new float[Nzz * Nx2 * Ny];

        float *atan_xy = new float[Ny * Nx];
        float *sin_atan_xy = new float[Ny * Nx];
        float *cos_atan_xy = new float[Ny * Nx];
#pragma omp parallel for num_threads(threadNumber)
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                float x_norm = i;
                float y_norm = (j >= Nyh) ? (j - Ny) : (j);
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
                atan_xy[j * Nx + i] = alpha;
                sin_atan_xy[j * Nx + i] = sin(2 * alpha);
                cos_atan_xy[j * Nx + i] = cos(2 * alpha);
            }
        }
        int *mak_pre = new int[1 << 8];
#pragma omp parallel for num_threads(threadNumber)
        for (int i = 0; i < (1 << 8); i++) {
            mak_pre[i] = 0;
            int now = 0;
            for (int j = 7; j >= 0; j--) {
                if ((i >> j) & 1) {
                    now <<= 1;
                    now++;
                    now <<= 1;
                    now++;
                } else {
                    now <<= 2;
                }
            }
            mak_pre[i] = now;
        }


        for (int n = 0; n < stack_orig.getNz(); n++) {
            t1 = GetTime();
            t2 = GetTime();
            cout << "Image " << n << ":" << endl;
            float theta_rad = theta[n] / 180 * M_PI;
            double theta_rad_cos = cos(theta_rad);
            double theta_rad_sin = sin(theta_rad);
            float *image_now = new float[stack_orig.getNx() * stack_orig.getNy()];
            stack_orig.read2DIm_32bit(image_now, n);
            t3 = GetTime();
            cost0 += t3 - t2;

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
                        stack_recon[j * Ny * h + i] += recon_now[i];
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
                            stack_recon[j * Nx * h + i] += recon_now[i];
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
                    int n_zz = 0;
                    fftwf_plan_with_nthreads(64);

                    // weighting
                    if (!skip_weighting) {
                        cout << "\tStart weighting..." << endl;

                        filter_weighting_1D_many(image_now, stack_orig.getNx(), stack_orig.getNy(), weighting_radial,
                                                 weighting_sigma);
                        cout << "\tDone" << endl;
                    }
                    t3 = GetTime();
                    cost1 += t3 - t2;

                    // 3D-CTF correction
                    //hotspot 1
                    cout << "\tStart 3D-CTF correction..." << endl;

                    //把第一波fft操作拿出来，预处理好放到bufc_pre中
                    float *bufc_pre = new float[(Nx + 2 - Nx % 2) * Ny];
                    t2 = GetTime();

                    fftwf_plan plan_fft_pre = fftwf_plan_dft_r2c_2d(Ny, Nx, (float *) bufc_pre,
                                                                    reinterpret_cast<fftwf_complex *>(bufc_pre),
                                                                    FFTW_ESTIMATE);
                    buf2fft(image_now, bufc_pre, Nx, Ny);
                    fftwf_execute(plan_fft_pre);
                    fftwf_destroy_plan(plan_fft_pre);
                    fftwf_plan_with_nthreads(1);

                    t3 = GetTime();
                    cost2 += t3 - t2;
                    t2 = GetTime();
                    int n_zz_thread[threadNumber];
                    for (int i = 0; i < threadNumber; i++) {
                        n_zz_thread[i] = 0;
                    }
                    int zl = -int(h_tilt_max / 2);
                    int zr = int(h_tilt_max / 2);


                    t3 = GetTime();
                    cost2_5 += t3 - t2;
                    t2 = GetTime();

                    int nz_range = int(h_tilt_max / defocus_step) + 1;

                    t3 = GetTime();
                    cost3 += t3 - t2;
                    t2 = GetTime();
                    CTF ctf = ctf_para[n];

                    float defocus1 = ctf.getDefocus1();
                    float defocus2 = ctf.getDefocus2();
                    float astig = ctf.getAstigmatism();
                    float lambda = ctf.getLambda();
                    float phase_shift = ctf.getPhaseShift();
                    float w_cos = ctf.getW();
                    float w_sin = sqrt(1 - w_cos * w_cos);
                    float pix = ctf.getPixelSize();
                    float Cs = ctf.getCs();


                    float *atan_xy_cal = new float[Nx * Ny];
#pragma omp parallel for num_threads(threadNumber)
                    for (int j = 0; j < Ny; j++) {
                        for (int i = 0; i < Nx; i++) {
                            float alpha = atan_xy[j * Nx + i];
                            atan_xy_cal[j * Nx + i] = (defocus1 - defocus2) * (cos(2 * (alpha - astig)));
                        }
                    }

                    __m512 zero_con = _mm512_set1_ps(0);
                    __m512d zero_cond = _mm512_set1_pd(0);
                    __m512 ofive = _mm512_set1_ps(0.5);
                    __m512 two_con = _mm512_set1_ps(2);
                    __m512 neg_one = _mm512_set1_ps(-1);
                    __m512i si_con = _mm512_set1_epi32(16);
                    __m512i xor_neg = _mm512_set1_epi32(0x80000000);
                    __m512 Nx_con = _mm512_set1_ps(Nx);
                    __m512 div_pix = _mm512_set1_ps(1 / pix);
//
#pragma omp parallel for num_threads(threadNumber)
                    for (int zz = zl; zz < zr; zz += defocus_step) {

                        int thread_id = omp_get_thread_num();
                        int n_z = (zz + int(h_tilt_max / 2)) / defocus_step;
                        float *image = stack_corrected + n_z * Nx2 * Ny;
                        float z_offset = float(zz) + float(defocus_step - 1) / 2;
                        memcpy(image, bufc_pre, sizeof(float) * (Nx2 * Ny));


                        float A = (defocus1 + defocus2 - 2 * z_offset * pix);
                        float sin2ast = sin(2 * astig);
                        float cos2ast = cos(2 * astig);
                        double a_w_cos = acos(w_cos);

                        float div_Nx = 1.0 / float(Nx) * (1 / pix);

                        __m512 a_con = _mm512_set1_ps(A);
                        __m512 d1_d2 = _mm512_set1_ps(defocus1 - defocus2);
                        __m512 div_con = _mm512_set1_ps(div_Nx);

                        __m512 astig_con = _mm512_set1_ps(astig);
                        __m512 lambda_con = _mm512_set1_ps(lambda);
                        __m512d M_PI_2_con = _mm512_set1_pd(M_PI_2);
                        __m512d M_PI_23_con = _mm512_set1_pd(M_PI_2 * 3);
                        __m512d M_PI_22con = _mm512_set1_pd(2 * M_PI);
                        __m512d M_PI_22con_div = _mm512_set1_pd(1.0 / 2 / M_PI);
                        __m512d a_w_cos_con = _mm512_set1_pd(a_w_cos);


                        //TODO
                        __m512 mul_tmp1 = _mm512_set1_ps(M_PI * lambda);
                        __m512 mul_tmp2 = _mm512_set1_ps(M_PI_2 * Cs * lambda * lambda * lambda);
                        __m512d mul_tmp1d = _mm512_set1_pd(M_PI * lambda);
                        __m512d mul_tmp2d = _mm512_set1_pd(M_PI_2 * Cs * lambda * lambda * lambda);


                        __m512d phase_shift_cond = _mm512_set1_pd(phase_shift);
                        __m512 phase_shift_con = _mm512_set1_ps(phase_shift);

                        __m512 sin2ast_con = _mm512_set1_ps(sin2ast);
                        __m512 cos2ast_con = _mm512_set1_ps(cos2ast);


                        for (int j = 0; j < Ny; j++) {

                            float y_norm = (j >= Nyh) ? (j - Ny) : (j);
                            float y_real = float(y_norm) / float(Ny) * (1 / pix);


                            __m512i idx = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
                            __m512i jNx_con = _mm512_set1_epi32(j * Nx);
                            __m512i jNx2_con = _mm512_set1_epi32(j * Nx2);
                            __m512 y_real_con = _mm512_set1_ps(y_real);

                            int i = 0;
                            for (; i + 16 <= Nxh; i += 16) {
                                //float x_norm = i;
                                __m512 x_norm = _mm512_cvtepi32_ps(idx);

                                //float x_real = x_norm / float(Nx) * (1 / pix);
                                __m512 x_real = _mm512_mul_ps(_mm512_div_ps(x_norm, Nx_con), div_pix);

                                //float alpha = atan_xy[j * Nx + i];
                                __m512 alpha = _mm512_load_ps(atan_xy + j * Nx + i);

//
                                //float freq2 = x_real * x_real + y_real * y_real;
                                __m512 freq2 = _mm512_add_ps(_mm512_mul_ps(x_real, x_real),
                                                             _mm512_mul_ps(y_real_con, y_real_con));

//                                float df_now = (A + (defocus1 - defocus2) * mycos(2 * (alpha - astig))) / 2.0;
//                                float df_now = (A + (defocus1 - defocus2) * (cos_atan_xy[j * Nx + i] * cos2ast +
//                                                                             sin_atan_xy[j * Nx + i] * sin2ast)) / 2.0;

//                                __m512 cos_tmp = _mm512_cos_ps(
//                                        _mm512_mul_ps(two_con, _mm512_sub_ps(alpha, astig_con)));
//                                __m512 cos_tmp = zero_con;

//                                cos_atan_xy[j * Nx + i] * cos2ast + sin_atan_xy[j * Nx + i] * sin2ast
//                                __m512 cos_tmp = _mm512_add_ps(
//                                        _mm512_mul_ps(_mm512_load_ps(cos_atan_xy + j * Nx + i), cos2ast_con),
//                                        _mm512_mul_ps(_mm512_load_ps(sin_atan_xy + j * Nx + i), sin2ast_con));

                                //float df_now = (A + atan_xy_cal[j * Nx + i]) / 2.0;
                                __m512 df_now = _mm512_mul_ps(
                                        _mm512_add_ps(a_con, _mm512_load_ps(atan_xy_cal + j * Nx + i)), ofive);

//                                __m512 df_now = _mm512_mul_ps(_mm512_add_ps(a_con, _mm512_mul_ps(d1_d2, cos_tmp)),
//                                ofive);
                                __m256 f0 = _mm512_extractf32x8_ps(df_now, 0);
                                __m256 f1 = _mm512_extractf32x8_ps(df_now, 1);
                                __m512d df_now1 = _mm512_cvtps_pd(f0);
                                __m512d df_now2 = _mm512_cvtps_pd(f1);
                                f0 = _mm512_extractf32x8_ps(freq2, 0);
                                f1 = _mm512_extractf32x8_ps(freq2, 1);
                                __m512d freq21 = _mm512_cvtps_pd(f0);
                                __m512d freq22 = _mm512_cvtps_pd(f1);
                                __m512d tmp11 = _mm512_mul_pd(mul_tmp1d, _mm512_mul_pd(df_now1, freq21));
                                __m512d tmp12 = _mm512_mul_pd(mul_tmp1d, _mm512_mul_pd(df_now2, freq22));

                                __m512d tmp21 = _mm512_mul_pd(mul_tmp2d, _mm512_mul_pd(freq21, freq21));
                                __m512d tmp22 = _mm512_mul_pd(mul_tmp2d, _mm512_mul_pd(freq22, freq22));
                                __m512d chi1d = _mm512_add_pd(_mm512_sub_pd(tmp11, tmp21), phase_shift_cond);
                                __m512d chi2d = _mm512_add_pd(_mm512_sub_pd(tmp12, tmp22), phase_shift_cond);

                                __m512d chi1 = _mm512_cvtps_pd(_mm512_cvtpd_ps(chi1d));
                                __m512d chi2 = _mm512_cvtps_pd(_mm512_cvtpd_ps(chi2d));


//                                __m512 df_now = _mm512_mul_ps(a_con, ofive);

//                                    float chi = M_PI * lambda * df_now * freq2 -
//                                                M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2
//                                                + phase_shift;
                                //M_PI * lambda * df_now * freq2
//                                __m512 tmp1 = _mm512_mul_ps(mul_tmp1, _mm512_mul_ps(df_now, freq2));

                                //M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2
//                                __m512 tmp2 = _mm512_mul_ps(mul_tmp2, _mm512_mul_ps(freq2, freq2));

                                //cal

//                                __m512 chi = _mm512_add_ps(_mm512_sub_ps(tmp1, tmp2), phase_shift_con);

                                //double rrr = chi - a_w_cos;

//                                __m256 chi0 = _mm512_extractf32x8_ps(chi, 0);
//                                __m256 chi1 = _mm512_extractf32x8_ps(chi, 1);
                                __m512d rr0 = _mm512_abs_pd(_mm512_sub_pd(chi1, a_w_cos_con));
                                __m512d rr1 = _mm512_abs_pd(_mm512_sub_pd(chi2, a_w_cos_con));

                                //int mul_base = 1;
//                                rrr -= int(rrr / (2 * M_PI)) * 2 * M_PI;
//                                if (rrr > M_PI_2 && rrr < M_PI_2 * 3)mul_base = -1;


//
                                rr0 = _mm512_sub_pd(rr0, _mm512_mul_pd(_mm512_cvt_roundepi64_pd(
                                        _mm512_cvt_roundpd_epi64(
                                                _mm512_mul_pd(rr0,
                                                              M_PI_22con_div),
                                                _MM_FROUND_TO_NEG_INF |
                                                _MM_FROUND_NO_EXC),
                                        _MM_FROUND_TO_NEG_INF |
                                        _MM_FROUND_NO_EXC), M_PI_22con));

                                rr1 = _mm512_sub_pd(rr1, _mm512_mul_pd(_mm512_cvt_roundepi64_pd(
                                        _mm512_cvt_roundpd_epi64(
                                                _mm512_mul_pd(rr1,
                                                              M_PI_22con_div),
                                                _MM_FROUND_TO_NEG_INF |
                                                _MM_FROUND_NO_EXC),
                                        _MM_FROUND_TO_NEG_INF |
                                        _MM_FROUND_NO_EXC), M_PI_22con));
                                __mmask8 mul_base0 = 0;
                                __mmask8 mul_base1 = 0;
//
                                mul_base0 = _mm512_cmp_pd_mask(rr0, M_PI_2_con, _CMP_GT_OS) &
                                            _mm512_cmp_pd_mask(rr0, M_PI_23_con, _CMP_LT_OS);
                                mul_base1 = _mm512_cmp_pd_mask(rr1, M_PI_2_con, _CMP_GT_OS) &
                                            _mm512_cmp_pd_mask(rr1, M_PI_23_con, _CMP_LT_OS);


//                                double tmpd1[8], tmpd2[8];
//                                _mm512_store_pd(tmpd1, rr0);
//                                _mm512_store_pd(tmpd2, rr1);
//                                for (int ii = 0; ii < 8; ii++) {
//                                    double rrr = tmpd1[ii];
//                                    int mul_base = 1;
//                                    rrr -= int(rrr / (2 * M_PI)) * 2 * M_PI;
//                                    if (rrr > M_PI_2 && rrr < M_PI_2 * 3)mul_base = -1;
//                                    if (mul_base < 0)mul_base0 |= (1 << ii);
//                                }
//
//                                for (int ii = 0; ii < 8; ii++) {
//                                    double rrr = tmpd2[ii];
//                                    int mul_base = 1;
//                                    rrr -= int(rrr / (2 * M_PI)) * 2 * M_PI;
//                                    if (rrr > M_PI_2 && rrr < M_PI_2 * 3)mul_base = -1;
//                                    if (mul_base < 0)mul_base1 |= (1 << ii);
//                                }
//                                mul_base0 = _mm512_cmp_pd_mask(_mm512_cos_pd(rr0), zero_cond, _CMP_LT_OS);
//                                mul_base1 = _mm512_cmp_pd_mask(_mm512_cos_pd(rr1), zero_cond, _CMP_LT_OS);



                                if (flip_contrast) {
                                    mul_base0 = !mul_base0;
                                    mul_base1 = !mul_base1;
                                }


//                                __mmask16 mul_base = mul_base0 + ((int(mul_base1)) << 8);

                                __mmask16 big_base0 = mak_pre[mul_base0];
                                __mmask16 big_base1 = mak_pre[mul_base1];

                                __m512 img1 = _mm512_load_ps(image + j * Nx2 + i * 2);
                                __m512 img2 = _mm512_load_ps(image + j * Nx2 + i * 2 + 16);
                                img1 = _mm512_mask_xor_ps(img1, big_base0, img1, __m512(xor_neg));
                                img2 = _mm512_mask_xor_ps(img2, big_base1, img2, __m512(xor_neg));
//                                img1 = _mm512_mask_sub_ps(img1, big_base0, zero_con, img1);
//                                img2 = _mm512_mask_sub_ps(img2, big_base1, zero_con, img2);
                                _mm512_store_ps(image + j * Nx2 + i * 2, img1);
                                _mm512_store_ps(image + j * Nx2 + i * 2 + 16, img2);

//                                for (int ii = 0; ii < 16; ii++) {
//                                    int idi = i * 2 + ii;
//                                    if ((big_base0 >> ii) & 1) {
//                                        image[idi + j * Nx2] *= -1;
//                                    }
//                                }
//                                for (int ii = 0; ii < 16; ii++) {
//                                    int idi = i * 2 + ii + 16;
//                                    if ((big_base1 >> ii) & 1) {
//                                        image[idi + j * Nx2] *= -1;
//                                    }
//                                }

//                                for (int ii = 0; ii < 16; ii++) {
//                                    int idi = i + ii;
//                                    if ((mul_base >> ii) & 1) {
//                                        image[idi * 2 + j * Nx2] *= -1;
//                                        image[(idi * 2 + 1) + j * Nx2] *= -1;
//                                    }
//                                }

                                idx = _mm512_add_epi32(idx, si_con);

                            }


                            for (; i < Nxh; i++) {
                                float x_norm = i;
                                float x_real = float(x_norm) / float(Nx) * (1 / pix);

                                float freq2 = x_real * x_real + y_real * y_real;
                                float df_now = (A + atan_xy_cal[j * Nx + i]) / 2.0;
//                                float df_now = (A + (defocus1 - defocus2) * (cos(2 * (alpha - astig)))) / 2.0;

//                                float df_now = (A + (defocus1 - defocus2) * mycos(2 * (alpha - astig))) / 2.0;
                                float chi = (M_PI) * lambda * df_now * freq2 -
                                            (M_PI_2) * Cs * lambda * lambda * lambda * freq2 * freq2
                                            + phase_shift;

                                double rrr = chi - a_w_cos;
                                int mul_base = 1;
                                rrr = fabs(rrr);
                                rrr -= int(rrr / (2 * M_PI)) * 2 * M_PI;
                                if (rrr > M_PI_2 && rrr < M_PI_2 * 3)mul_base = -1;
                                if (flip_contrast) {
                                    mul_base = -mul_base;
                                }
                                if (mul_base < 0) {
                                    image[i * 2 + j * Nx2] *= -1;
                                    image[(i * 2 + 1) + j * Nx2] *= -1;
                                }
                            }
                        }
                        fftwf_plan plan_ifft;
//#pragma omp critical
                        {
                            plan_ifft = fftwf_plan_dft_c2r_2d(Ny, Nx, reinterpret_cast<fftwf_complex *>(image),
                                                              (float *) image, FFTW_ESTIMATE);
                        }
                        fftwf_execute(plan_ifft);
//#pragma omp critical
                        {
                            fftwf_destroy_plan(plan_ifft);
                        }
                        for (int i = 0; i < Nx2 * Ny; i++)image[i] = image[i] / (Nx * Ny);
                        n_zz_thread[thread_id]++;
                    }
//
//#pragma omp parallel for num_threads(threadNumber)
//                    for (int zz = zl; zz < zr; zz += defocus_step) {
//
//                        int thread_id = omp_get_thread_num();
//
//
//                        int n_z = (zz + int(h_tilt_max / 2)) / defocus_step;
//                        float *image = stack_corrected + n_z * Nx2 * Ny;
//                        CTF ctf = ctf_para[n];
//                        float z_offset = float(zz) + float(defocus_step - 1) / 2;
//                        memcpy(image, bufc_pre, sizeof(float) * (Nx2 * Ny));
//
//                        float defocus1 = ctf.defocus1;
//                        float defocus2 = ctf.defocus2;
//                        float astig = ctf.astig;
//                        float lambda = ctf.lambda;
//                        float phase_shift = ctf.phase_shift;
//                        float w_sin = ctf.w_sin;
//                        float w_cos = ctf.w_cos;
//                        float pix = ctf.pix;
//                        float Cs = ctf.Cs;
//                        float A = (defocus1 + defocus2 - 2 * z_offset * pix);
//                        float sin2ast = sin(2 * astig);
//                        float cos2ast = cos(2 * astig);
//
//                        float div_Nx = 1.0 / float(Nx) * (1 / pix);
//
//                        for (int j = 0; j < Ny; j++) {
//
//                            float y_norm = (j >= Nyh) ? (j - Ny) : (j);
//                            float y_real = float(y_norm) / float(Ny) * (1 / pix);
//
//                            for (int i = 0; i < Nxh; i++) {
//
//                                float x_norm = i;
//                                float x_real = float(x_norm) * div_Nx;
//                                float alpha = atan_xy[j * Nx + i];
//
//                                float freq2 = x_real * x_real + y_real * y_real;
//                                float df_now = (A + (defocus1 - defocus2) * (cos_atan_xy[j * Nx + i] * cos2ast +
//                                                                             sin_atan_xy[j * Nx + i] * sin2ast)) / 2.0;
////                                float df_now = (A + (defocus1 - defocus2) * (cos(2 * (alpha - astig)))) / 2.0;
////
////                                float df_now = (A + (defocus1 - defocus2) * mycos(2 * (alpha - astig))) / 2.0;
//                                float chi = (M_PI) * lambda * df_now * freq2 -
//                                            (M_PI_2) * Cs * lambda * lambda * lambda * freq2 * freq2
//                                            + phase_shift;
//                                double rrr = chi - a_w_cos;
//                                int mul_base = 1;
//                                //TODO check if rrr will bigger than 2*pi
//                                rrr = fabs(rrr);
//                                rrr -= int(rrr / (2 * M_PI)) * 2 * M_PI;
//                                if (rrr > M_PI_2 && rrr < M_PI_2 * 3)mul_base = -1;
////                                if (cos(rrr) < 0)mul_base = -1;
//                                if (flip_contrast) {
//                                    mul_base = -mul_base;
//                                }
//                                if (mul_base < 0) {
//                                    image[i * 2 + j * Nx2] *= -1;
//                                    image[(i * 2 + 1) + j * Nx2] *= -1;
//                                }
//                            }
//                        }
//                        fftwf_plan plan_ifft;
//
////#pragma omp critical
//                        {
//                            plan_ifft = fftwf_plan_dft_c2r_2d(Ny, Nx, reinterpret_cast<fftwf_complex *>(image),
//                                                              (float *) image, FFTW_ESTIMATE);
//                        }
//                        fftwf_execute(plan_ifft);
////#pragma omp critical
//                        {
//                            fftwf_destroy_plan(plan_ifft);
//                        }
//                        for (int i = 0; i < Nx2 * Ny; i++)image[i] = image[i] / (Nx * Ny);
//                        n_zz_thread[thread_id]++;
//                    }


                    for (int i = 0; i < threadNumber; i++) {
                        n_zz += n_zz_thread[i];
                    }


                    t3 = GetTime();
                    cost4 += t3 - t2;
                    t2 = GetTime();
                    delete[]atan_xy_cal;

                    t3 = GetTime();
                    cost5 += t3 - t2;
                    t2 = GetTime();
                    delete[] bufc_pre;


                    t3 = GetTime();
                    cost6 += t3 - t2;
                    cout << "\tDone!" << endl;



                    // recontruction

                    cout << "\tPerform reconstruction:" << endl;
                    t2 = GetTime();

                    // loop: Ny (number of xz-slices)


#pragma omp parallel for num_threads(threadNumber)
                    for (int j = 0; j < Ny; j++) {
                        __m512i nxy_con = _mm512_set1_epi32(Nx2 * Ny);
                        __m512i nx2_con = _mm512_set1_epi32(Nx2);
                        __m512i j_con = _mm512_set1_epi32(j);
                        __m512i si_con = _mm512_set1_epi32(16);
                        __m512d cos_con = _mm512_set1_pd(theta_rad_cos);
                        __m512d sin_con = _mm512_set1_pd(theta_rad_sin);
                        __m512 offset_con = _mm512_set1_ps(x_orig_offset);
                        __m512 eps_con = _mm512_set1_ps(1 - eps);
                        __m512 one_con = _mm512_set1_ps(1);
                        __m512 d_step_con = _mm512_set1_ps(1.0 / defocus_step);

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

                            __m512 c_con = _mm512_set1_ps(C);
                            __m512i idx = _mm512_set_epi32(li + 15, li + 14, li + 13, li + 12, li + 11, li + 10,
                                                           li + 9,
                                                           li + 8, li + 7, li + 6, li + 5, li + 4, li + 3, li + 2,
                                                           li + 1, li + 0);


//                            static int cnt = 100;
//                            bool checcc = 0;
                            for (; i + 16 <= ri + 1; i += 16) {



                                //float x_orig = (i - x_orig_offset) * theta_rad_cos - A;
                                //float z_orig = (i - x_orig_offset) * theta_rad_sin + B;


                                //i - x_orig_offset
                                __m512 sub_tmp = _mm512_sub_ps(_mm512_cvtepi32_ps(idx), offset_con);



                                //select pre 8 float to f0, last 8 to f1
                                __m256 f0 = _mm512_extractf32x8_ps(sub_tmp, 0);
                                __m256 f1 = _mm512_extractf32x8_ps(sub_tmp, 1);
                                __m512d d0 = _mm512_cvtps_pd(f0);
                                __m512d d1 = _mm512_cvtps_pd(f1);
                                // sub_tmp * theta_rad_cos - A
                                __m512d x_tmp0 = _mm512_sub_pd(_mm512_mul_pd(d0, cos_con), a_con);
                                __m512d x_tmp1 = _mm512_sub_pd(_mm512_mul_pd(d1, cos_con), a_con);

                                //merge 8+8 double to 16 float
                                __m512 x_ori;
                                x_ori = _mm512_insertf32x8(x_ori, _mm512_cvtpd_ps(x_tmp0), 0);
                                x_ori = _mm512_insertf32x8(x_ori, _mm512_cvtpd_ps(x_tmp1), 1);



                                // sub_tmp * theta_rad_sin + B
                                __m512d z_tmp0 = _mm512_add_pd(_mm512_mul_pd(d0, sin_con), b_con);
                                __m512d z_tmp1 = _mm512_add_pd(_mm512_mul_pd(d1, sin_con), b_con);
                                //merge 8+8 double to 16 float
                                __m512 z_ori;
                                z_ori = _mm512_insertf32x8(z_ori, _mm512_cvtpd_ps(z_tmp0), 0);
                                z_ori = _mm512_insertf32x8(z_ori, _mm512_cvtpd_ps(z_tmp1), 1);

                                //int x1 = int(x_orig);
                                __m512i x1 = _mm512_cvt_roundps_epi32(x_ori,
                                                                      _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

                                //int x2 = int(x_orig + 1 - eps);
                                __m512i x2 = _mm512_cvt_roundps_epi32(_mm512_add_ps(x_ori, eps_con),
                                                                      _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);



                                //float coeff = x_orig - x1;
                                __m512 coeff = _mm512_sub_ps(x_ori, _mm512_cvtepi32_ps(x1));

                                //int n_z = int((z_orig + C) / defocus_step);
                                //TODO change div to mul
                                __m512 nz_tmp = _mm512_mul_ps(_mm512_add_ps(z_ori, c_con), d_step_con);


                                __m512i n_z = _mm512_cvt_roundps_epi32(nz_tmp,
                                                                       _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);




                                //cal  n_z * Nx2 * Ny + j * Nx2
                                __m512i tt0 = _mm512_mullo_epi32(n_z, nxy_con);
                                __m512i tt1 = _mm512_mullo_epi32(j_con, nx2_con);


                                __m512i ids_base = _mm512_add_epi32(tt0, tt1);


                                //cal  n_z * Nx2 * Ny + j * Nx2 + x1
                                __m512i ids1 = _mm512_add_epi32(ids_base, x1);
                                //cal  n_z * Nx2 * Ny + j * Nx2 + x2
                                __m512i ids2 = _mm512_add_epi32(ids_base, x2);


                                __m512 cc0 = _mm512_i32gather_ps(ids1, stack_corrected, 4);

                                __m512 cc1 = _mm512_i32gather_ps(ids2, stack_corrected, 4);

//                                cout << "111" << endl;

                                //(1 - coeff) * stack_corrected[n_z][j * Nx2 + x1]
                                __m512 c0 = _mm512_mul_ps(_mm512_sub_ps(one_con, coeff), cc0);


                                //(coeff) * stack_corrected[n_z][j * Nx2 + x2]
                                __m512 c1 = _mm512_mul_ps(coeff, cc1);


//                                float *p_now = stack_recon[j] + k * Nx + i;
                                float *p_now = stack_recon + j * Nx * h + k * Nx + i;
                                __m512 tmplod = _mm512_load_ps(p_now);

                                //last baba
                                __m512 resss = _mm512_add_ps(_mm512_add_ps(c0, c1), tmplod);
                                _mm512_store_ps(p_now, resss);


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
                                stack_recon[j * Nx * h + i + k * Nx] +=
                                        (1 - coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x1] +
                                        (coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x2];
                            }


                        }
                    }



//#pragma omp parallel for num_threads(threadNumber)
//                    for (int j = 0; j < Ny; j++) {
//
//                        // loop: Nx*h (whole xz-slice)
//                        for (int k = 0; k < h; k++) {
//                            double l = 0;
//                            double r = Nx - 1;
//
//                            l = max(l, ((k - z_orig_offset) * theta_rad_sin - x_orig_offset)
//                                       / theta_rad_cos + x_orig_offset);
//                            r = min(r, ((k - z_orig_offset) * theta_rad_sin - x_orig_offset + Nx - 1)
//                                       / theta_rad_cos + x_orig_offset);
//                            if (theta_rad_cos < 0)swap(l, r);
//
//                            double ll = 0;
//                            double rr = Nx - 1;
//                            ll = max(ll, ((z_orig_offset - k) * theta_rad_cos - int(h_tilt_max / 2))
//                                         / theta_rad_sin + x_orig_offset);
//
//                            rr = min(rr, ((z_orig_offset - k) * theta_rad_cos - int(h_tilt_max / 2) +
//                                          n_zz * defocus_step) / theta_rad_sin + x_orig_offset);
//                            if (theta_rad_sin < 0)swap(ll, rr);
//                            l = max(l, ll);
//                            r = min(r, rr);
//
//                            int li = ceil(l);
//                            int ri = floor(r);
//                            double A = (k - z_orig_offset) * theta_rad_sin - x_orig_offset;
//                            double B = (k - z_orig_offset) * theta_rad_cos + z_orig_offset;
//                            float C = -z_orig_offset + int(h_tilt_max / 2);
//
//
//                            for (int i = li; i <= ri; i++)   // loop for the xz-plane to perform BP
//                            {
//                                float x_orig = (i - x_orig_offset) * theta_rad_cos - A;
//                                float z_orig = (i - x_orig_offset) * theta_rad_sin + B;
//                                int x1 = int(x_orig);
//                                int x2 = int(x_orig + 1 - eps);
//                                float coeff = x_orig - x1;
//                                int n_z = int((z_orig + C) / defocus_step);
//                                // the num in the corrected stack for the current height
//                                stack_recon[j * Nx * h + i + k * Nx] +=
//                                        (1 - coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x1] +
//                                        (coeff) * stack_corrected[n_z * Nx2 * Ny + j * Nx2 + x2];
//                            }
//
//
//                        }
//                    }

                    t3 = GetTime();
                    printf("reconstruction cost %.3f\n", t3 - t2);
                    cost7 += t3 - t2;
                    cout << "\tDone" << endl;
                }
            }
            delete[] image_now;
            printf("image %d cost %.3f\n", n, GetTime() - t1);
        }

        t2 = GetTime();
        delete stack_corrected;
        t3 = GetTime();
        cost8 += t3 - t2;


        printf("main for cost %.3f\n", GetTime() - t0);
        printf("init and read cost %.3f\n", cost0);
        printf("weight cost %.3f\n", cost1);
        printf("pre bufc cost %.3f\n", cost2);
        printf("pre atant cost %.3f\n", cost2_5);
        printf("fftw malloc cost %.3f\n", cost3);
        printf("hotspots1 cost %.3f\n", cost4);
        printf("fftw free cost %.3f\n", cost5);
        printf("bufc free cost %.3f\n", cost6);
        printf("rebu cost %.3f\n", cost7);
        printf("stack_corrected free cost %.3f\n", cost8);


        double t_write = GetTime();
        double tw = GetTime();
        // write out final result
        cout << "Wrtie out final reconstruction result:" << endl;
//        string outdir = "/dev/shm/" + output_mrc;
//        MRC stack_final(outdir.c_str(), "wb");
        MRC stack_final(output_mrc.c_str(), "wb");

        stack_final.createMRC_empty(stack_orig.getNx(), h, stack_orig.getNy(), 2); // (x,z,y)
        // loop: Ny (number of xz-slices)
//        printf("outdir %s\n", outdir.c_str());

        printf("write t1 cost %.3f\n", GetTime() - tw);
        tw = GetTime();


        float min_thread, max_thread;
        double mean_thread;
        min_thread = stack_recon[0];
        max_thread = stack_recon[0];
        mean_thread = 0.0;

#pragma omp parallel for num_threads(threadNumber)
        for (int j = 0; j < stack_orig.getNy(); j++) {
            double mean_now = 0.0;
            float mis = stack_recon[j * Nx * h + 0];
            float mxs = stack_recon[j * Nx * h + 0];

            for (int i = 0; i < stack_orig.getNx() * h; i++) {
                mean_now += stack_recon[j * Nx * h + i];
                mis = min(mis, stack_recon[j * Nx * h + i]);
                mxs = max(mxs, stack_recon[j * Nx * h + i]);
            }
#pragma omp critical
            {
                min_thread = min(min_thread, mis);
                max_thread = max(max_thread, mxs);
                mean_thread += (mean_now / (stack_orig.getNx() * h));
            }

        }
        float min_all = min_thread;
        float max_all = max_thread;
        double mean_all = mean_thread;


        printf("write t4 cost %.3f\n", GetTime() - tw);
        tw = GetTime();
        mean_all /= stack_orig.getNy();
        stack_final.computeHeader(pix, false, min_all, max_all, float(mean_all));

        printf("write t5 cost %.3f\n", GetTime() - tw);
        tw = GetTime();
//        for (int j = 0; j < stack_orig.getNy(); j++) {
//            delete[] stack_recon[j];
//        }
        printf("offset %d\n", stack_final.getSymdatasize());
        size_t ImSize = (size_t) stack_final.getImSize() * Ny;
        size_t offset = 1024 + stack_final.getSymdatasize();

        printf("ImSize %lld\n", 1ll * stack_final.getImSize() * Ny);
        printf("ImSize %zu\n", ImSize);

//        if (fseek(stack_final.getFp(), offset, SEEK_SET) != 0) {
//            printf("GG\n");
//        }
//
//        int res = fwrite(stack_recon, 1, ImSize, stack_final.getFp());
//        printf("res %d\n", res);
//        delete[] stack_recon;

        printf("write t6 cost %.3f\n", GetTime() - tw);
        tw = GetTime();
        stack_final.close();

        printf("write close cost %.3f\n", GetTime() - tw);
        tw = GetTime();

        //optimize fclose with aio

        {
            double t00 = GetTime();

            int N = ImSize;
            aio_context_t ctx;
            struct iocb cb;
            struct iocb *cbs[1];
            struct io_event events[1];
            int ret;
            int fd;
            int i;

            printf("p1 cost %lf\n", GetTime() - t00);
            t00 = GetTime();

            fd = open(output_mrc.c_str(), O_RDWR | O_CREAT, S_IRWXU);
            if (fd < 0) {
                perror("open error");
            }
            ctx = 0;
            ret = io_setup(128, &ctx);
            printf("after io_setup ctx:%ld\n", ctx);
            if (ret < 0) {
                perror("io_setup error");
            } /* setup I/O control block */
            memset(&cb, 0, sizeof(cb));
            cb.aio_fildes = fd;
            cb.aio_lio_opcode = IOCB_CMD_PWRITE;/* command-specific options */
            cb.aio_buf = (uint64_t) stack_recon;
            cb.aio_offset = 1024;
            cb.aio_nbytes = N * sizeof(float);
            cbs[0] = &cb;
            ret = io_submit(ctx, 1, cbs);
            if (ret != 1) {
                if (ret < 0)
                    perror("io_submit error");
                else
                    fprintf(stderr, "could not sumbit IOs");
            } /* get the reply */
            printf("p2 cost %lf\n", GetTime() - t00);
            t00 = GetTime();

            ret = io_destroy(ctx);
            if (ret < 0) {
                perror("io_destroy error");
            }


            printf("p3 cost %lf\n", GetTime() - t00);
            t00 = GetTime();

            delete[] stack_recon;
            printf("p4 cost %lf\n", GetTime() - t00);
            t00 = GetTime();
        }


        printf("write t7 cost %.3f\n", GetTime() - tw);
        tw = GetTime();
        cout << "Done" << endl;


        printf("write part cost %.3f\n", GetTime() - t_write);

    }

    stack_orig.close();

    cout << endl << "Finish reconstruction successfully!" <<
         endl;
    cout << "All results save in: " << path << endl <<
         endl;

    printf("function cost %.3f\n", GetTime() - t_start);

}