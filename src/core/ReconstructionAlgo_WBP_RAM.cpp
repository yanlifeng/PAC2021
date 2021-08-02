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
#include <omp.h>
#include <chrono>

#include <immintrin.h>

typedef std::chrono::high_resolution_clock Clock;
#define TDEF(x_) chrono::high_resolution_clock::time_point x_##_t0, x_##_t1;
#define TSTART(x_) x_##_t0 = Clock::now();
#define TEND(x_) x_##_t1 = Clock::now();
#define TINT(x_) chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()
#define TPRINT(x_, str) cout << endl << "        " << str << "  " << chrono::duration_cast<chrono::microseconds>(x_##_t1 - x_##_t0).count()/1e6 << "s" << endl << endl;

static void buf2fft(float *buf, float *fft, int nx, int ny) {
    int nxb = nx + 2 - nx % 2;
    int i;
    for (i = 0; i < nxb * ny; i++) {
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
    for (i = 0; i < (nx_final + 2 - nx_final % 2) * ny_final; i++) {
        fft[i] = 0.0;
    }
#pragma omp parallel for
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
    for (i = 0; i < nx_orig * ny_orig; i++) {
        buf[i] = 0.0;
    }
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
#pragma omp parallel for
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
#pragma omp parallel for
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

static void ctf_correction(float *image,int Nx,int Ny,CTF ctf,bool flip_contrast,float z_offset)   // z_offset in pixels
{
    fftwf_plan plan_fft,plan_ifft;
    float *bufc=new float[(Nx+2-Nx%2)*Ny];
    plan_fft=fftwf_plan_dft_r2c_2d(Ny,Nx,(float*)bufc,reinterpret_cast<fftwf_complex*>(bufc),FFTW_ESTIMATE);
    plan_ifft=fftwf_plan_dft_c2r_2d(Ny,Nx,reinterpret_cast<fftwf_complex*>(bufc),(float*)bufc,FFTW_ESTIMATE);
    buf2fft(image,bufc,Nx,Ny);
    fftwf_execute(plan_fft);

    // loop: Ny (all Fourier components for y-axis)
    for(int j=0;j<Ny;j++)
    {
        // loop: Nx+2-Nx%2 (all Fourier components for x-axis)
        for(int i=0;i<(Nx+2-Nx%2);i+=2)
        {
            float ctf_now=ctf.computeCTF2D(i/2,j,Nx,Ny,true,flip_contrast,z_offset);

            bufc[i+j*(Nx+2-Nx%2)]*=ctf_now;
            bufc[(i+1)+j*(Nx+2-Nx%2)]*=ctf_now;
        }
    }

    fftwf_execute(plan_ifft);
    fft2buf(image,bufc,Nx,Ny);
    for(int i=0;i<Nx*Ny;i++)   // normalization
    {
        image[i]=image[i]/(Nx*Ny);
    }
    fftwf_destroy_plan(plan_fft);
    fftwf_destroy_plan(plan_ifft);
    delete [] bufc;
}


static void
ctf_correction_perbufc(int Nx, int Ny,float defocus1,float defocus2,float astig,float lambda,float phase_shift,float w_sin,float w_cos,float pix,float Cs, bool flip_contrast, float z_offset,float *bufc)   // z_offset in pixels
{

    // loop: Ny (all Fourier components for y-axis)

//#pragma omp parallel for
    int Nxc = Nx + 2 -(Nx&1);
    for (int j = 0; j < Ny; j++) {
        // loop: Nx+2-Nx%2 (all Fourier components for x-axis)
        float y=j;
        float y_norm = (y >= int(ceil(float(Ny + 1) / 2))) ? (y - Ny) : (y);
        float y_real = float(y_norm) / float(Ny) * (1 / pix);
        float y_real_2 = y_real*y_real;
        // 特殊处理 x_norm = 0
        {
            float ctf_now;
            float x_norm = 0;
            float x_real = 0;
            float alpha;
            if (y_norm > 0) {
                alpha = M_PI_2;
            } else if (y_norm < 0) {
                alpha = -M_PI_2;
            } else {
                alpha = 0.0;
            }
            float freq2 = y_real_2;
            float df_now =
                    ((defocus1 + defocus2 - 2 * z_offset * pix) +
                     (defocus1 - defocus2) * cos(2 * (alpha - astig))) / 2.0;
            float chi = M_PI * lambda * df_now * freq2 -
                        M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2 + phase_shift;
            ctf_now = w_sin * sin(chi) + w_cos * cos(chi);
            if (flip_contrast) {
                ctf_now = -ctf_now;
            }
            ctf_now = 1-(int)(((*((unsigned int*)&ctf_now)) >> 31) << 1);
            bufc[j * Nxc] *= ctf_now;
            bufc[1 + j * Nxc] *= ctf_now;
        }
        for (int i = 2; i < Nxc; i += 2) {
            //float ctf_now = ctf.computeCTF2D(i / 2, j, Nx, Ny, true, flip_contrast, z_offset);
            float ctf_now;
            {
                float x = i / 2;
                float x_norm = (x >= int(ceil(float(Nx + 1) / 2))) ? (x - Nx) : (x);
                float x_real = float(x_norm) / float(Nx) * (1 / pix);
                float alpha = atan(y_real / x_real);
                float freq2 = x_real * x_real + y_real_2;
                float df_now =
                        ((defocus1 + defocus2 - 2 * z_offset * pix) +
                         (defocus1 - defocus2) * cos(2 * (alpha - astig))) / 2.0;
                float chi = M_PI * lambda * df_now * freq2 -
                            M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2 + phase_shift;
                ctf_now = w_sin * sin(chi) + w_cos * cos(chi);
                if (flip_contrast) {
                    ctf_now = -ctf_now;
                }
                ctf_now = 1-(int)(((*((unsigned int*)&ctf_now)) >> 31) << 1);
            }
            bufc[i + j * Nxc] *= ctf_now;
            bufc[(i + 1) + j * Nxc] *= ctf_now;
        }
    }
    fftwf_plan plan_ifft;
#pragma omp critical
    {
        plan_ifft = fftwf_plan_dft_c2r_2d(Ny, Nx, reinterpret_cast<fftwf_complex *>(bufc), (float *) bufc,
                                          FFTW_ESTIMATE);
    }
    fftwf_execute(plan_ifft);
  //  fft2buf(image, bufc, Nx, Ny);
#pragma omp simd
    for (int i = 0; i < Nxc * Ny; i++)   // normalization
    {
        bufc[i] = bufc[i] / (Nx * Ny);
    }
#pragma omp critical
    {
        fftwf_destroy_plan(plan_ifft);
    }
    //delete[] bufc;
}

ReconstructionAlgo_WBP_RAM::~ReconstructionAlgo_WBP_RAM() {

}

void ReconstructionAlgo_WBP_RAM::doReconstruction(map<string, string> &inputPara, map<string, string> &outputPara) {
    TDEF(INPUT)
    TSTART(INPUT)
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
    TEND(INPUT)
    TPRINT(INPUT, "Input file and check time is")

    /*
     * stack_orig.getNx()  682
     * stack_orig.getNy()  960
     * stack_orig.getNz()  69
     */
    TDEF(Reconstruction_ALL)
    TSTART(Reconstruction_ALL)
    int Nx = stack_orig.getNx();
    int Ny = stack_orig.getNy();
    int Nz = stack_orig.getNz();
    int Nxc = Nx + 2 - Nx % 2;
    int filter_weighting_1D_many_all = 0 ;
    int ctf_correction_all = 0;
    int recon_now_all = 0;
    // Reconstruction
    cout << endl << "Reconstruction with (W)BP in RAM:" << endl << endl;
    if (!unrotated_stack)    // input rotated stack (y-axis as tilt axis)
    {
        cout << "Using rotated stack" << endl;

        float *stack_recon[Ny]; // (x,z,y)
        for (int j = 0; j < Ny; j++) {
            stack_recon[j] = new float[Nx * h];
#pragma omp parallel for
            for (int i = 0; i < Nx * h; i++) {
                stack_recon[j][i] = 0.0;
            }
        }
        /*
         *  pre_atan
         */

        float *atan_pre_2 = new float[(Nxc/2)*Ny];
#pragma omp parallel for
        for (int i=0;i<Ny;i++){
            int y = i;
            float y_norm = (y >= int(ceil(float(Ny + 1) / 2))) ? (y - Ny) : (y);
#pragma omp simd
            for (int j=0;j<Nxc/2;j++){
                atan_pre_2[j+i*(Nxc/2)] = 2 * atan(y_norm / j * Nx / Ny);
            }
        }



        cout << "Start reconstruction:" << endl;
        float x_orig_offset = float(stack_orig.getNx()) / 2.0, z_orig_offset = float(h) / 2.0;
        // loop: Nz (number of images)
        TDEF(n_stack_orig_getNz)
        TSTART(n_stack_orig_getNz)
        float *stack_corrected[int(h_tilt_max / defocus_step) + 1]; // 第一维遍历不同高度，第二维x，第三维y
        for (int i=0;i<(int(h_tilt_max / defocus_step) + 1);i++) stack_corrected[i] = new float[Nxc*Ny];

        fftwf_init_threads();

        for (int n = 0; n < Nz; n++)   // loop for every micrograph
        {
            TDEF(Image)
            TSTART(Image)
            cout << "Image " << n << ":" << endl;
            float theta_rad = theta[n] / 180 * M_PI;
            double cos_theta_rad = cos(theta_rad);
            double sin_theta_rad = sin(theta_rad);
            float *image_now = new float[Nx * Ny];
            stack_orig.read2DIm_32bit(image_now, n);
            float *image_now_backup = new float[Nx * Ny];
            memcpy(image_now_backup, image_now, sizeof(float) * Nx * Ny);

            if (skip_ctfcorrection)  // no correction, simple (W)BP
            {
                cout << "\tSkip CTF correction, use raw micrograph!" << endl;
                // weighting
                if (!skip_weighting) {
                    cout << "\tStart weighting: " << endl;
                    filter_weighting_1D_many(image_now, Nx, stack_orig.getNy(), weighting_radial,
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

                    int n_zz = 0;
                    TDEF(filter_weighting_1D_many)
                    TSTART(filter_weighting_1D_many)
                    /*  fftw线程数     运行时间（s）
                     *      1           7.06
                     *      2           7.01
                     *      4           6.39
                     *      8           6.27
                     *      16及以上线程会导致出粗
                     */
                    // weighting
                    fftwf_plan_with_nthreads(4);
                    if (!skip_weighting) {
                        cout << "\tStart weighting..." << endl;

                        filter_weighting_1D_many(image_now, Nx, Ny, weighting_radial,
                                                 weighting_sigma);
                        cout << "\tDone" << endl;
                    }
                    TEND(filter_weighting_1D_many)
                    TPRINT(filter_weighting_1D_many, "filter_weighting_1D_many time is")
                    filter_weighting_1D_many_all += TINT(filter_weighting_1D_many) ;
                    // 3D-CTF correction
                    cout << "\tStart 3D-CTF correction..." << endl;


                    TDEF(ctf_correction)
                    TSTART(ctf_correction)

                    //loop: number of blocks (for 3D-CTF correction)
                    /*
                     *  这里加openmp 会出错
                     *  浮点数例外
                     */

                    int zz_l = -int(h_tilt_max / 2);
                    int zz_r = int(h_tilt_max / 2);
                    fftwf_plan plan_fft;
                    float *bufc = new float[Nxc * Ny];
                    plan_fft = fftwf_plan_dft_r2c_2d(Ny, Nx, (float *) bufc, reinterpret_cast<fftwf_complex *>(bufc),FFTW_ESTIMATE);
                    buf2fft(image_now, bufc, Nx, Ny);
                    fftwf_execute(plan_fft);
                    fftwf_plan_with_nthreads(1);
                    float defocus1 =  ctf_para[n].getDefocus1();
                    float defocus2 =  ctf_para[n].getDefocus2();
                    float astig =  ctf_para[n].getAstigmatism();
                    float astig_2 = 2*astig;
                    __m512 m_astig_2 = _mm512_set1_ps(astig_2);
                    float lambda =  ctf_para[n].getLambda();
                    float phase_shift =  ctf_para[n].getPhaseShift();
                    float w_sin =  ctf_para[n].getWSin();
                    float w_cos =  ctf_para[n].getWCos();
                    float A = sqrt(w_sin*w_sin+w_cos*w_cos);
                    float B = atan(w_cos/w_sin);
                    float phase_shift_B = B + phase_shift;
                    __m512 m_phase_shift_B = _mm512_set1_ps(phase_shift_B);
                    float pix =  ctf_para[n].getPixelSize();
                    float Cs =  ctf_para[n].getCs();
                    float x_real_l = 1.0 / (Nx * pix);
                    __m512 m_x_real_l=_mm512_set1_ps(x_real_l);
                    float chi_A = M_PI * lambda;
                    __m512 m_chi_A   = _mm512_set1_ps(chi_A);
                    float chi_B = M_PI_2 * Cs * lambda * lambda * lambda;
                    __m512 m_chi_B   = _mm512_set1_ps(-chi_B);
                    float df_now_B = (defocus1 - defocus2)*0.5;
                    __m512 m_df_now_B= _mm512_set1_ps(df_now_B);
                    n_zz += ((zz_r-zz_l-1)/defocus_step+1);
                    fftwf_plan *plan_ifft =  new fftwf_plan[n_zz];
                    for (int n_z=0;n_z<n_zz;n_z++) {
                        plan_ifft[n_z] = fftwf_plan_dft_c2r_2d(Ny, Nx, reinterpret_cast<fftwf_complex *>(stack_corrected[n_z]), (float *) stack_corrected[n_z],
                                                          FFTW_ESTIMATE);
                    }
#pragma omp parallel for
                    for (int n_z=0;n_z<n_zz;n_z++) {
                        memcpy(stack_corrected[n_z],bufc,Nxc * Ny*sizeof(float));
                    }
                    // defocus1,defocus2,astig,lamba,phase_shift,w_sin,w_cos,pix,Cs
#pragma omp parallel for
                    for (int zz = zz_l;
                         zz < zz_r; zz += defocus_step)    // loop over every height (correct with different defocus)
                    {
                        int n_z = (zz + zz_r) / defocus_step;
//                        memcpy(stack_corrected[n_z],bufc,Nxc * Ny*sizeof(float));
//                        fftwf_plan plan_ifft;
//#pragma omp critical
//                        {
//                            plan_ifft = fftwf_plan_dft_c2r_2d(Ny, Nx, reinterpret_cast<fftwf_complex *>(stack_corrected[n_z]), (float *) stack_corrected[n_z],
//                                                              FFTW_ESTIMATE);
//                        }
//                        ctf_correction_perbufc(Nx, Ny, defocus1,defocus2,astig,lambda,phase_shift,w_sin,w_cos,pix,Cs,
//                                               flip_contrast, float(zz) + float(defocus_step - 1) / 2,stack_corrected[n_z]);
                        float z_offset = float(zz) + float(defocus_step - 1) * 0.5;
                        float df_now_A = (defocus1 + defocus2 - 2 * z_offset * pix) * 0.5;
                        __m512 m_df_now_A= _mm512_set1_ps(df_now_A);
                        for (int j = 0; j < Ny; j++) {
                            // loop: Nx+2-Nx%2 (all Fourier components for x-axis)
                            float y=j;
                            float y_norm = (y >= int(ceil(float(Ny + 1) / 2))) ? (y - Ny) : (y);
                            float y_real = float(y_norm) / float(Ny) * (1 / pix);
                           __m512 m_y_real  = _mm512_set1_ps(y_real);
                            // 特殊处理 x_norm = 0
                            {
                                float ctf_now;
                                float x_norm = 0;
                                float x_real = 0;
                                float alpha;
                                if (y_norm > 0) {
                                    alpha = M_PI_2;
                                } else if (y_norm < 0) {
                                    alpha = -M_PI_2;
                                } else {
                                    alpha = 0.0;
                                }
                                float freq2 = y_real * y_real;
                                float df_now = df_now_A + df_now_B * cos(2 * (alpha - astig));
//                                float chi = M_PI * lambda * df_now * freq2 -
//                                            M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2 + phase_shift;
//                                float chi = M_PI * lambda * df_now * freq2 -
//                                            M_PI_2 * Cs * lambda * lambda * lambda * freq2 * freq2 + phase_shift_B;
                                float chi = chi_A * df_now * freq2 -
                                            chi_B * freq2 * freq2 + phase_shift_B;
//                                ctf_now = w_sin * sin(chi) + w_cos * cos(chi);
//                                ctf_now = A*sin(chi);
                                ctf_now = sin(chi);
                                if (flip_contrast) {
                                    ctf_now = -ctf_now;
                                }
                                if (ctf_now >= 0){
                                    stack_corrected[n_z][j * Nxc] *= 1;
                                    stack_corrected[n_z][1 + j * Nxc] *= 1;
                                } else {
                                    stack_corrected[n_z][j * Nxc] *= -1;
                                    stack_corrected[n_z][1 + j * Nxc] *= -1;
                                }
//                                ctf_now = 1-(int)(((*((unsigned int*)&ctf_now)) >> 31) << 1);
//                                stack_corrected[n_z][j * Nxc] *= ctf_now;
//                                stack_corrected[n_z][1 + j * Nxc] *= ctf_now;
                            }


                            int i = 1;
                            for (;i<=(Nxc/2-1)&15;i++){
                                //float ctf_now = ctf.computeCTF2D(i / 2, j, Nx, Ny, true, flip_contrast, z_offset);
                                float ctf_now;
                                float x_real = x_real_l * i;
                                float alpha_astig = atan_pre_2[i+j*(Nxc/2)] - astig_2;
                                float freq2 = x_real * x_real + y_real*y_real;
                                float df_now = df_now_A + df_now_B * cos(alpha_astig);
                                float chi = chi_A * df_now * freq2 - chi_B * freq2 * freq2 + phase_shift_B;
//                              ctf_now = w_sin * sin(chi) + w_cos * cos(chi);
//                              ctf_now = A*sin(chi);
                                ctf_now = sin(chi);
                                if (flip_contrast) {
                                    ctf_now = -ctf_now;
                                }
                                if (ctf_now >= 0){
                                    stack_corrected[n_z][ i*2 + j * Nxc] *= 1;
                                    stack_corrected[n_z][(i*2 + 1) + j * Nxc] *= 1;
                                } else {
                                    stack_corrected[n_z][i*2 + j * Nxc] *= -1;
                                    stack_corrected[n_z][(i*2 + 1) + j * Nxc] *= -1;
                                }
//                                ctf_now = 1-(int)(((*((unsigned int*)&ctf_now)) >> 31) << 1);
//                                stack_corrected[n_z][i * 2 + j * Nxc] *= ctf_now;
//                                stack_corrected[n_z][(i * 2 + 1) + j * Nxc] *= ctf_now;
                            }
                            for (; i <= Nxc/2 -1; i +=16) {
                                //float ctf_now = ctf.computeCTF2D(i / 2, j, Nx, Ny, true, flip_contrast, z_offset);

                                __m512 m_x_real  = _mm512_mul_ps(m_x_real_l,_mm512_set_ps(i+15,i+14,i+13,i+12,i+11,i+10,i+9,i+8,i+7,i+6,i+5,i+4,i+3,i+2,i+1,i));
                                //float x_real = x_real_l * i;

                                __m512 m_alpha   = _mm512_loadu_ps(atan_pre_2+i+j*(Nxc/2));
                                __m512 m_alpha_astig = _mm512_mul_ps(m_alpha,m_astig_2);
                                //float alpha_astig = atan_pre_2[i+j*(Nx/2)] - astig_2;

                                __m512 m_freq2   = _mm512_add_ps(_mm512_mul_ps(m_x_real,m_x_real),_mm512_mul_ps(m_y_real,m_y_real));
                                //float freq2 = x_real * x_real + y_real*y_real;
                                m_alpha_astig    = _mm512_cos_ps(m_alpha_astig);
                                __m512 m_df_now  = _mm512_add_ps(m_df_now_A,_mm512_mul_ps(m_df_now_B,m_alpha_astig));
//                                float df_now = df_now_A + df_now_B * cos(alpha_astig);
                                __m512 m_chi_A_now          = _mm512_mul_ps(m_chi_A,_mm512_mul_ps(m_df_now,m_freq2));
                                __m512 m_chi_B_now          = _mm512_mul_ps(m_chi_B,_mm512_mul_ps(m_freq2,m_freq2));
                                m_chi_A_now          = _mm512_add_ps(m_chi_A_now,m_chi_B_now);
                                __m512 m_chi     = _mm512_add_ps(m_chi_A_now,m_phase_shift_B);
//                                float chi = chi_A * df_now * freq2 - chi_B * freq2 * freq2 + phase_shift_B;
//                              ctf_now = w_sin * sin(chi) + w_cos * cos(chi);
//                              ctf_now = A*sin(chi);
                                __m512 m_ctf_now = _mm512_sin_ps(m_chi);
//                                ctf_now = sin(chi);
                                if (flip_contrast) {
                                    m_ctf_now    = _mm512_mul_ps(m_ctf_now,_mm512_set1_ps(-1));
//                                    ctf_now = -ctf_now;
                                }
                                float *ctf_now_m = new float[16];
                                _mm512_store_ps(ctf_now_m,m_ctf_now);
                                for (int ii=0;ii<16;ii++){
                                    if (ctf_now_m[ii] >= 0){
                                        stack_corrected[n_z][ (i+ii)*2 + j * Nxc] *= 1;
                                        stack_corrected[n_z][((i+ii)*2 + 1) + j * Nxc] *= 1;
                                    } else {
                                        stack_corrected[n_z][(i+ii)*2 + j * Nxc] *= -1;
                                        stack_corrected[n_z][((i+ii)*2 + 1) + j * Nxc] *= -1;
                                    }
                                }

//                                ctf_now = 1-(int)(((*((unsigned int*)&ctf_now)) >> 31) << 1);
//                                stack_corrected[n_z][i * 2 + j * Nxc] *= ctf_now;
//                                stack_corrected[n_z][(i * 2 + 1) + j * Nxc] *= ctf_now;
                            }


//
//                            for (int i=1; i <= Nxc/2 -1; i ++) {
//                                //float ctf_now = ctf.computeCTF2D(i / 2, j, Nx, Ny, true, flip_contrast, z_offset);
//                                float ctf_now;
//                                float x_real = x_real_l * i;
//                                float alpha_astig = atan_pre_2[i+j*(Nx/2)] - astig_2;
//                                float freq2 = x_real * x_real + y_real*y_real;
//                                float df_now = df_now_A + df_now_B * cos(alpha_astig);
//                                float chi = chi_A * df_now * freq2 - chi_B * freq2 * freq2 + phase_shift_B;
////                              ctf_now = w_sin * sin(chi) + w_cos * cos(chi);
////                              ctf_now = A*sin(chi);
//                                ctf_now = sin(chi);
//                                if (flip_contrast) {
//                                    ctf_now = -ctf_now;
//                                }
//                                if (ctf_now >= 0){
//                                    stack_corrected[n_z][ i*2 + j * Nxc] *= 1;
//                                    stack_corrected[n_z][(i*2 + 1) + j * Nxc] *= 1;
//                                } else {
//                                    stack_corrected[n_z][i*2 + j * Nxc] *= -1;
//                                    stack_corrected[n_z][(i*2 + 1) + j * Nxc] *= -1;
//                                }
////                                ctf_now = 1-(int)(((*((unsigned int*)&ctf_now)) >> 31) << 1);
////                                stack_corrected[n_z][i * 2 + j * Nxc] *= ctf_now;
////                                stack_corrected[n_z][(i * 2 + 1) + j * Nxc] *= ctf_now;
//                            }
//


                        }

                        fftwf_execute(plan_ifft[n_z]);
                        //  fft2buf(image, bufc, Nx, Ny);
//#pragma omp simd
//                        for (int i = 0; i < Nxc * Ny; i++)   // normalization
//                        {
//                            stack_corrected[n_z][i] = stack_corrected[n_z][i] / (Nx * Ny);
//                        }
//#pragma omp critical
//                        {
//                            fftwf_destroy_plan(plan_ifft);
//                        }
                    }
#pragma omp parallel for
                    for (int n_z=0;n_z<n_zz;n_z++){
//#pragma omp simd
                        for (int i = 0; i < Nxc * Ny; i++)   // normalization
                        {
                            stack_corrected[n_z][i] = stack_corrected[n_z][i] / (Nx * Ny);
                        }
                       // fftwf_destroy_plan(plan_ifft[n_z]);
                    }

                    for (int n_z=0;n_z<n_zz;n_z++){
                        fftwf_destroy_plan(plan_ifft[n_z]);
                    }
                    delete [] bufc;
                    fftwf_destroy_plan(plan_fft);
                    TEND(ctf_correction)
                    TPRINT(ctf_correction, "ctf_correction time is")
                    ctf_correction_all +=TINT(ctf_correction);
                    cout << "\tDone!" << endl;



                    // recontruction

                    cout << "\tPerform reconstruction:" << endl;


                    TDEF(recon_now)
                    TSTART(recon_now)
                    // loop: Ny (number of xz-slices)
#pragma omp parallel for
                    for (int j = 0; j < Ny; j++) {
                        for (int k = 0; k < h; k++) {
                            float x_orig_k=(float(k) - z_orig_offset) * sin(theta_rad) - x_orig_offset;
                            float z_orig_k=(float(k) - z_orig_offset) * cos(theta_rad);
                            int l=0,r=Nx-1;
                            int l1,r1;
                            if (cos_theta_rad > 0 ){
                                l1 = ceil(x_orig_k / cos_theta_rad + x_orig_offset);
                                r1 = floor(((Nx-1)+x_orig_k)/cos_theta_rad + x_orig_offset);
                            }else{
                                r1 = ceil(x_orig_k / cos_theta_rad + x_orig_offset);
                                l1 = floor(((Nx-1)+x_orig_k)/cos_theta_rad + x_orig_offset);
                            }
                            l=max(l,l1); r=min(r,r1);
                            if (sin_theta_rad > 0 ){
                                l1 = ceil(-(z_orig_k+ int(h_tilt_max/2))/sin_theta_rad + x_orig_offset);
                                r1 = floor(((n_zz * defocus_step) - z_orig_k -int(h_tilt_max/2)) / sin_theta_rad + x_orig_offset - 1);
                            } else {
                                r1 = ceil(-(z_orig_k+ int(h_tilt_max/2))/sin_theta_rad + x_orig_offset);
                                l1 = floor(((n_zz * defocus_step) - z_orig_k -int(h_tilt_max/2)) / sin_theta_rad + x_orig_offset - 1);
                            }
                            l=max(l,l1); r=min(r,r1);

                            int i = l;
                            for (; i <= l -1 +(r - l + 1)&15 ; i++)   // loop for the xz-plane to perform BP
                            {
                                float x_orig = (float(i) - x_orig_offset) * cos_theta_rad - (float(k) - z_orig_offset) * sin_theta_rad + x_orig_offset;
                                float z_orig = (float(i) - x_orig_offset) * sin_theta_rad + (float(k) - z_orig_offset) * cos_theta_rad + z_orig_offset;
                                int x_orig_l = floor(x_orig);
                                int x_orig_r = ceil(x_orig);
                                int n_z = (int)((z_orig - z_orig_offset + int(h_tilt_max / 2) ) / defocus_step);    // the num in the corrected stack for the current height
                                // 上面这行代码一动就精度大问题 例如 int(h_tilt_max / 2) -> zz_r
                                float coeff = x_orig - x_orig_l;
                                stack_recon[j][i + k * Nx] += (1 - coeff) * stack_corrected[n_z][j * Nxc + x_orig_l] + (coeff)*stack_corrected[n_z][j*Nxc + x_orig_r];
                            }
                            for (; i <= r ; i+= 16)   // loop for the xz-plane to perform BP
                            {
                                __m512 m_cos_theta_rad = _mm512_set1_ps(cos_theta_rad);
                                __m512 m_sin_theta_rad = _mm512_set1_ps(sin_theta_rad);
                                __m512 m_x_orig_offset = _mm512_set1_ps(x_orig_offset);
                                __m512 m_z_orig_offset = _mm512_set1_ps(z_orig_offset);
                                __m512 m_orig_a = _mm512_set_ps(i+15 ,i+14,i+13,i+12,i+11,i+10,i+9,i+8,i+7,i+6,i+5,i+4,i+3,i+2,i+1,i);
                                m_orig_a         = _mm512_sub_ps(m_orig_a,m_x_orig_offset);
                                __m512 m_orig_b1 = _mm512_set1_ps(float(k) - z_orig_offset);
                                __m512 m_orig_b2 = _mm512_set1_ps(-float(k) + z_orig_offset);
                                __m512 m_x_orig  = _mm512_mul_ps(m_orig_a,m_cos_theta_rad);
                                m_x_orig         = _mm512_add_ps(m_x_orig,_mm512_mul_ps(m_orig_b2,m_sin_theta_rad));
                                m_x_orig         = _mm512_add_ps(m_x_orig,m_x_orig_offset);
//                                float x_orig = (float(i) - x_orig_offset) * cos_theta_rad - (float(k) - z_orig_offset) * sin_theta_rad + x_orig_offset;
                                __m512 m_z_orig  = _mm512_mul_ps(m_orig_a,m_sin_theta_rad);
                                m_z_orig         = _mm512_add_ps(m_z_orig,_mm512_mul_ps(m_orig_b1,m_cos_theta_rad));
                                m_z_orig         = _mm512_add_ps(m_z_orig,m_z_orig_offset);
//                                float z_orig = (float(i) - x_orig_offset) * sin_theta_rad + (float(k) - z_orig_offset) * cos_theta_rad + z_orig_offset;
                                __m512 m_x_orig_l = _mm512_floor_ps(m_x_orig);
//                                int x_orig_l = floor(x_orig);
                                __m512 m_x_orig_r = _mm512_ceil_ps(m_x_orig);
//                                int x_orig_r = ceil(x_orig);
                                __m512 m_h_tilt_max_1_2 = _mm512_set1_ps(int(h_tilt_max / 2) + 0.0 - z_orig_offset);
                                __m512 m_defocus_step   = _mm512_set1_ps(defocus_step);
                                __m512 m_nz_float       = _mm512_div_ps(_mm512_add_ps(m_z_orig,m_h_tilt_max_1_2),m_defocus_step);
                                __m512 m_nz_int        = _mm512_floor_ps(m_nz_float);
                                // __m512i _mm512_castps_si512 (__m512 a)
                                // int n_z = (int)((z_orig - z_orig_offset + int(h_tilt_max / 2) ) / defocus_step);    // the num in the corrected stack for the current height
                                // 上面这行代码一动就精度大问题 例如 int(h_tilt_max / 2) -> zz_r

                                __m512 m_coeff  = _mm512_sub_ps(m_x_orig,m_x_orig_l);
//                                float coeff = x_orig - x_orig_l;
                                float *coeff = new float[16];
                                _mm512_store_ps(coeff,m_coeff);
                                float *n_z = new float[16];
                                _mm512_store_ps(n_z,m_nz_int);
                                float *x_orig_l = new float[16];
                                _mm512_store_ps(x_orig_l,m_x_orig_l);
                                float *x_orig_r = new float[16];
                                _mm512_store_ps(x_orig_r,m_x_orig_r);
                                for (int ii=0;ii<16;ii++){
                                    stack_recon[j][i+ii + k * Nx] += (1 - coeff[ii]) * stack_corrected[(int)n_z[ii]][j * Nxc + (int)x_orig_l[ii]] + (coeff[ii])*stack_corrected[(int)n_z[ii]][j*Nxc + (int)x_orig_r[ii]];
                                }

                            }
//                            printf("\t\ti : %d\t r : %d\n",i,r);


//                            for (int i = l; i <= r ; i++)   // loop for the xz-plane to perform BP
//                            {
//                                float x_orig = (float(i) - x_orig_offset) * cos_theta_rad - (float(k) - z_orig_offset) * sin_theta_rad + x_orig_offset;
//                                float z_orig = (float(i) - x_orig_offset) * sin_theta_rad + (float(k) - z_orig_offset) * cos_theta_rad + z_orig_offset;
//                                int x_orig_l = floor(x_orig);
//                                int x_orig_r = ceil(x_orig);
//                                int n_z = (int)((z_orig - z_orig_offset + zz_r ) / defocus_step);    // the num in the corrected stack for the current height
//                                // 上面这行代码一动就精度大问题 例如 int(h_tilt_max / 2) -> zz_r
//                                float coeff = x_orig - x_orig_l;
//                                stack_recon[j][i + k * Nx] += (1 - coeff) * stack_corrected[n_z][j * Nxc + x_orig_l] + (coeff)*stack_corrected[n_z][j*Nxc + x_orig_r];
//                            }
                        }
                    }
                    TEND(recon_now)
                    TPRINT(recon_now, "recon_now time is")
                    recon_now_all += TINT(recon_now);

                    cout << "\tDone" << endl;

//                    for (int n_z = 0; n_z < n_zz; n_z++) {
//                        delete[] stack_corrected[n_z];
//                    }
                }
            }
            delete[] image_now;
            delete[] image_now_backup;
            TEND(Image)
            TPRINT(Image, "Image one time is ")
        }
        TEND(n_stack_orig_getNz)
        TPRINT(n_stack_orig_getNz,"n_stack_orig_getNz is ");

        // write out final result
        cout << "Wrtie out final reconstruction result:" << endl;
        MRC stack_final(output_mrc.c_str(), "wb");
        stack_final.createMRC_empty(stack_orig.getNx(), h, Ny, 2); // (x,z,y)
        // loop: Ny (number of xz-slices)
        for (int j = 0; j < Ny; j++) {
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
        for (int j = 0; j < Ny; j++) {
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
        mean_all /= Ny;
        stack_final.computeHeader(pix, false, min_all, max_all, float(mean_all));

        for (int j = 0; j < Ny; j++) {
            delete[] stack_recon[j];
        }

        stack_final.close();
        cout << "Done" << endl;
        TEND(Reconstruction_ALL)
        TPRINT(Reconstruction_ALL, "Reconstruction_ALL time is")

    }

    stack_orig.close();

    cout << endl << "Finish reconstruction successfully!" << endl;
    cout << "All results save in: " << path << endl << endl;

    printf("\t\tfilter weighting 1D many all time is %f s\n",(float)filter_weighting_1D_many_all/1e6);
    printf("\t\tctf correction all time is %f s\n",(float)ctf_correction_all/1e6);
    printf("\t\trecon now all time is %f s\n",(float)recon_now_all/1e6);


}