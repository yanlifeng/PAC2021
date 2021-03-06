/*******************************************************************
 *       Filename:  tomo.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/15/2020 05:48:48 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:
 *          Email:
 *        Company:
 *
 *******************************************************************/
#include <stdio.h>
#include "util.h"
#include "mrc.h"
#include "CTFAlgo.h"
#include "ReconstructionAlgo_WBP_RAM.h"
#include "time.h"
#include <sys/time.h>

double GetTime2() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

void getTime(time_t start_time) {
    time_t end_time;
    time(&end_time);

    double seconds_total = difftime(end_time, start_time);
    int hours = ((int) seconds_total) / 3600;
    int minutes = (((int) seconds_total) % 3600) / 60;
    int seconds = (((int) seconds_total) % 3600) % 60;

    cout << "Time elapsed: ";
    if (hours > 0) {
        cout << hours << "h ";
    }
    if (minutes > 0 || hours > 0) {
        cout << minutes << "m ";
    }
    cout << seconds << "s" << endl << endl;
}

int main(int argc, char **argv) {

    if (argc != 2) {
        fprintf(stderr, "\n  Usage: \n%6s%s /path/to/parameter/file\n\n", "", argv[0]);
        exit(1);
    }


    time_t start_time;
    time(&start_time);

    double t_start = GetTime2();

    map<string, string> inputPara;
    map<string, string> outputPara;
    const char *paraFileName = argv[1];
    readParaFile(inputPara, paraFileName);
    getAllParas(inputPara);

    bool do_alignment = 0, do_CTF = 0, do_reconstruction_WBP = 0, do_reconstruction_WBP_in_RAM = 0, do_reconstruction_SIRT = 0, do_reconstruction_SIRT_in_RAM = 0, do_reconstruction_FD = 0;
    map<string, string>::iterator it = inputPara.find("do_alignment");
    if (it != inputPara.end()) {
        do_alignment = atoi(it->second.c_str());
    }
    it = inputPara.find("do_CTF");
    if (it != inputPara.end()) {
        do_CTF = atoi(it->second.c_str());
    }

    it = inputPara.find("do_reconstruction_WBP_in_RAM");
    if (it != inputPara.end()) {
        do_reconstruction_WBP_in_RAM = atoi(it->second.c_str());
    }

    if (do_CTF) {
        cout << endl << "Do CTF!" << endl << endl;
        CTFBase *cb = new CTFAlgo();
        cb->doCTF(inputPara, outputPara);
        delete cb;
    }

    if (do_reconstruction_WBP_in_RAM) {
        cout << endl << "Do Reconstruction with WBP in RAM!" << endl << endl;
        ReconstructionBase *rb = new ReconstructionAlgo_WBP_RAM();
        rb->doReconstruction(inputPara, outputPara);
        delete rb;
    }


    cout << endl << "Finish!" << endl;
    getTime(start_time);
    printf("tot cost %.3f\n", GetTime2() - t_start);

    return 0;
}
