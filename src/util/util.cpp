/*******************************************************************
 *       Filename:  util.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/16/2020 10:26:57 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ruan Huabin
 *          Email:  ruanhuabin@tsinghua.edu.cn
 *        Company:  Dep. of CS, Tsinghua Unversity
 *
 *******************************************************************/
#include "util.h"

string &trim(std::string &s) {
    if (s.empty()) {
        return s;
    }

    s.erase(0, s.find_first_not_of(" "));
    s.erase(s.find_last_not_of(" ") + 1);
    return s;
}

void readParaFile(map <string, string> &paraMap, const char *fileName) {
    ifstream ifs(fileName, ios::in);

    if (!ifs) {
        fprintf(stderr, "Open parameterfile [ %s ] failed\n", fileName);
        abort();
    }

    string s;
    while (getline(ifs, s)) {
        size_t pos = s.find_first_of("=");
        if (pos == string::npos) {
            fprintf(stderr, "Error: No '=' found in [%s] in parameter file [%s]\n", s.c_str(), fileName);
        }
        string key = s.substr(0, pos);
        string value = s.substr(pos + 1);
        trim(key);
        trim(value);
        paraMap[key] = value;
    }
}

string getPara(map <string, string> paraMap, string key) {
    return "";
}

void getAllParas(map <string, string> paraMap) {
    map<string, string>::iterator iter;
    for (iter = paraMap.begin(); iter != paraMap.end(); iter++) {
        cout << iter->first << "=" << iter->second << endl;
    }
}
