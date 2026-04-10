#include "type.h"
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
#define ll long long

bool readFileTrans(const string &path, const string &dataid) {
    ifstream ifs(path + dataid + ".txt");
    if (!ifs.is_open()) {
        cerr << "Cannot open file: " << path + dataid + ".txt" << "\n";
        return false;
    }
    vector<pair<VertexID, VertexID>> ed;
    map<VertexID, VertexID> mp, hashMp;
    ui n, m, u, v;
    ul t;
    vector<pair<ul, pair<VertexID, VertexID>>> nowedge;
    while (ifs >> u >> v >> t) {
        nowedge.push_back({t, {u, v}});
        mp[u] = 1;
        mp[v] = 1;
    }
    sort(nowedge.begin(), nowedge.end());
    ul startTime = nowedge[0].first;

    ui idx = 0;
    m = nowedge.size();
    n = mp.size();
    ui *inDeg = new ui[n + 10]();
    ui *outDeg = new ui[n + 10]();
    ui *labels = new ui[n + 10]();

    for (auto it : mp) {
        labels[idx] = it.first;
        hashMp[it.first] = idx;
        idx++;
    }

    for (ui i = 0; i < m; i++) {
        u = hashMp[nowedge[i].second.first];
        v = hashMp[nowedge[i].second.second];
        nowedge[i].second.first = u;
        nowedge[i].second.second = v;
        inDeg[v]++;
        outDeg[u]++;
        nowedge[i].first = 0;
    }

    string outFile = "../dataset/transform/" + dataid + "_tran.txt";
    ofstream ofs(outFile);
    ofs << idx << " " << m << "\n";
    for (ui i = 0; i < idx; i++) {
        ofs << "v " << i << " " << labels[i] << " " << outDeg[i] << " "
            << inDeg[i] << "\n";
    }
    for (ui j = 0; j < m; j++) {

        ofs << "e " << nowedge[j].second.first << " "
            << nowedge[j].second.second << " \n";
    }
    cout << "transOver\n";
    return 1;
}

bool readFileTrans_NOTime(const string &path, const string &dataid) {
    ifstream ifs(path + dataid + ".txt");
    if (!ifs.is_open()) {
        cerr << "Cannot open file: " << path + dataid + ".txt" << "\n";
        return false;
    }
    vector<pair<VertexID, VertexID>> ed;
    map<VertexID, VertexID> mp, hashMp;
    ui n, m, u, v;
    ul t;
    vector<pair<ul, pair<VertexID, VertexID>>> nowedge;
    while (ifs >> u >> v) {
        nowedge.push_back({0, {u, v}});
        mp[u] = 1;
        mp[v] = 1;
    }
    sort(nowedge.begin(), nowedge.end());
    ul startTime = nowedge[0].first;

    ui idx = 0;
    m = nowedge.size();
    n = mp.size();
    ui *inDeg = new ui[n + 10]();
    ui *outDeg = new ui[n + 10]();
    ui *labels = new ui[n + 10]();

    for (auto it : mp) {
        labels[idx] = it.first;
        hashMp[it.first] = idx;
        idx++;
    }

    for (ui i = 0; i < m; i++) {
        u = hashMp[nowedge[i].second.first];
        v = hashMp[nowedge[i].second.second];
        nowedge[i].second.first = u;
        nowedge[i].second.second = v;
        inDeg[v]++;
        outDeg[u]++;
        nowedge[i].first = 0;
    }

    string outFile = "../dataset/transform/" + dataid + "_tran.txt";
    ofstream ofs(outFile);
    ofs << idx << " " << m << "\n";
    for (ui i = 0; i < idx; i++) {
        ofs << "v " << i << " " << labels[i] << " " << outDeg[i] << " "
            << inDeg[i] << "\n";
    }
    for (ui j = 0; j < m; j++) {

        ofs << "e " << nowedge[j].second.first << " "
            << nowedge[j].second.second << " \n";
    }
    cout << "transOver\n";
    return 1;
}

map<string, size_t> get_index_mem() {
    FILE *fp = fopen("/proc/self/status", "r");
    char line[128];
    map<string, size_t> res;
    while (fgets(line, 128, fp) != NULL) {
        if (strncmp(line, "VmPeak:", 7) == 0) {
            string p = line;
            res["pk"] = size_t(stoull(p.substr(7)));
            // cout << line << endl;
            // printf("当前进程占用虚拟内存大小为：%d KB\n", atoi(line + 6));
        }
    }
    fclose(fp);
    return res;
}
void get_code_mem(string dataid) {
    map<string, size_t> mem = get_index_mem();
    if (mem.find("pk") != mem.end()) {
        cout << "VmPeak: " << mem["pk"] << " KB" << endl;
        cout << "VmPeak: " << mem["pk"] / 1024.0 << " MB" << endl; // 转 MB
        ofstream ans("../answer_MEM.csv", ios::app);
        ans << dataid << "," << mem["pk"] / 1024.0 << "MB\n";
    } else {
        cout << "VmPeak not found" << endl;
    }
}