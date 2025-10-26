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

    while (ifs >> u >> v) {
        ed.push_back({u, v});
        mp[u] = 1;
        mp[v] = 1;
    }
    ui idx = 0;
    m = ed.size();
    n = mp.size();
    ui *inDeg = new ui[n + 10];
    ui *outDeg = new ui[n + 10];
    ui *labels = new ui[n + 10];

    for (auto it : mp) {
        labels[idx] = it.first;
        hashMp[it.first] = idx;
        idx++;
    }

    for (ui i = 0; i < m; i++) {
        u = hashMp[ed[i].first];
        v = hashMp[ed[i].second];
        ed[i].first = u;
        ed[i].second = v;
        inDeg[v]++;
        outDeg[u]++;
    }

    string outFile = "../dataset/transform/" + dataid + "_tran.txt";
    ofstream ofs(outFile);
    ofs << idx << " " << m << "\n";
    for (ui i = 0; i < idx; i++) {
        ofs << "v " << i << " " << labels[i] << " " << outDeg[i] << " "
            << inDeg[i] << "\n";
    }
    for (ui j = 0; j < m; j++) {
        ofs << "e " << ed[j].first << " " << ed[j].second << "\n";
    }
    cout << "transOver\n";
    return 1;
}
