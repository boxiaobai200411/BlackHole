#include "Graph.h"
#include "tools.hpp"
#include <bits/stdc++.h>
#include <chrono>
using namespace std;
string path = "../dataset/";
string dataid;
const int N = 1e5 + 100;
string output_dst_file_path;
int k;

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
int main(int argc, char **argv) {

    ui l = -1;
    ui r = -1;
    int maxset = 0x3f3f3f3f;

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "-d" && i + 1 < argc) {
            dataid = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            l = atoi(argv[++i]);
        } else if (arg == "-r" && i + 1 < argc) {
            r = atoi(argv[++i]);
        } else if (arg == "-m" && i + 1 < argc) {
            maxset = atoi(argv[++i]);
        } else {
            cerr << "Unknown or incomplete argument: " << arg << endl;
            return 1;
        }
    }
    if (dataid.empty() || l == -1 || r == -1) {
        cerr << "Usage: ./program -d <dataset> -l <left> -r <right> [-m "
                "<minBHSize>]"
             << endl;
        return 1;
    }
    readFileTrans(path, dataid);
    Graph graph;
    bool build_flag = 0;
    build_flag = graph.loadIndex(dataid);
    if (build_flag == 0) {
        graph.loadGraphFromFile("../dataset/transform/" + dataid + "_tran.txt");
        graph.buildIndex(maxset);
        graph.saveIndex(dataid);
    }
    cout << "overbuild\n";

    for (ui i = l; i <= r; i++) {

        ofstream ans("../answer.csv", ios::app);
        double sum_ser = 0;
        double sum_op = 0;
        ui sum_anss = 0;
        for (ui j = 0; j <= 10; j++) {
            auto t0 = std::chrono::steady_clock::now();
            ul sum_ans = graph.queryKBlackHoleParallel(i, "/dev/null"); //
            auto t1 = std::chrono::steady_clock::now();
            double elapsed_s = std::chrono::duration<double>(t1 - t0).count();
            sum_op += elapsed_s;
            auto t2 = std::chrono::steady_clock::now();
            ul sum_ans2 = graph.queryKBlackHole(i, "/dev/null");
            auto t3 = std::chrono::steady_clock::now();
            double elapsed_s2 = std::chrono::duration<double>(t3 - t2).count();
            sum_ser += elapsed_s2;
            sum_anss = sum_ans;
        }

        cout.setf(std::ios::fixed);
        cout << setprecision(6);
        cout << dataid << "," << i << ",ans:" << sum_anss << ","
             << sum_op / 10.0 << ",";
        cout << "" << sum_ser / 10.0 << "\n";
        if (sum_anss != 0) {
            ans << dataid << "," << i << "," << sum_anss << "," << sum_op / 10.0
                << "," << sum_ser / 10.0 << "\n";
        }

        // ofstream ans("../BaseLine_answer.csv", ios::app);
        // auto t0 = std::chrono::steady_clock::now();
        // ul sum_ans = graph.baseline(i, "/dev/null");
        // auto t1 = std::chrono::steady_clock::now();
        // double elapsed_s2 = std::chrono::duration<double>(t1 - t0).count();
        // cout << i << " " << elapsed_s2 << "\n";
        // ans << dataid << "," << i << "," << elapsed_s2 << "\n";

        ans.close();
    }

    map<string, size_t> mem = get_index_mem();
    if (mem.find("pk") != mem.end()) {
        cout << "VmPeak: " << mem["pk"] << " KB" << endl;
        cout << "VmPeak: " << mem["pk"] / 1024.0 << " MB" << endl; // 转 MB
        ofstream ans("../answer_MEM.csv", ios::app);
        ans << dataid << "," << mem["pk"] / 1024.0 << "MB\n";
    } else {
        cout << "VmPeak not found" << endl;
    }
    return 0;
}
