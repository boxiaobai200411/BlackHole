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

int main(int argc, char **argv) {

    ui l = -1;
    ui r = -1;
    ui z = -1;
    int maxy = 24;
    double up = -1;
    int maxset = 0x3f3f3f3f;
    int readT = 0;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "-d" && i + 1 < argc) {
            dataid = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            l = atoi(argv[++i]);
        } else if (arg == "-r" && i + 1 < argc) {
            r = atoi(argv[++i]);
        } else if (arg == "-e" && i + 1 < argc) {
            z = atoi(argv[++i]);
        }else if (arg == "-y" && i + 1 < argc) {
            maxy = atoi(argv[++i]);
        } 
        else if (arg == "-t" && i + 1 < argc) {
            readT = stoi(argv[++i]);

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
    if (readT == 0) {
        readFileTrans_NOTime(path, dataid);
    } else if (readT == 1) {
        readFileTrans(path, dataid);
    }

    Graph graph;
    maxset = r;
    bool loadIdex = graph.loadIndex(dataid, maxset);
    if (!loadIdex) {
        auto buildTimeStart = std::chrono::steady_clock::now();
        graph.loadGraphFromFile("../dataset/transform/" + dataid + "_tran.txt");
        graph.buildIndex(maxset);
        auto buildTimeEnd = std::chrono::steady_clock::now();
        double buildTime =
            std::chrono::duration<double>(buildTimeEnd - buildTimeStart)
                .count();
        graph.saveIndex(dataid, maxset);
        ofstream indexTime("../index/buildTime.csv", ios::app);
        indexTime << dataid << "," << maxset << "," << buildTime << "\n";
        indexTime.close();
    }
    cout << "overbuild\n";
    for (ui i = l; i <= r; i++) {
        ofstream ans("../answer.csv", ios::app);
        for (ui eachin = 1; eachin <= z; eachin++) {
            bool endflag = 0;
            for (double y = eachin; y <= maxy; y += 1) {
                if (y < eachin || endflag == 1) {
                    continue;
                }
                double sum_elapsed_s = 0, sum_elapsed_s2 = 0;
                ul sum_ans;
                for (int repeat = 0; repeat < 10; repeat++) {

                    auto t0 = std::chrono::steady_clock::now();
                    sum_ans = graph.queryKBlackHoleParallel(i, y, eachin,
                                                            "/dev/null");
                    auto t1 = std::chrono::steady_clock::now();
                    sum_elapsed_s2 +=
                        std::chrono::duration<double>(t1 - t0).count();

                    auto t2 = std::chrono::steady_clock::now();
                    sum_ans = graph.queryKBlackHole(i, y, eachin,
                    "/dev/null"); auto t3 = std::chrono::steady_clock::now();
                    sum_elapsed_s +=
                        std::chrono::duration<double>(t3 - t2).count();
                }

                cout.setf(std::ios::fixed);
                cout << setprecision(6);

                if (sum_ans != 0)
                    cout << dataid << ",k:" << i << ",eachin:" << eachin
                         << ",y:" << y << ",ans:" << sum_ans
                         << ",time1:" << sum_elapsed_s / 10.0
                         << "s,time2:" << sum_elapsed_s2 / 10.0 << "\n";
                if (sum_ans != 0) {
                    ans << dataid << "," << i << "," << eachin << "," << y
                        << "," << sum_ans << "," << sum_elapsed_s / 10.0 <<
                        ","
                        << sum_elapsed_s2 / 10.0 << "\n";
                }

                if (sum_ans == 0) {
                    cout << dataid << ",k:" << i << ",eachin:" << eachin
                         << ",y:" << y << " No Answer\n";
                    endflag = 1;
                }
            }
        }
    }
    get_code_mem(dataid);
    return 0;
}

// }
