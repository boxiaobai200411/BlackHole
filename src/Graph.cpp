#include "Graph.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <queue>
#include <set>
#include <unordered_set>
#include <vector>
#ifdef _OPENMP


#endif
using namespace std;


#define epsBuild 0


using Clock = std::chrono::steady_clock;
using Seconds = std::chrono::duration<double>;
void Graph::loadGraphFromFile(const std::string &file_path) {
    std::ifstream infile(file_path);
    if (!infile.is_open()) {
        std::cout << "Can not open the graph file " << file_path << " ."
                  << std::endl;
        exit(-1);
    }

    char type;
    infile >> vertices_count_ >> edges_count_;
    offsets_ = new ul[vertices_count_ + 1];
    offsets_[0] = 0;
    neighbors_ = new VertexID[edges_count_];

    reverse_offsets_ = new ul[vertices_count_ + 1];
    reverse_offsets_[0] = 0;
    reverse_neighbors_ = new VertexID[edges_count_];

    labels_ = new LabelID[vertices_count_];

    std::vector<ui> neighbors_offset(vertices_count_, 0);
    std::vector<ui> reverse_neighbors_offset(vertices_count_, 0);

    while (infile >> type) {
        if (type == 'v') { // Read vertex.
            VertexID id;
            LabelID label;
            ui outDegree, inDegree;
            infile >> id >> label >> outDegree >> inDegree;
            labels_[id] = label;
            offsets_[id + 1] = offsets_[id] + outDegree;
            reverse_offsets_[id + 1] = reverse_offsets_[id] + inDegree;
        } else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
            ui nowtime;
            infile >> begin >> end;
            ui offset = offsets_[begin] + neighbors_offset[begin];
            neighbors_[offset] = end;
            neighbors_offset[begin] += 1;

            ui reverse_offset =
                reverse_offsets_[end] + reverse_neighbors_offset[end];
            reverse_neighbors_[reverse_offset] = begin;
            reverse_neighbors_offset[end] += 1;
        }
    }
    infile.close();

    std::cout << "lord_over\n";
    ui count;
}

void Graph::buildSCCGraph() {

    cout << "buildSCCGraph\n";
    std::vector<ui> dfn(vertices_count_, 0), low(vertices_count_, 0), stk,
        dfs_stk;
    std::vector<bool> inStack(vertices_count_, false);
    std::vector<ui> next_neighbor(vertices_count_, 0);

    scc_id_ = new ui[vertices_count_];
    std::fill(scc_id_, scc_id_ + vertices_count_, (ui)-1);

    ui timer = 0, scc_cnt = 0;

    for (ui start = 0; start < vertices_count_; start++) {
        if (dfn[start])
            continue;

        dfs_stk.push_back(start);
        next_neighbor[start] = 0;

        while (!dfs_stk.empty()) {
            ui u = dfs_stk.back();

            if (dfn[u] == 0) {

                dfn[u] = low[u] = ++timer;
                stk.push_back(u);
                inStack[u] = true;
            }

            ui deg;
            const ui *nbr = getVertexNeighbors(u, deg);
            bool found_unvisited = false;

            for (ui i = next_neighbor[u]; i < deg; i++) {
                ui v = nbr[i];
                next_neighbor[u] = i + 1;

                if (dfn[v] == 0) {

                    dfs_stk.push_back(v);
                    next_neighbor[v] = 0;
                    found_unvisited = true;
                    break;
                } else if (inStack[v]) {
                    low[u] = std::min(low[u], dfn[v]);
                }
            }

            if (!found_unvisited) {

                dfs_stk.pop_back();

                if (!dfs_stk.empty()) {
                    ui parent = dfs_stk.back();
                    low[parent] = std::min(low[parent], low[u]);
                }

                if (dfn[u] == low[u]) {
                    while (true) {
                        ui v = stk.back();
                        stk.pop_back();
                        inStack[v] = false;
                        scc_id_[v] = scc_cnt;
                        if (v == u)
                            break;
                    }
                    scc_cnt++;
                }
            }
        }
    }
    scc_count_ = scc_cnt;

    std::vector<std::vector<ui>> new_adj(scc_count_), new_rev_adj(scc_count_);
    for (ui u = 0; u < vertices_count_; u++) {
        ui u_scc = scc_id_[u];
        ui deg;
        const ui *nbr = getVertexNeighbors(u, deg);
        for (ui i = 0; i < deg; i++) {
            ui v = nbr[i];
            ui v_scc = scc_id_[v];
            if (u_scc != v_scc) {
                new_adj[u_scc].push_back(v_scc);
                new_rev_adj[v_scc].push_back(u_scc);
            }
        }
    }
    auto dedup = [](std::vector<ui> &vec) {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    };
    for (ui i = 0; i < scc_count_; i++) {
        dedup(new_adj[i]);
        dedup(new_rev_adj[i]);
    }

    scc_offsets_ = new ul[scc_count_ + 1];
    scc_offsets_[0] = 0;
    ui total_edges = 0;
    for (ui i = 0; i < scc_count_; i++) {
        total_edges += new_adj[i].size();
        scc_offsets_[i + 1] = total_edges;
    }
    scc_neighbors_ = new VertexID[total_edges];
    ui idx = 0;
    for (ui i = 0; i < scc_count_; i++) {
        for (ui v : new_adj[i]) {
            scc_neighbors_[idx++] = v;
        }
    }

    scc_rev_offsets_ = new ul[scc_count_ + 1];
    scc_rev_offsets_[0] = 0;
    total_edges = 0;
    for (ui i = 0; i < scc_count_; i++) {
        total_edges += new_rev_adj[i].size();
        scc_rev_offsets_[i + 1] = total_edges;
    }
    scc_rev_neighbors_ = new VertexID[total_edges];
    idx = 0;
    for (ui i = 0; i < scc_count_; i++) {
        for (ui v : new_rev_adj[i]) {
            scc_rev_neighbors_[idx++] = v;
        }
    }

    std::vector<std::vector<ui>> scc_vertices_tmp(scc_count_);
    for (ui v = 0; v < vertices_count_; v++) {
        scc_vertices_tmp[scc_id_[v]].push_back(v);
    }

    scc_vertices_offsets_ = new ul[scc_count_ + 1];
    scc_vertices_offsets_[0] = 0;
    ui total_v = 0;
    for (ui i = 0; i < scc_count_; i++) {
        total_v += scc_vertices_tmp[i].size();
        scc_vertices_offsets_[i + 1] = total_v;
    }
    scc_vertices_ = new VertexID[total_v];
    idx = 0;
    for (ui i = 0; i < scc_count_; i++) {
        for (ui v : scc_vertices_tmp[i]) {
            scc_vertices_[idx++] = v;
        }
    }

    cout << "before:" << vertices_count_ << " after:" << scc_count_ << "\n";
}

void Graph::getTheMinBlackHoleBFS(int maxset) {

    vector<pair<pair<ui, ui>, vector<ui>>> answer;
    vector<ui> vis(scc_count_);
    ui nowtime = 0;
    for (ui i = 0; i < scc_count_; i++) {
        if (i % (scc_count_ / 10) == 0) {
            std::cout << "Find MinBlackHole over:"
                      << (double)i / scc_count_ * 100.0 << "%\n";
        }
        queue<ui> que;
        vector<ui> nowans;
        nowans.reserve(maxset + 10);
        que.push(i);
        nowtime++;
        vis[i] = nowtime;
        ui flag = 1;
        while (!que.empty()) {
            VertexID u = que.front();
            ui count;
            que.pop();
            const ui *neighbors = getSCCVertices(u, count);
            for (ui j = 0; j < count; j++) {
                nowans.push_back(neighbors[j]);
                if (nowans.size() > maxset)
                    break;
            }
            if (nowans.size() > maxset)
                break;
            neighbors = getSCCNeighbors(u, count);
            for (ui j = 0; j < count; j++) {
                VertexID v = neighbors[j];
                if (vis[v] != nowtime) {
                    vis[v] = nowtime;
                    que.push(v);
                    flag++;
                }
                if (flag > maxset)
                    break;
            }
            if (nowans.size() > maxset || flag > maxset)
                break;
        }
        if (nowans.size() > maxset || flag > maxset) {
            continue;
        }


        ul minIndeg = 0x3f3f3f3f;
        for (auto ver : nowans) {
            ui count;
            ul Indeg = 0;
            const ui *neighborsVer = getVertexReverseNeighbors(ver, count);
            for (int j = 0; j < count; j++) {

                if (vis[scc_id_[neighborsVer[j]]] == nowtime)
                    continue;
                Indeg++;
            }
            minIndeg = min(minIndeg, Indeg);
            if (minIndeg <= 0) {
                break;
            }
        }
#if epsBuild == 1
        if (minIndeg == 0)
            continue;
#endif
        

        sort(nowans.begin(), nowans.end());

        answer.push_back(
            {{nowans.size(), minIndeg},
             nowans}); // first.first=size,first.second=Indeg,second=vec
    }

    std::sort(answer.begin(), answer.end(),
              [](const pair<pair<ui, ui>, vector<ui>> &A,
                 const pair<pair<ui, ui>, vector<ui>> &B) {
                  if (A.first.first != B.first.first)
                      return A.first.first > B.first.first;
                  return A.first.second < B.first.second;
              });
    scc_count_ = answer.size();

    partBlackHole_offsets_ = new ul[scc_count_ + 1];
    partBlackHole_minIndeg_ = new ui[scc_count_ + 1];
    partBlackHole_offsets_[0] = 0;

    for (ui i = 1; i <= scc_count_; i++) {
        partBlackHole_offsets_[i] =
            partBlackHole_offsets_[i - 1] + answer[i - 1].first.first;
        partBlackHole_minIndeg_[i - 1] = answer[i - 1].first.second;
    }

    partBlackHole_neighbors_ = new VertexID[partBlackHole_offsets_[scc_count_]];

    for (ui i = 0; i < scc_count_; i++) {
        ui count = answer[i].first.first;
        ul start = partBlackHole_offsets_[i];
        for (ui j = 0; j < count; j++) {
            partBlackHole_neighbors_[partBlackHole_offsets_[i] + j] =
                answer[i].second[j];
        }
    }

    cout << partBlackHole_offsets_[scc_count_] << " size\n";
    cout << "Find MinBlackHole_over\n";
}

void Graph::connectSCC(int maxset) {

    vector<vector<ui>> scc_neighbors(scc_count_);
    for (ui i = 0; i < scc_count_; i++) {
        ui len;
        const ui *neighbors = getByPartBHVerticesNeighbors(i, len);
        scc_neighbors[i].assign(neighbors, neighbors + len);
    }

    vector<vector<ui>> edge(scc_count_);

    for (ui i = 0; i < scc_count_; i++) {
        if (i % (scc_count_ / 10) == 0) {
            std::cout << "connect SCC over:" << (double)i / scc_count_ * 100.0
                      << "%\n";
        }
        if (scc_neighbors[i].size() > maxset || scc_neighbors[i].size() == 0) {
            cout << "sbsbsbsbsb\n";
            continue;
        }

        if (scc_neighbors[i].size() == 1)
            break;

        for (ui j = i + 1; j < scc_count_; j++) {

            if (scc_neighbors[j].size() == 1)
                break;

            ui flag = checkRelation(scc_neighbors[i], scc_neighbors[j]);
            if (flag == 1) {
                {
                    edge[i].push_back(j);
                }
            }
        }
    }
    ui cu = 0;
    partBlackHoleConnect_offsets_ = new ul[scc_count_ + 1];
    partBlackHoleConnect_offsets_[0] = 0;
    for (ui i = 1; i <= scc_count_; i++) {
        partBlackHoleConnect_offsets_[i] =
            partBlackHoleConnect_offsets_[i - 1] + edge[i - 1].size();
    }
    partBlackHoleConnect_neighbors_ =
        new VertexID[partBlackHoleConnect_offsets_[scc_count_]];
    for (ui i = 0; i < scc_count_; i++) {
        ul start = partBlackHoleConnect_offsets_[i];
        ui count = edge[i].size();
        for (ui j = 0; j < count; j++) {

            partBlackHoleConnect_neighbors_[start + j] = edge[i][j];
            cu++;
        }
    }

    cout << "Connect BlackHole edge:" << cu << "\n";
}

void Graph::buildIndex(int maxset) {
    buildSCCGraph();
    getTheMinBlackHoleBFS(maxset);
    connectSCC(maxset);
}
#ifdef COUTTIME
double TimeSelectAnswer = 0, combinset = 0, combinnei = 0, allreadyTime = 0,
       createTime = 0;
#endif

void Graph::searchAnswer(ui level, ui start, const ui k, ui answerCount,
                         VertexID *answer_, bool *visPoint, ui NeighborCount,
                         VertexID *NeighborSet, bool *VisNeighbor,
                         ui upPointCount, VertexID *upPoint, double y,
                         int eachin, std::ofstream &ofs, uint64_t &kSizeCount,
                         BloomFilter &localDistinct) {

    if (answerCount == k) {
#ifdef COUTTIME
        auto tans1ex = std::chrono::steady_clock::now();
#endif
        bool alreadyExists = localDistinct.checkAndInsertArray(
            answer_, answerCount); // 看遇到过吗
#ifdef COUTTIME
        auto tans2ex = std::chrono::steady_clock::now();
        allreadyTime +=
            std::chrono::duration<double>(tans2ex - tans1ex).count();
#endif

        if (!alreadyExists) {
#ifdef COUTTIME
            auto tans1 = std::chrono::steady_clock::now();
#endif

            ui outCount = 0;
            bool eachin_flag = 1;
            for (ui i = 0; i < answerCount; i++) {
                ui count;
                const ui *neighbor =
                    getVertexReverseNeighbors(answer_[i], count);
                ui nowCount = 0;
                for (ui j = 0; j < count; j++) {
                    VertexID v = neighbor[j];
                    if (visPoint[v])
                        continue;
                    nowCount++;
                    outCount++;
                }
                if (nowCount < eachin)
                    eachin_flag = 0;
            }
            double nowy = outCount / (double)k;
            if (nowy >= y && eachin_flag == 1) {
                kSizeCount++;
            }
#ifdef COUTTIME
            auto tans2 = std::chrono::steady_clock::now();
            TimeSelectAnswer +=
                std::chrono::duration<double>(tans2 - tans1).count();
#endif
        }

        // for(ui i=0;i<answerCount;i++){
        //     ofs<<answer_[i]<<" ";
        // }
        // ofs<<"\n";
        return;
    } else if (answerCount >= k) {
        return;
    }

    for (ui i = level; i < NeighborCount; i++) {
        ui addPoint = 0;
        ui addPointNeighbor = 0;
        ui count;
        VertexID u = NeighborSet[i];
#ifdef COUTTIME
        auto tans1 = std::chrono::steady_clock::now();
#endif

        const ui *neighbor = getByPartBHVerticesNeighbors(u, count);
        for (ui j = 0; j < count; j++) {
            VertexID v = neighbor[j];
            if (visPoint[v])
                continue;
            visPoint[v] = true;
            answer_[answerCount + addPoint] = v; // 当前点放进去
            addPoint++;
        }

        bool flagInDegree = 1;
        for (ui j = 0; j < addPoint; j++) {
            ui nowVertexinCount = 0;

            const ui *neighbornow =
                getVertexReverseNeighbors(answer_[answerCount + j],
                                          count); // 这个点的反向邻居

            for (ui inNei = 0; inNei < count; inNei++) {
                if (visPoint[neighbornow[inNei]])
                    continue;
                nowVertexinCount++;
            }

            if (nowVertexinCount < eachin) {
                flagInDegree = 0;
                break;
            }
        }

        if (flagInDegree == 0) {
            // cout << "NO\n";
            for (ui j = 0; j < addPoint; j++) {
                visPoint[answer_[answerCount + j]] = false;
            }
            // cout << "YES\n";
            continue;
        }
#ifdef COUTTIME
        auto tans2 = std::chrono::steady_clock::now();
        combinset += std::chrono::duration<double>(tans2 - tans1).count();
        auto tans11 = std::chrono::steady_clock::now();
#endif

        neighbor = getPartBlackHoleConnectNeighbors(u, count);

        for (ui j = 0; j < count; j++) {
            VertexID v = neighbor[j];
            if (VisNeighbor[v] || v <= start )
                continue;

#if epsBuild
            if(partBlackHole_minIndeg_[v] < eachin) continue;
#endif

            VisNeighbor[v] = true;
            NeighborSet[NeighborCount + addPointNeighbor] = v;
            addPointNeighbor++;
        }
#ifdef COUTTIME
        auto tans22 = std::chrono::steady_clock::now();
        combinnei += std::chrono::duration<double>(tans22 - tans11).count();
#endif

        upPoint[upPointCount] = u;

        if (addPoint != 0)
            searchAnswer(i + 1, start, k, answerCount + addPoint, answer_,
                         visPoint, NeighborCount + addPointNeighbor,
                         NeighborSet, VisNeighbor, upPointCount + 1, upPoint, y,
                         eachin, ofs, kSizeCount, localDistinct);

#ifdef COUTTIME
        tans1 = std::chrono::steady_clock::now();
#endif

        for (ui j = 0; j < addPoint; j++) {
            visPoint[answer_[answerCount + j]] = false;
        }
#ifdef COUTTIME
        tans2 = std::chrono::steady_clock::now();
        combinset += std::chrono::duration<double>(tans2 - tans1).count();
        tans11 = std::chrono::steady_clock::now();
#endif

        for (ui j = 0; j < addPointNeighbor; j++) {
            VisNeighbor[NeighborSet[NeighborCount + j]] = false;
        }

#ifdef COUTTIME
        tans22 = std::chrono::steady_clock::now();
        combinnei += std::chrono::duration<double>(tans2 - tans1).count();
#endif
    }
}
uint64_t Graph::queryKBlackHoleParallel(const ui k, double y, int eachin,
                                        const std::string answer_path) {

    std::ofstream ofs(answer_path);
    uint64_t totalKSizeCount = 0;

    VertexID **localAnswer = new VertexID *[threads_count_];
    VertexID **localNeighborSet = new VertexID *[threads_count_];
    bool **localVisPoint = new bool *[threads_count_]();
    bool **localVisNeighbor = new bool *[threads_count_]();
    VertexID **localUpPoint = new VertexID *[threads_count_];
    for (ui i = 0; i < threads_count_; i++) {
        localAnswer[i] = new VertexID[vertices_count_];
        localNeighborSet[i] = new VertexID[scc_count_];
        localVisPoint[i] = new bool[vertices_count_]();
        localVisNeighbor[i] = new bool[scc_count_]();
        localUpPoint[i] = new VertexID[scc_count_];
    }

#pragma omp parallel for schedule(dynamic) reduction(+ : totalKSizeCount)
    for (int i = 0; i < (int)scc_count_; ++i) {
        if (getByPartBHVerticesDegree(i) == 1 ||
            getByPartBHVerticesDegree(i) > k)
            continue;

        BloomFilter localDistinct(1000000, 0.01);
        int thread_id = omp_get_thread_num();
        if (thread_id > threads_count_) {
            cout << "nowThread_id:" << thread_id
                 << " maxThreads_count_:" << threads_count_ << "\n";
            cout << "Please use setThreads_count_\n";
            exit(0);
        }
        ui upPointCount = 0;
        ui answer_cnt_ = 0;
        ui NeighborCount = 0;
        uint64_t localKSizeCount = 0;

        localUpPoint[thread_id][upPointCount++] = i;
        localVisNeighbor[thread_id][i] = true;

        ui count = 0;
        ui countNei = 0;

        const ui *neighbors = getByPartBHVerticesNeighbors(i, count);
        for (ui j = 0; j < count; j++) {
            localAnswer[thread_id][answer_cnt_++] = neighbors[j];
            localVisPoint[thread_id][neighbors[j]] = true;
        }

        const ui *neighborsN = getPartBlackHoleConnectNeighbors(i, countNei);

        for (ui j = 0; j < countNei; j++) {
            if (neighborsN[j] <= i)
                continue;
            localNeighborSet[thread_id][NeighborCount++] = neighborsN[j];
            localVisNeighbor[thread_id][neighborsN[j]] = true;
        }

        if (answer_cnt_ <= k) {
            searchAnswer(0, i, k, answer_cnt_, localAnswer[thread_id],
                         localVisPoint[thread_id], NeighborCount,
                         localNeighborSet[thread_id],
                         localVisNeighbor[thread_id], upPointCount,
                         localUpPoint[thread_id], y, eachin, ofs,
                         localKSizeCount, localDistinct);
        }
        for (ui j = 0; j < count; j++) {
            localVisPoint[thread_id][neighbors[j]] = false;
        }
        for (ui j = 0; j < countNei; j++) {
            if (neighborsN[j] <= i)
                continue;
            localVisNeighbor[thread_id][neighborsN[j]] = false;
        }
        totalKSizeCount += localKSizeCount;
    }

    ofs.close();
    for (ui i = 0; i < threads_count_; i++) {
        delete[] localAnswer[i];
        delete[] localNeighborSet[i];
        delete[] localVisPoint[i];
        delete[] localVisNeighbor[i];
        delete[] localUpPoint[i];
    }
    delete[] localAnswer;
    delete[] localNeighborSet;
    delete[] localVisPoint;
    delete[] localVisNeighbor;
    delete[] localUpPoint;
    return totalKSizeCount;
}

uint64_t Graph::queryKBlackHole(const ui k, double y, int eachin,
                                const std::string answer_path) {
#ifdef COUTTIME
    combinnei = combinset = TimeSelectAnswer = allreadyTime = 0;
#endif

    VertexID *answer_ = new VertexID[vertices_count_];
    VertexID *NeighborSet = new VertexID[scc_count_];
    bool *visPoint = new bool[vertices_count_]();
    bool *VisNeighbor = new bool[scc_count_]();
    VertexID *upPoint = new VertexID[scc_count_];

    ui upPointCount = 0;
    ui answer_cnt_ = 0;
    ui NeighborCount = 0;
    auto tansall1 = std::chrono::steady_clock::now();

#ifdef COUTTIME
    auto tans1ex = std::chrono::steady_clock::now();
#endif

    BloomFilter globalDistinct(1000000, 0.01);

#ifdef COUTTIME
    auto tans2ex = std::chrono::steady_clock::now();
    allreadyTime += std::chrono::duration<double>(tans2ex - tans1ex).count();
#endif

    std::ofstream ofs(answer_path);
    globalDistinct.clear();
    uint64_t kSizeCount = 0;
    for (ui i = 0; i < scc_count_; i++) {
        if (getByPartBHVerticesDegree(i) == 1 ||
            getByPartBHVerticesDegree(i) > k)
            continue;
        upPoint[upPointCount++] = i;
        VisNeighbor[i] = true;

#ifdef COUTTIME
        tans1ex = std::chrono::steady_clock::now();
#endif

        ui count = 0;

        const ui *neighbors = getByPartBHVerticesNeighbors(i, count);
        for (ui j = 0; j < count; j++) {
            answer_[answer_cnt_++] = neighbors[j];
            visPoint[neighbors[j]] = true;
        }

#ifdef COUTTIME
        tans2ex = std::chrono::steady_clock::now();
        combinset += std::chrono::duration<double>(tans2ex - tans1ex).count();
#endif
        ui countNei = 0;
        const ui *neighborsN = getPartBlackHoleConnectNeighbors(i, countNei);

#ifdef COUTTIME
        tans1ex = std::chrono::steady_clock::now();
#endif

        for (ui j = 0; j < countNei; j++) {
            if (neighborsN[j] <= i)
                continue;
            NeighborSet[NeighborCount++] = neighborsN[j];
            VisNeighbor[neighborsN[j]] = true;
        }
#ifdef COUTTIME
        tans2ex = std::chrono::steady_clock::now();
        combinnei += std::chrono::duration<double>(tans2ex - tans1ex).count();
#endif

        searchAnswer(0, i, k, answer_cnt_, answer_, visPoint, NeighborCount,
                     NeighborSet, VisNeighbor, upPointCount, upPoint, y, eachin,
                     ofs, kSizeCount, globalDistinct);

        VisNeighbor[i] = false;
        upPointCount = 0;
#ifdef COUTTIME
        tans1ex = std::chrono::steady_clock::now();
#endif

        for (ui j = 0; j < count; j++) {
            visPoint[neighbors[j]] = false;
        }
#ifdef COUTTIME
        tans2ex = std::chrono::steady_clock::now();
        combinset += std::chrono::duration<double>(tans2ex - tans1ex).count();
#endif

        answer_cnt_ = 0;
        NeighborCount = 0;
#ifdef COUTTIME
        tans1ex = std::chrono::steady_clock::now();
#endif

        for (ui j = 0; j < countNei; j++) {
            VisNeighbor[neighborsN[j]] = false;
        }
#ifdef COUTTIME
        tans2ex = std::chrono::steady_clock::now();
        combinnei += std::chrono::duration<double>(tans2ex - tans1ex).count();
#endif
    }
    ofs.close();
    auto tansall2 = std::chrono::steady_clock::now();
    double TimeAll = std::chrono::duration<double>(tansall2 - tansall1).count();
    delete[] answer_;
    delete[] NeighborSet;
    delete[] visPoint;
    delete[] VisNeighbor;
    delete[] upPoint;
#ifdef COUTTIME
    cout << "nei:" << combinnei << "s   set:" << combinset
         << "s, ansTime:" << TimeSelectAnswer << "s  exitTime:" << allreadyTime
         << "s TimeALL:" << TimeAll << "\n";
#endif

    return kSizeCount;
}

void Graph::saveIndex(std::string dataid, int max) {
    std::string filename =
        "../index/" + dataid + "_" + std::to_string(max) + "_index.bin";
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "Cannot open file for saving index: " << filename << "\n";
        return;
    }

    // 基本信息
    ofs.write((char *)&vertices_count_, sizeof(ui));
    ofs.write((char *)&edges_count_, sizeof(ui));
    ofs.write((char *)&scc_count_, sizeof(ui));

    // 原图正向边
    ofs.write((char *)offsets_, sizeof(ul) * (vertices_count_ + 1));
    ofs.write((char *)neighbors_, sizeof(VertexID) * edges_count_);

    // 原图反向边
    ofs.write((char *)reverse_offsets_, sizeof(ul) * (vertices_count_ + 1));
    ofs.write((char *)reverse_neighbors_, sizeof(VertexID) * edges_count_);

    // labels
    ofs.write((char *)labels_, sizeof(LabelID) * vertices_count_);

    // scc图
    ofs.write((char *)scc_id_, sizeof(ui) * vertices_count_);

    ul scc_edge_total = scc_offsets_[scc_count_];
    ofs.write((char *)scc_offsets_, sizeof(ul) * (scc_count_ + 1));
    ofs.write((char *)scc_neighbors_, sizeof(VertexID) * scc_edge_total);

    ul scc_rev_edge_total = scc_rev_offsets_[scc_count_];
    ofs.write((char *)scc_rev_offsets_, sizeof(ul) * (scc_count_ + 1));
    ofs.write((char *)scc_rev_neighbors_,
              sizeof(VertexID) * scc_rev_edge_total);

    ul scc_v_total = scc_vertices_offsets_[scc_count_];
    ofs.write((char *)scc_vertices_offsets_, sizeof(ul) * (scc_count_ + 1));
    ofs.write((char *)scc_vertices_, sizeof(VertexID) * scc_v_total);

    // partBlackHole（最小黑洞点集）
    ul bh_total = partBlackHole_offsets_[scc_count_];
    ofs.write((char *)partBlackHole_offsets_, sizeof(ul) * (scc_count_ + 1));
    ofs.write((char *)partBlackHole_minIndeg_, sizeof(ui) * scc_count_);
    ofs.write((char *)partBlackHole_neighbors_, sizeof(VertexID) * bh_total);

    // partBlackHoleConnect（黑洞间连接）
    ul conn_total = partBlackHoleConnect_offsets_[scc_count_];
    ofs.write((char *)partBlackHoleConnect_offsets_,
              sizeof(ul) * (scc_count_ + 1));
    ofs.write((char *)partBlackHoleConnect_neighbors_,
              sizeof(VertexID) * conn_total);

    ofs.close();
    std::cout << "saveIndex done: " << filename << "\n";
}

bool Graph::loadIndex(std::string dataid, int max) {
    std::string filename =
        "../index/" + dataid + "_" + std::to_string(max) + "_index.bin";
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.is_open()) {
        std::cerr << "Cannot open index file: " << filename << "\n";
        return false;
    }

    ifs.read((char *)&vertices_count_, sizeof(ui));
    ifs.read((char *)&edges_count_, sizeof(ui));
    ifs.read((char *)&scc_count_, sizeof(ui));

    // 原图正向边
    offsets_ = new ul[vertices_count_ + 1];
    neighbors_ = new VertexID[edges_count_];
    ifs.read((char *)offsets_, sizeof(ul) * (vertices_count_ + 1));
    ifs.read((char *)neighbors_, sizeof(VertexID) * edges_count_);

    // 原图反向边
    reverse_offsets_ = new ul[vertices_count_ + 1];
    reverse_neighbors_ = new VertexID[edges_count_];
    ifs.read((char *)reverse_offsets_, sizeof(ul) * (vertices_count_ + 1));
    ifs.read((char *)reverse_neighbors_, sizeof(VertexID) * edges_count_);

    // labels
    labels_ = new LabelID[vertices_count_];
    ifs.read((char *)labels_, sizeof(LabelID) * vertices_count_);

    // scc图
    scc_id_ = new ui[vertices_count_];
    ifs.read((char *)scc_id_, sizeof(ui) * vertices_count_);

    scc_offsets_ = new ul[scc_count_ + 1];
    ifs.read((char *)scc_offsets_, sizeof(ul) * (scc_count_ + 1));
    ul scc_edge_total = scc_offsets_[scc_count_];
    scc_neighbors_ = new VertexID[scc_edge_total];
    ifs.read((char *)scc_neighbors_, sizeof(VertexID) * scc_edge_total);

    scc_rev_offsets_ = new ul[scc_count_ + 1];
    ifs.read((char *)scc_rev_offsets_, sizeof(ul) * (scc_count_ + 1));
    ul scc_rev_edge_total = scc_rev_offsets_[scc_count_];
    scc_rev_neighbors_ = new VertexID[scc_rev_edge_total];
    ifs.read((char *)scc_rev_neighbors_, sizeof(VertexID) * scc_rev_edge_total);

    scc_vertices_offsets_ = new ul[scc_count_ + 1];
    ifs.read((char *)scc_vertices_offsets_, sizeof(ul) * (scc_count_ + 1));
    ul scc_v_total = scc_vertices_offsets_[scc_count_];
    scc_vertices_ = new VertexID[scc_v_total];
    ifs.read((char *)scc_vertices_, sizeof(VertexID) * scc_v_total);

    // partBlackHole（最小黑洞点集）
    partBlackHole_offsets_ = new ul[scc_count_ + 1];
    ifs.read((char *)partBlackHole_offsets_, sizeof(ul) * (scc_count_ + 1));
    partBlackHole_minIndeg_ = new ui[scc_count_];
    ifs.read((char *)partBlackHole_minIndeg_, sizeof(ui) * scc_count_);
    ul bh_total = partBlackHole_offsets_[scc_count_];
    partBlackHole_neighbors_ = new VertexID[bh_total];
    ifs.read((char *)partBlackHole_neighbors_, sizeof(VertexID) * bh_total);

    // partBlackHoleConnect（黑洞间连接）
    partBlackHoleConnect_offsets_ = new ul[scc_count_ + 1];
    ifs.read((char *)partBlackHoleConnect_offsets_,
             sizeof(ul) * (scc_count_ + 1));
    ul conn_total = partBlackHoleConnect_offsets_[scc_count_];
    partBlackHoleConnect_neighbors_ = new VertexID[conn_total];
    ifs.read((char *)partBlackHoleConnect_neighbors_,
             sizeof(VertexID) * conn_total);

    ifs.close();
    std::cout << "loadIndex done: " << filename << "\n";
    return true;
}




void Graph::baselineDfs(ui level, const ui k, double y, int eachin, int root,
                        VertexID *NeighborSet, ui neighborCount, int *visPoint,
                        int *visAnswer, VertexID *answer_, ui answerCount,
                        uint64_t &kSizeCount, std::string answer_path) {
    if (answerCount == k) {
        ui sumIndeg = 0;
        ui minIndeg = 0x3f3f3f3f;
        for (int i = 0; i < answerCount; i++) {
            ui count;
            const ui *neighbor = getVertexNeighbors(answer_[i], count);
            for (int j = 0; j < count; j++) {
                ui nowPoint = neighbor[j];
                if (visAnswer[nowPoint] == root)
                    continue;
                return; // 保证每个子图没有出边。
            }

            neighbor = getVertexReverseNeighbors(answer_[i], count);
            ui inDeg = 0;
            for (int j = 0; j < count; j++) {
                ui nowPoint = neighbor[j];
                if (visAnswer[nowPoint] == root) // 查看有多少入边。
                    continue;
                inDeg++;
            }
            minIndeg = min(minIndeg, inDeg);
            sumIndeg += inDeg;
        }
        if (minIndeg < eachin)
            return;                          // 满足每个点至少有eachin个入边。
        double yinli = (double)sumIndeg / k; // 求平均每个点有多少入边
        if (yinli < y)
            return;
        kSizeCount++;
        return;
    }
    for (ui i = level; i < neighborCount; i++) {
        ui nowPoint = NeighborSet[i];
        answer_[answerCount] = nowPoint;
        visAnswer[nowPoint] = root;
        ui count;
        const ui *neighbor = getVertexNeighbors(nowPoint, count);
        ui addNeighbor = 0;
        for (ui j = 0; j < count; j++) {
            ui u = neighbor[j];
            if (visPoint[u] == root || u <= root)
                continue;
            visPoint[u] = root;
            NeighborSet[neighborCount + addNeighbor] = u;
            addNeighbor++;
        }
        neighbor = getVertexReverseNeighbors(nowPoint, count);
        for (ui j = 0; j < count; j++) {
            ui u = neighbor[j];
            if (visPoint[u] == root || u <= root)
                continue;
            visPoint[u] = root;
            NeighborSet[neighborCount + addNeighbor] = u;
            addNeighbor++;
        }
        baselineDfs(i + 1, k, y, eachin, root, NeighborSet,
                    neighborCount + addNeighbor, visPoint, visAnswer, answer_,
                    answerCount + 1, kSizeCount, answer_path);

        for (ui j = 0; j < addNeighbor; j++) {
            visPoint[NeighborSet[neighborCount + j]] = -1;
        }
        visAnswer[nowPoint] = -1;
    }
}


uint64_t Graph::baseline(const ui k, double y, int eachin,
                         const std::string answer_path) {
    VertexID *answer_ = new VertexID[vertices_count_];
    VertexID *NeighborSet = new VertexID[vertices_count_];
    int *visPoint = new int[vertices_count_];
    int *visAnswer = new int[vertices_count_];

    for (ui i = 0; i < vertices_count_; i++) {
        visPoint[i] = -1;
        visAnswer[i] = -1;
    }
    uint64_t kSizeCount = 0;
    for (ui i = 0; i < vertices_count_; i++) {
        ui neighborCount = 0;
        ui count = 0;
        const ui *neighbor = getVertexNeighbors(i, count);
        visPoint[i] = i;
        visAnswer[i] = i;
        for (int j = 0; j < count; j++) {
            if (neighbor[j] > i)
                NeighborSet[neighborCount++] = neighbor[j];
        }
        neighbor = getVertexReverseNeighbors(i, count);
        for (int j = 0; j < count; j++) {
            if (neighbor[j] > i)
                NeighborSet[neighborCount++] = neighbor[j];
        }
        answer_[0] = i;
        baselineDfs(0, k, y, eachin, i, NeighborSet, neighborCount, visPoint,
                    visAnswer, answer_, 1, kSizeCount, answer_path);
    }
    return kSizeCount;
}


void Graph::baselineDfs2(ui level, ui *answer, ui answerCount, ul &KsizeCount,
                        ui k, ui *ind, ui *selectset, ui setSize,
                        BloomFilter &globalDistinct) {

    if (level == setSize)
        return;
    if (answerCount == k) {

        // for (ui i = 0; i < answerCount; i++) {
        //     cout << answer[i] << " ";
        // }
        // cout << "\n";

        UnionFind connect_find(k);
        bool valid = true;
        for (ui i = 0; i < k; ++i) {
            ui count;
            const ui *neighbors = getVertexNeighbors(answer[i], count);
            for (ui j = 0; j < count; ++j) {
                ui neighbor = neighbors[j];
                if (ind[neighbor] != -1) {
                    connect_find.unite(i, ind[neighbor]);
                } else {
                    valid = false;
                    break;
                }
            }
            if (!valid) {
                break;
            }
        }
        ui cou = connect_find.countRoots();
        bool alreadyExists =
            globalDistinct.checkAndInsertArray(answer, answerCount);
        if (cou == 1 && valid == true && !alreadyExists) {
            // KsizeCount++;
        }
        return;
    }

    VertexID v = selectset[level];
    answer[answerCount] = v;
    ind[v] = answerCount;
    baselineDfs2(level + 1, answer, answerCount + 1, KsizeCount, k, ind,
                selectset, setSize, globalDistinct);
    ind[v] = -1;

    baselineDfs2(level + 1, answer, answerCount, KsizeCount, k, ind, selectset,
                setSize, globalDistinct);
}
uint64_t Graph::baseline2(const ui k, string answer_path) {
    ui *answer = new ui[k];
    ui *ind = new ui[vertices_count_];
    for (ui i = 0; i < vertices_count_; i++)
        ind[i] = -1;
    ul Kisze = 0;
    ui nowTime = 0;
    ui *vis = new ui[vertices_count_]();
    ui *selectSet = new ui[vertices_count_];
    ui setSize = 0;
    BloomFilter globalDistinct(1000000, 0.01);
    for (ui i = 0; i < vertices_count_; i++) {
        globalDistinct.clear();
        queue<pair<ui, ui>> que;
        nowTime++;
        que.push({i, 0});
        setSize = 1;
        selectSet[0] = i;
        vis[i] = nowTime;
        while (!que.empty()) {
            ui step = que.front().second;
            ui id = que.front().first;
            que.pop();
            ui count;
            const ui *neighbors = getVertexNeighbors(i, count);
            for (ui j = 0; j < count; j++) {
                VertexID v = neighbors[j];
                if (vis[v] != nowTime && step < k - 1) {
                    vis[v] = nowTime;
                    selectSet[setSize] = v;
                    setSize++;
                    que.push({v, step + 1});
                }
            }
            neighbors = getVertexReverseNeighbors(i, count);
            for (ui j = 0; j < count; j++) {
                VertexID v = neighbors[j];
                if (vis[v] != nowTime && step < k - 1) {
                    vis[v] = nowTime;
                    selectSet[setSize] = v;
                    setSize++;
                    que.push({v, step + 1});
                }
            }
        }
        baselineDfs2(0, answer, 0, Kisze, k, ind, selectSet, setSize,
                    globalDistinct);
    }

    delete[] answer;
    delete[] ind;
    return Kisze;
}
