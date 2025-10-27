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
ul countiCount = 0, useCount = 0;

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
            labels_vertices_[label] = id;
        } else if (type == 'e') { // Read edge.
            VertexID begin;
            VertexID end;
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

    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(neighbors_ + offsets_[i], neighbors_ + offsets_[i + 1]);
    }
    for (ui i = 0; i < vertices_count_; ++i) {
        std::sort(reverse_neighbors_ + reverse_offsets_[i],
                  reverse_neighbors_ + reverse_offsets_[i + 1]);
    }
}

void Graph::buildSCCGraph() {
    cout << "buildSCC\n";
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

    vector<vector<ui>> answer;
    vector<ui> vis(scc_count_);
    ui nowtime = 0;
    for (ui i = 0; i < scc_count_; i++) { // 每个点
        if (i % (scc_count_ / 10) == 0) {
            std::cout << "create minimum over:" << (double)i / scc_count_
                      << "\n";
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

        sort(nowans.begin(), nowans.end());
        answer.push_back(nowans);
    }

    std::sort(answer.begin(), answer.end(),
              [](const std::vector<ui> &A, const std::vector<ui> &B) {
                  return A.size() > B.size();
              });
    scc_count_ = answer.size();

    partBlackHole_offsets_ = new ul[scc_count_ + 1];
    partBlackHole_offsets_[0] = 0;

    for (ui i = 1; i <= scc_count_; i++) {
        partBlackHole_offsets_[i] =
            partBlackHole_offsets_[i - 1] + answer[i - 1].size();
    }
    partBlackHole_neighbors_ = new VertexID[partBlackHole_offsets_[scc_count_]];

    for (ui i = 0; i < scc_count_; i++) {
        ui count = answer[i].size();
        ul start = partBlackHole_offsets_[i];
        for (ui j = 0; j < count; j++) {
            partBlackHole_neighbors_[partBlackHole_offsets_[i] + j] =
                answer[i][j];
        }
    }

    cout << partBlackHole_offsets_[scc_count_] << " size\n";
    cout << "getTheMinBlackHole_over\n";
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
            std::cout << "connect SCC over:" << (double)i / scc_count_ << "\n";
        }
        if (scc_neighbors[i].size() > maxset || scc_neighbors[i].size() == 0)
            continue;
        if (scc_neighbors[i].size() == 1)
            break;

        for (ui j = i + 1; j < scc_count_; j++) {

            if (scc_neighbors[j].size() == 1)
                continue;

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

    cout << "edge:" << cu << "\n";
    cout << "connectSCC_over\n";
}
void Graph::buildIndex(int maxset) {
    buildSCCGraph();
    getTheMinBlackHoleBFS(maxset);
    connectSCC(maxset);
}

void Graph::searchAnswer(ui level, ui start, const ui k, ui answerCount,
                         VertexID *answer_, bool *visPoint, ui NeighborCount,
                         VertexID *NeighborSet, bool *VisNeighbor,
                         ui upPointCount, VertexID *upPoint, std::ofstream &ofs,
                         uint64_t &kSizeCount, BloomFilter &localDistinct) {

    if (answerCount == k) {

        bool alreadyExists =
            localDistinct.checkAndInsertArray(answer_, answerCount);
        if (!alreadyExists) {
            kSizeCount++;
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

        const ui *neighbor = getByPartBHVerticesNeighbors(u, count);
        for (ui j = 0; j < count; j++) {
            VertexID v = neighbor[j];
            if (visPoint[v])
                continue;
            visPoint[v] = true;

            answer_[answerCount + addPoint] = v;
            addPoint++;
        }

        neighbor = getPartBlackHoleConnectNeighbors(u, count);

        for (ui j = 0; j < count; j++) {
            VertexID v = neighbor[j];
            if (VisNeighbor[v] || v <= start)
                continue;
            VisNeighbor[v] = true;
            NeighborSet[NeighborCount + addPointNeighbor] = v;
            addPointNeighbor++;
        }
        upPoint[upPointCount] = u;

        if (addPoint != 0)
            searchAnswer(i + 1, start, k, answerCount + addPoint, answer_,
                         visPoint, NeighborCount + addPointNeighbor,
                         NeighborSet, VisNeighbor, upPointCount + 1, upPoint,
                         ofs, kSizeCount, localDistinct);

        for (ui j = 0; j < addPoint; j++) {
            visPoint[answer_[answerCount + j]] = false;
        }

        for (ui j = 0; j < addPointNeighbor; j++) {
            VisNeighbor[NeighborSet[NeighborCount + j]] = false;
        }
    }
}

uint64_t Graph::queryKBlackHoleParallel(const ui k,
                                        const std::string answer_path) {

    std::ofstream ofs(answer_path);
    uint64_t totalKSizeCount = 0;
    countiCount = 0, useCount = 0;

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

        if (answer_cnt_ != k) {
            searchAnswer(0, i, k, answer_cnt_, localAnswer[thread_id],
                         localVisPoint[thread_id], NeighborCount,
                         localNeighborSet[thread_id],
                         localVisNeighbor[thread_id], upPointCount,
                         localUpPoint[thread_id], ofs, localKSizeCount,
                         localDistinct);

        } else {
            localKSizeCount++;
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
    return totalKSizeCount;
}

uint64_t Graph::queryKBlackHole(const ui k, const std::string answer_path) {

    VertexID *extendFather = new VertexID[scc_count_];
    VertexID *answer_ = new VertexID[vertices_count_];
    VertexID *NeighborSet = new VertexID[scc_count_];
    bool *visPoint = new bool[vertices_count_]();
    bool *VisNeighbor = new bool[scc_count_]();
    VertexID *upPoint = new VertexID[scc_count_];

    ui upPointCount = 0;
    ui answer_cnt_ = 0;
    ui NeighborCount = 0;

    BloomFilter globalDistinct(1000000, 0.01);

    std::ofstream ofs(answer_path);
    uint64_t kSizeCount = 0;
    countiCount = 0, useCount = 0;
    for (ui i = 0; i < scc_count_; i++) {
        if (getByPartBHVerticesDegree(i) == 1 ||
            getByPartBHVerticesDegree(i) > k)
            continue;
        globalDistinct.clear();
        upPoint[upPointCount++] = i;
        VisNeighbor[i] = true;

        ui count = 0;
        ui countNei = 0;

        const ui *neighbors = getByPartBHVerticesNeighbors(i, count);
        for (ui j = 0; j < count; j++) {
            answer_[answer_cnt_++] = neighbors[j];
            visPoint[neighbors[j]] = true;
        }

        const ui *neighborsN = getPartBlackHoleConnectNeighbors(i, countNei);

        for (ui j = 0; j < countNei; j++) {
            if (neighborsN[j] <= i)
                continue;
            extendFather[NeighborCount] = i;
            NeighborSet[NeighborCount++] = neighborsN[j];

            VisNeighbor[neighborsN[j]] = true;
        }

        if (answer_cnt_ != k) {
            searchAnswer(0, i, k, answer_cnt_, answer_, visPoint, NeighborCount,
                         NeighborSet, VisNeighbor, upPointCount, upPoint, ofs,
                         kSizeCount, globalDistinct);
        } else {
            kSizeCount++;
        }

        VisNeighbor[i] = false;
        upPointCount = 0;
        for (ui j = 0; j < count; j++) {
            visPoint[neighbors[j]] = false;
        }

        answer_cnt_ = 0;
        NeighborCount = 0;

        for (ui j = 0; j < countNei; j++) {
            VisNeighbor[neighborsN[j]] = false;
        }
    }
    ofs.close();

    delete[] answer_;
    delete[] NeighborSet;
    delete[] visPoint;
    delete[] VisNeighbor;
    delete[] upPoint;

    return kSizeCount;
}

void Graph::baselineDfs(ui level, ui *answer, ui answerCount, ul &KsizeCount,
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
    baselineDfs(level + 1, answer, answerCount + 1, KsizeCount, k, ind,
                selectset, setSize, globalDistinct);
    ind[v] = -1;

    baselineDfs(level + 1, answer, answerCount, KsizeCount, k, ind, selectset,
                setSize, globalDistinct);
}
uint64_t Graph::baseline(const ui k, string answer_path) {
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
        baselineDfs(0, answer, 0, Kisze, k, ind, selectSet, setSize,
                    globalDistinct);
    }

    delete[] answer;
    delete[] ind;
    return Kisze;
}

void Graph::saveIndex(std::string dataid) {
    std::string filename = "../index/" + dataid + "_index";
    std::ofstream out(filename, std::ios::binary);

    if (!out.is_open()) {
        std::cout << "Cannot create index file: " << filename << std::endl;
        return;
    }

    out.write(reinterpret_cast<const char *>(&vertices_count_),
              sizeof(vertices_count_));
    out.write(reinterpret_cast<const char *>(&edges_count_),
              sizeof(edges_count_));

    if (labels_ != nullptr) {
        out.write(reinterpret_cast<const char *>(labels_),
                  sizeof(LabelID) * vertices_count_);
    }

    if (offsets_ != nullptr) {
        out.write(reinterpret_cast<const char *>(offsets_),
                  sizeof(ul) * (vertices_count_ + 1));
    }
    if (neighbors_ != nullptr) {
        out.write(reinterpret_cast<const char *>(neighbors_),
                  sizeof(VertexID) * edges_count_);
    }

    if (reverse_offsets_ != nullptr) {
        out.write(reinterpret_cast<const char *>(reverse_offsets_),
                  sizeof(ul) * (vertices_count_ + 1));
    }
    if (reverse_neighbors_ != nullptr) {
        out.write(reinterpret_cast<const char *>(reverse_neighbors_),
                  sizeof(VertexID) * edges_count_);
    }

    out.write(reinterpret_cast<const char *>(&scc_count_), sizeof(scc_count_));

    if (scc_count_ > 0) {
        if (partBlackHole_offsets_ != nullptr) {
            out.write(reinterpret_cast<const char *>(partBlackHole_offsets_),
                      sizeof(ul) * (scc_count_ + 1));
            out.write(reinterpret_cast<const char *>(partBlackHole_neighbors_),
                      sizeof(VertexID) * partBlackHole_offsets_[scc_count_]);
        }

        if (partBlackHoleConnect_offsets_ != nullptr) {
            out.write(
                reinterpret_cast<const char *>(partBlackHoleConnect_offsets_),
                sizeof(ul) * (scc_count_ + 1));
            out.write(
                reinterpret_cast<const char *>(partBlackHoleConnect_neighbors_),
                sizeof(VertexID) * partBlackHoleConnect_offsets_[scc_count_]);
        }
    }

    size_t labels_vertices_size = labels_vertices_.size();
    out.write(reinterpret_cast<const char *>(&labels_vertices_size),
              sizeof(labels_vertices_size));
    for (const auto &pair : labels_vertices_) {
        LabelID label = pair.first;
        ui vertex_id = pair.second;
        out.write(reinterpret_cast<const char *>(&label), sizeof(label));
        out.write(reinterpret_cast<const char *>(&vertex_id),
                  sizeof(vertex_id));
    }

    if (scc_offsets_ != nullptr) {
        out.write(reinterpret_cast<const char *>(scc_offsets_),
                  sizeof(ul) * (scc_count_ + 1));
        out.write(reinterpret_cast<const char *>(scc_neighbors_),
                  sizeof(VertexID) * scc_offsets_[scc_count_]);
    }

    if (scc_id_ != nullptr) {
        out.write(reinterpret_cast<const char *>(scc_id_),
                  sizeof(ui) * vertices_count_);
    }

    if (scc_vertices_offsets_ != nullptr) {
        out.write(reinterpret_cast<const char *>(scc_vertices_offsets_),
                  sizeof(ul) * (scc_count_ + 1));
        out.write(reinterpret_cast<const char *>(scc_vertices_),
                  sizeof(VertexID) * scc_vertices_offsets_[scc_count_]);
    }

    if (scc_rev_offsets_ != nullptr) {
        out.write(reinterpret_cast<const char *>(scc_rev_offsets_),
                  sizeof(ul) * (scc_count_ + 1));
        out.write(reinterpret_cast<const char *>(scc_rev_neighbors_),
                  sizeof(VertexID) * scc_rev_offsets_[scc_count_]);
    }

    out.close();
    std::cout << "index load: " << filename << std::endl;
}

bool Graph::loadIndex(std::string dataid) {
    std::string filename = "../index/" + dataid + "_index";
    std::ifstream in(filename, std::ios::binary);

    if (!in.is_open()) {
        std::cout << "NO index file: " << filename << std::endl;
        return 0;
    }

    if (!in.read(reinterpret_cast<char *>(&vertices_count_),
                 sizeof(vertices_count_)) ||
        !in.read(reinterpret_cast<char *>(&edges_count_),
                 sizeof(edges_count_))) {
        cout << "read index error " << endl;
        in.close();
        exit(0);
    }

    if (vertices_count_ > 0) {
        labels_ = new LabelID[vertices_count_];
        in.read(reinterpret_cast<char *>(labels_),
                sizeof(LabelID) * vertices_count_);
    }

    if (vertices_count_ > 0) {
        offsets_ = new ul[vertices_count_ + 1];
        in.read(reinterpret_cast<char *>(offsets_),
                sizeof(ul) * (vertices_count_ + 1));
    }

    if (edges_count_ > 0) {
        neighbors_ = new VertexID[edges_count_];
        in.read(reinterpret_cast<char *>(neighbors_),
                sizeof(VertexID) * edges_count_);
    }

    if (vertices_count_ > 0) {
        reverse_offsets_ = new ul[vertices_count_ + 1];
        in.read(reinterpret_cast<char *>(reverse_offsets_),
                sizeof(ul) * (vertices_count_ + 1));
    }

    if (edges_count_ > 0) {
        reverse_neighbors_ = new VertexID[edges_count_];
        in.read(reinterpret_cast<char *>(reverse_neighbors_),
                sizeof(VertexID) * edges_count_);
    }

    if (!in.read(reinterpret_cast<char *>(&scc_count_), sizeof(scc_count_))) {
        cout << "read index error " << endl;
        in.close();
        exit(0);
    }

    // 读取黑洞相关结构
    if (scc_count_ > 0) {
        partBlackHole_offsets_ = new ul[scc_count_ + 1];
        if (!in.read(reinterpret_cast<char *>(partBlackHole_offsets_),
                     sizeof(ul) * (scc_count_ + 1))) {
            cout << "read index error " << endl;
            in.close();
            exit(0);
        }

        partBlackHole_neighbors_ =
            new VertexID[partBlackHole_offsets_[scc_count_]];
        if (!in.read(reinterpret_cast<char *>(partBlackHole_neighbors_),
                     sizeof(VertexID) * partBlackHole_offsets_[scc_count_])) {
            cout << "read index error " << endl;
            in.close();
            exit(0);
        }

        partBlackHoleConnect_offsets_ = new ul[scc_count_ + 1];
        if (!in.read(reinterpret_cast<char *>(partBlackHoleConnect_offsets_),
                     sizeof(ul) * (scc_count_ + 1))) {
            cout << "read index error " << endl;
            in.close();
            exit(0);
        }

        partBlackHoleConnect_neighbors_ =
            new VertexID[partBlackHoleConnect_offsets_[scc_count_]];
        if (!in.read(reinterpret_cast<char *>(partBlackHoleConnect_neighbors_),
                     sizeof(VertexID) *
                         partBlackHoleConnect_offsets_[scc_count_])) {
            cout << "read index error " << endl;
            in.close();
            exit(0);
        }
    }

    size_t labels_vertices_size;
    if (!in.read(reinterpret_cast<char *>(&labels_vertices_size),
                 sizeof(labels_vertices_size))) {
        cout << "read index error " << endl;
        in.close();
        exit(0);
    }
    labels_vertices_.clear();
    for (size_t i = 0; i < labels_vertices_size; i++) {
        LabelID label;
        ui vertex_id;
        if (!in.read(reinterpret_cast<char *>(&label), sizeof(label)) ||
            !in.read(reinterpret_cast<char *>(&vertex_id), sizeof(vertex_id))) {
            cout << "read index error " << endl;
            in.close();
            exit(0);
        }
        labels_vertices_[label] = vertex_id;
    }
    if (scc_count_ > 0) {
        scc_offsets_ = new ul[scc_count_ + 1];
        in.read(reinterpret_cast<char *>(scc_offsets_),
                sizeof(ul) * (scc_count_ + 1));

        scc_neighbors_ = new VertexID[scc_offsets_[scc_count_]];
        in.read(reinterpret_cast<char *>(scc_neighbors_),
                sizeof(VertexID) * scc_offsets_[scc_count_]);

        scc_id_ = new ui[vertices_count_];
        in.read(reinterpret_cast<char *>(scc_id_),
                sizeof(ui) * vertices_count_);

        scc_vertices_offsets_ = new ul[scc_count_ + 1];
        in.read(reinterpret_cast<char *>(scc_vertices_offsets_),
                sizeof(ul) * (scc_count_ + 1));

        scc_vertices_ = new VertexID[scc_vertices_offsets_[scc_count_]];
        in.read(reinterpret_cast<char *>(scc_vertices_),
                sizeof(VertexID) * scc_vertices_offsets_[scc_count_]);

        scc_rev_offsets_ = new ul[scc_count_ + 1];
        in.read(reinterpret_cast<char *>(scc_rev_offsets_),
                sizeof(ul) * (scc_count_ + 1));

        scc_rev_neighbors_ = new VertexID[scc_rev_offsets_[scc_count_]];
        in.read(reinterpret_cast<char *>(scc_rev_neighbors_),
                sizeof(VertexID) * scc_rev_offsets_[scc_count_]);
    }
    in.close();

    std::cout << "Index load over" << endl;
    return 1;
}
