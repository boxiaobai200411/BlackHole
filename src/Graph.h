#ifndef GRAPH_H
#define GRAPH_H

#include "type.h"

#include <string>
#include <unordered_map>
#include <vector>
class Graph {
  private:
    ui vertices_count_; // 点个数
    ui edges_count_;    // 边个数
    ui threads_count_ = 20;
    std::string dataid_;
    LabelID *labels_; // 标签数组 大小vertices_count_

    ul *offsets_;         // 每个点邻居偏移量 大小vertices_count_+1
    VertexID *neighbors_; // 邻居数组 大小neighbors_[vertices_count_]

    ul *partBlackHole_offsets_; // 黑洞的偏移 大小ByPartBHVertices_count_+1
    VertexID *
        partBlackHole_neighbors_; // 每个点的黑洞邻居 上方点连接下方
                                  // 大小partBlackHole_offsets_[ByPartBHVertices_count_]
    ul *partBlackHoleConnect_offsets_;
    VertexID *partBlackHoleConnect_neighbors_; // 最小黑洞之间的连接

    ul *reverse_offsets_;         // 反向边的偏移量  vertices_count_+1
    VertexID *reverse_neighbors_; // 反向邻居 reverse_offsets_[vertices_count_]

    std::unordered_map<LabelID, ui>
        labels_vertices_;       // 一共有大小vertices_count_个对
    BloomFilter distinctAnswer; // 去重不用save和lord

    ui scc_count_; // 缩点后的点数

    ul *scc_offsets_;         // 缩点后邻居偏移量
    VertexID *scc_neighbors_; // 缩点后邻居
    ui *scc_id_;              // 每个原始点所属的SCC编号

    ul *scc_vertices_offsets_; // 每个SCC包含的点的偏移量
    VertexID *scc_vertices_;   // 每个SCC对应的原始点集合

    ul *scc_rev_offsets_;         // 缩点后反向邻居偏移量
    VertexID *scc_rev_neighbors_; // 缩点后反向邻居

  public:
    Graph() {
        vertices_count_ = 0;
        edges_count_ = 0;

        offsets_ = NULL;
        neighbors_ = NULL;
        labels_ = NULL;

        reverse_offsets_ = NULL;
        reverse_neighbors_ = NULL;

        partBlackHole_offsets_ = NULL;
        partBlackHole_neighbors_ = NULL;

        labels_vertices_.clear();
    }
    ~Graph() {
        delete[] offsets_;
        delete[] neighbors_;
        delete[] labels_;

        delete[] reverse_offsets_;
        delete[] reverse_neighbors_;
        delete[] partBlackHole_offsets_;
        delete[] partBlackHole_neighbors_;
    }

  public:
    void loadGraphFromFile(const std::string &file_path);

    uint64_t queryKBlackHoleParallel(const ui k, const std::string answer_path);
    uint64_t queryKBlackHole(const ui k, const std::string answer_path);
    uint64_t baseline(const ui k, const std::string answer_path);
    void searchAnswer(ui level, ui start, const ui k, ui answerCount,
                      VertexID *answer_, bool *visPoint, ui NeighborCount,
                      VertexID *NeighborSet, bool *VisNeighbor, ui upPointCount,
                      VertexID *upPoint, std::ofstream &ofs,
                      uint64_t &kSizeCount, BloomFilter &localDistinct);

    void saveIndex(std::string dataid);
    bool loadIndex(std::string dataid);
    double getTotalMemoryMB();

    void buildIndex(int maxset);
    void buildSCCGraph();
    void getTheMinBlackHole();
    void getTheMinBlackHoleBFS(int maxset);
    void connectSCC(int maxset);
    bool isSCCGraphWeaklyConnected();

    void getVerticePartBlackHole(const VertexID id,
                                 std::vector<VertexID> &PartBlackHole);
    void getVerticesPartBlackHole();
    void baselineDfs(ui level, ui *answer, ui answerCount, ul &KsizeCount, ui k,
                     ui *ind, ui *selectset, ui setSize);

  public:
    void setThreads_count_(ui count) { threads_count_ = count; }
    const ui getVerticesCount() const { return vertices_count_; }
    const ui getEdgesCount() const { return edges_count_; }
    const ul getVertexOutDegree(const VertexID id) const {
        return offsets_[id + 1] - offsets_[id];
    }
    const ui *getVertexNeighbors(const VertexID id, ui &count) const {
        count = offsets_[id + 1] - offsets_[id];
        return neighbors_ + offsets_[id];
    } // 每个点的邻居

    const ul getVertexInDegree(const VertexID id) const {
        return reverse_offsets_[id + 1] - reverse_offsets_[id];
    }
    const ui *getVertexReverseNeighbors(const VertexID id, ui &count) const {
        count = reverse_offsets_[id + 1] - reverse_offsets_[id];
        return reverse_neighbors_ + reverse_offsets_[id];
    } // 反向邻居

    const ul getPartBlackHoleConnectDegree(const VertexID id) const {
        return partBlackHoleConnect_offsets_[id + 1] -
               partBlackHoleConnect_offsets_[id];
    }
    const ui *getPartBlackHoleConnectNeighbors(const VertexID id,
                                               ui &count) const {
        count = partBlackHoleConnect_offsets_[id + 1] -
                partBlackHoleConnect_offsets_[id];
        return partBlackHoleConnect_neighbors_ +
               partBlackHoleConnect_offsets_[id];
    } // 最小黑洞之间的连接

    const ul getByPartBHVerticesDegree(const VertexID id) const {
        return partBlackHole_offsets_[id + 1] - partBlackHole_offsets_[id];
    }
    const ui *getByPartBHVerticesNeighbors(const VertexID id, ui &count) const {
        count = partBlackHole_offsets_[id + 1] - partBlackHole_offsets_[id];
        return partBlackHole_neighbors_ + partBlackHole_offsets_[id]; // 上连下
    } // 每个最小黑洞的点

    const ui getLabelsVertices(const LabelID label) const {
        return labels_vertices_.find(label) == labels_vertices_.end()
                   ? -1
                   : labels_vertices_.at(label);
    }
    const LabelID getVertexLabel(const VertexID id) const {
        return labels_[id];
    }

    // 获取缩点后的某个点的邻居
    const ui *getSCCNeighbors(const ui id, ui &count) const {
        count = scc_offsets_[id + 1] - scc_offsets_[id];
        return scc_neighbors_ + scc_offsets_[id];
    }

    // 获取某个SCC包含的原始点
    const ui *getSCCVertices(const ui id, ui &count) const {
        count = scc_vertices_offsets_[id + 1] - scc_vertices_offsets_[id];
        return scc_vertices_ + scc_vertices_offsets_[id];
    }

    // 获取缩点后的某个点的反向邻居
    const ui *getSCCReverseNeighbors(const ui id, ui &count) const {
        count = scc_rev_offsets_[id + 1] - scc_rev_offsets_[id];
        return scc_rev_neighbors_ + scc_rev_offsets_[id];
    }
};

#endif