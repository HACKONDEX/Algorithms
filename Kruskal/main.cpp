#include <algorithm>
#include <iostream>
#include <vector>

class DisjointSetUnion {
  public:
    DisjointSetUnion(int n) : vertex_count(n) {
        for (int i = 0; i < n; ++i)
            roots.emplace_back(i, 1);
    }

    ~DisjointSetUnion() = default;

    bool AreInSameComponent(int v, int u);
    void Join(int u, int v);
    int GetVertexCount() const;

  private:
    int GetRoot(int v);

  private:
    std::vector<std::pair<int, int> > roots;
    int vertex_count;

};

int DisjointSetUnion::GetVertexCount() const {
    return vertex_count;
}

bool DisjointSetUnion::AreInSameComponent(int v, int u) {
    return GetRoot(v) == GetRoot(u);
}

int DisjointSetUnion::GetRoot(int v) {
    if (v == roots[v].first) {
        return v;
    }
    int current_root = GetRoot(roots[v].first);
    roots[v].first = current_root;
    return current_root;
}

void DisjointSetUnion::Join(int u, int v) {
    int head = GetRoot(u), tail = GetRoot(v);
    if (roots[head].second >= roots[tail].second) {
        std::swap(head, tail);
    }
    roots[head].second += roots[tail].second;
    roots[tail].first = head;
}

struct Edge {
    int left;
    int right;
    int weight;

    explicit Edge(int l = 0, int r = 0, int w = 0) : left(l), right(r), weight(w) {}

    bool operator<(const Edge &other) const;
};

bool Edge::operator<(const Edge &other) const {
    if (weight < other.weight) {
        return true;
    } else if (weight == other.weight) {
        if (left < other.left) {
            return true;
        } else if (left == other.left) {
            return right < other.right;
        }
        return false;
    }
    return false;
}

bool Comparator(const Edge &first, const Edge &second) {
    return first < second;
}

void
GetSortedEdges(std::vector<Edge> &edges, int m, std::vector<int> &from, std::vector<int> &to, std::vector<int> &cost) {
    for (int i = 0; i < m; ++i) {
        edges.emplace_back(from[i], to[i], cost[i]);
    }
    std::sort(edges.begin(), edges.end(), Comparator);
}

long long Kruskal(std::vector<Edge> &edges, DisjointSetUnion &dsu) {
    long long minimum_weight_sum = 0;
    int n = dsu.GetVertexCount(), u, v, marked = 0;
    std::vector<bool> colour(n, false);
    for (int i = 0; i < edges.size(); ++i) {
        u = edges[i].left;
        v = edges[i].right;
        if (dsu.AreInSameComponent(v, u)) {
            continue;
        }
        minimum_weight_sum += edges[i].weight;
        dsu.Join(u, v);
        marked += (!colour[u] ? 1 : 0) + (!colour[v] ? 1 : 0);
        colour[u] = colour[v] = true;
        if (marked >= n - 1) {
            break;
        }
    }
    return minimum_weight_sum;
}

long long FindMSTWeight(int n, int m, std::vector<int> &from, std::vector<int> &to, std::vector<int> &cost) {
    std::vector<Edge> edges;
    GetSortedEdges(edges, m, from, to, cost);
    DisjointSetUnion dsu(n + 1);
    return Kruskal(edges, dsu);
}

int main() {
    int n = 0, m = 0;
    std::cin >> n >> m;
    std::vector<int> from(m);
    std::vector<int> to(m);
    std::vector<int> cost(m);
    for (int i = 0; i < m; ++i) {
        std::cin >> from[i] >> to[i] >> cost[i];
    }

    std::cout << FindMSTWeight(n, m, from, to, cost) << std::endl;
    return 0;
}
