#include <iostream>
#include <set>
#include <vector>

template<class T>
class Graph {
  public:
    explicit Graph(int vert_count = 0) : vertex_count(vert_count) {
        for (int i = 0; i < vert_count; ++i)
            graph.emplace_back();
    }
    ~Graph() = default;
    Graph(const Graph &) = default;
    Graph(Graph &&) = default;

    Graph &operator=(const Graph &) = delete;
    Graph &operator=(Graph &&) = delete;

    void AddVertex();
    void AddEdge(int from, T to);
    bool HasEdge(int from, T to) const;
    const std::vector<T> &GetNextVertices(int from) const;
    int GetVertexCount() const;

  private:
    int vertex_count;
    std::vector<std::vector<T> > graph;

};

template<class T>
void Graph<T>::AddVertex() {
    ++vertex_count;
}

template<class T>
void Graph<T>::AddEdge(int from, T to) {
    graph[from].push_back(to);
}

template<class T>
bool Graph<T>::HasEdge(int from, T to) const {
    for (int i = 0; i < graph[from].size(); ++i) {
        if (graph[from][i] == to) {
            return true;
        }
    }
    return false;
}

template<class T>
const std::vector<T> &Graph<T>::GetNextVertices(int from) const {
    return graph[from];
}

template<class T>
int Graph<T>::GetVertexCount() const {
    return vertex_count;
}

Graph<std::pair<int, int>> BuildGraph(int a, int b, int m) {
    Graph<std::pair<int, int>> graph(m + 1);
    for (long long i = 0; i < m; ++i) {
        graph.AddEdge(i, std::make_pair((i + 1) % m, a));
        graph.AddEdge(i, std::make_pair((i * i + 1) % m, b));
    }
    return graph;
}

void Dijkstra(Graph<std::pair<int, int>>& graph, int start, std::vector<int>& dist) {
    std::vector<std::pair<int, int>> next_vertices;
    std::set<std::pair<int, int>> bin_tree;
    int current, to, weight;
    dist[start] = 0;
    bin_tree.insert(std::make_pair(dist[start], start));

    while (!bin_tree.empty()) {
        current = bin_tree.begin()->second;
        bin_tree.erase(bin_tree.begin());

        next_vertices = graph.GetNextVertices(current);
        for (int i = 0; i < next_vertices.size(); ++i) {
            to = next_vertices[i].first;
            weight = next_vertices[i].second;
            if (dist[current] + weight < dist[to]) {
                bin_tree.erase(std::make_pair(dist[to], to));
                dist[to] = dist[current] + weight;
                bin_tree.insert(std::make_pair(dist[to], to));
            }
        }
    }
}

static const int MAX = 10e8;

int find_cost_of_teleportation(int a, int b, int m, int x, int y) {
    if (x == y) {
        return 0;
    }
    int ans = 0;
    Graph<std::pair<int, int>> graph = BuildGraph(a, b, m);

    std::vector<int> dist(m + 2);
    for (int i = 0; i <= m; ++i) {
        dist[i] = MAX;
    }

    Dijkstra(graph, x, dist);
    return dist[y];
}

int main() {
    int a;
    int b;
    int m;
    int x;
    int y;
    std::cin >> a >> b >> m >> x >> y;

    std::cout << find_cost_of_teleportation(a, b, m, x, y) << std::endl;
    return 0;
}