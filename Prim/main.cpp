#include <iostream>
#include <set>
#include <vector>

template<class T>
class Graph {
  public:
    explicit Graph(int vert_count = 0);
    ~Graph() = default;
    Graph(const Graph &) = delete;
    Graph(Graph &&) = delete;

    Graph &operator=(const Graph &) = delete;
    Graph &operator=(Graph &&) = delete;

    void AddVertex();
    void AddEdge(int from, const T &to);
    bool HasEdge(int from, const T &to) const;
    const std::vector<T> &GetNextVertices(int from) const;
    int GetVertexCount() const;

  private:
    int vertex_count;
    std::vector<std::vector<T>> graph;
};

template<typename T>
Graph<T>::Graph(int vert_count) : vertex_count(vert_count) {
    for (int i = 0; i < vert_count; ++i) {
        graph.emplace_back();
    }
}

template<class T>
void Graph<T>::AddVertex() {
    ++vertex_count;
}

template<class T>
void Graph<T>::AddEdge(int from, const T &to) {
    graph[from].push_back(to);
}

template<class T>
bool Graph<T>::HasEdge(int from, const T &to) const {
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

//// Solution
void BuildGraph(Graph<std::pair<int, int>> &graph, int m,
                std::vector<int> &from, std::vector<int> &to,
                std::vector<int> &cost) {
    for (int i = 0; i < m; ++i) {
        graph.AddEdge(from[i], std::make_pair(cost[i], to[i]));
        graph.AddEdge(to[i], std::make_pair(cost[i], from[i]));
    }
}

long long Prim(Graph<std::pair<int, int>> &graph) {
    long long sum = 0;
    int vertex_count = graph.GetVertexCount() - 1, current_vertex = 1;
    std::vector<bool> colour(vertex_count + 1);
    for (int i = 0; i <= vertex_count; ++i) {
        colour[i] = false;
    }

    std::set<std::pair<int, int>> min_heap;
    std::vector<std::pair<int, int>> next_vertices;
    std::pair<int, int> current_edge;

    while (vertex_count > 0) {
        colour[current_vertex] = true;
        next_vertices = graph.GetNextVertices(current_vertex);

        for(const auto& vert: next_vertices) {
            min_heap.insert(vert);
        }

        while (!min_heap.empty()) {
            current_edge = *min_heap.begin();
            min_heap.erase(min_heap.begin());

            if (!colour[current_edge.second]) {
                sum += current_edge.first;
                current_vertex = current_edge.second;
                break;
            }
        }
        --vertex_count;
    }

    return sum;
}

long long MSTWeight(int n, int m, std::vector<int> &from,
                    std::vector<int> &to, std::vector<int> &cost) {
    Graph<std::pair<int, int>> graph(n + 1);
    BuildGraph(graph, m, from, to, cost);
    return Prim(graph);
}

int main() {
    int n, m;
    std::cin >> n >> m;
    std::vector<int> from(m);
    std::vector<int> to(m);
    std::vector<int> cost(m);
    for (int i = 0; i < m; ++i) {
        std::cin >> from[i] >> to[i] >> cost[i];
    }

    std::cout << MSTWeight(n, m, from, to, cost) << std::endl;
    return 0;
}