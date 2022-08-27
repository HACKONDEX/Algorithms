#include <iostream>
#include <string>
#include <vector>

template<class T>
class Graph {
  public:
    explicit Graph(int vert_count = 0) : vertex_count(vert_count) {
        for (int i = 0; i < vert_count; ++i)
            graph.emplace_back();
    }

    ~Graph() = default;
    Graph(const Graph &) = delete;
    Graph(Graph &&) = delete;

    Graph &operator=(const Graph &) = delete;
    Graph &operator=(Graph &&) = delete;

    void AddVertex();
    void AddEdge(int from, const T &to);
    bool HasEdge(int from, const T &to) const;
    const std::vector<T> &GetNextVertices(int from) const;
    int GetVertexCount();
    bool HasEdge(int vertex) const;
    void Print() const;

  private:
    int vertex_count;
    std::vector<std::vector<T> > graph;

};

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
int Graph<T>::GetVertexCount() {
    return vertex_count;
}

template<class T>
void Graph<T>::Print() const {
    for (int i = 0; i < vertex_count; ++i) {
        std::cout << "( " << i << " ) -->";
        for (int j = 0; j < graph[i].size(); ++j) {
            std::cout << " " << graph[i][j];
        }
        std::cout << std::endl;
    }
}

template<class T>
bool Graph<T>::HasEdge(int vertex) const {
    return !graph[vertex].empty();
}

//// Solution

int GetMatrixBridgeCount(std::vector<std::vector<int>>& int_bridge_matrix,
                         int n, int m, std::vector<std::string>& bridge) {
    int number = 1;
    int share = 1;
    for (int i = 0; i < n; ++i) {
        share = (i % 2 == 0) ? 1 : -1;
        for (int j = 0; j < m; ++j) {
            int_bridge_matrix[i][j] = (bridge[i][j] == '*') ? share * (number++) : 0;
            share *= -1;
        }
    }
    return number;
}

void BuildGraph(Graph<int>& graph, int n, int m, std::vector<std::vector<int>>& bridge) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (bridge[i][j] >= 0) {
                continue;
            }
            if (i - 1 >= 0 && bridge[i - 1][j] != 0) {
                graph.AddEdge(-bridge[i][j], bridge[i - 1][j]);
            }
            if (i + 1 < n && bridge[i + 1][j] != 0) {
                graph.AddEdge(-bridge[i][j], bridge[i + 1][j]);
            }
            if (j - 1 >= 0 && bridge[i][j - 1] != 0) {
                graph.AddEdge(-bridge[i][j], bridge[i][j - 1]);
            }
            if (j + 1 < m && bridge[i][j + 1] != 0) {
                graph.AddEdge(-bridge[i][j], bridge[i][j + 1]);
            }
        }
    }
}

bool Kuhn(int vertex, Graph<int>& graph, std::vector<int> &pairs,
          std::vector<bool>& colour) {
    if (colour[vertex]) {
        return false;
    }
    colour[vertex] = true;
    std::vector<int> next_vertices = graph.GetNextVertices(vertex);
    int to = 0;
    for (int i = 0; i < next_vertices.size(); ++i) {
        to = next_vertices[i];
        if (pairs[to] == -1 || Kuhn(pairs[to], graph, pairs, colour)) {
            pairs[to] = vertex;
            return true;
        }
    }
    return false;
}

void GetMaximalPairs(Graph<int>& graph, std::vector<int> &pairs) {
    int size = pairs.size();
    std::vector<bool> colour(size);
    for (int j = 1; j < size; ++j) {
        if (!graph.HasEdge(j)) {
            continue;
        }
        for (int i = 0; i < size; ++i) {
            colour[i] = false;
        }
        Kuhn(j, graph, pairs, colour);
    }
}

long long GetMinimalTime(std::vector<int> &pairs, int a, int b) {
    long long answer;
    int count = 0;
    int size = pairs.size();
    for (int i = 1; i < size; ++i) {
        if (pairs[i] != -1) {
            ++count;
        }
    }
    answer = count * a + ((size - 1) - 2 * count) * b;
    return answer;
}

long long GetMinimalTime(int n, int m, int a, int b, std::vector<std::string>& bridge) {
    std::vector<std::vector<int>> int_bridge_matrix(n);
    for (int i = 0; i < n; ++i) {
        int_bridge_matrix.emplace_back(std::vector<int>(m));
    }

    long long answer = 0;
    int max_vertex_count = GetMatrixBridgeCount(int_bridge_matrix, n, m, bridge);
    Graph<int> graph(max_vertex_count);
    if (a >= b * 2) {
        answer = (max_vertex_count - 1) * b;
    } else {
        BuildGraph(graph, n, m, int_bridge_matrix);
        std::vector<int> maximal_pairs(max_vertex_count, -1);
        GetMaximalPairs(graph, maximal_pairs);
        answer = GetMinimalTime(maximal_pairs, a, b);
    }
    return answer;
}

int main() {
    int n, m, a, b;
    std::cin >> n >> m >> a >> b;
    std::vector<std::string> bridge(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> bridge[i];
    }

    std::cout << GetMinimalTime(n, m, a, b, bridge) << std::endl;
    return 0;
}