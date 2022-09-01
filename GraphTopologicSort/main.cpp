#include <algorithm>
#include <iostream>
#include <stack>
#include <vector>

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
    void AddEdge(int from, int to);
    [[nodiscard]] bool HasEdge(int from, int to) const;
    [[nodiscard]] const std::vector<int> &GetNextVertices(int from) const;
    [[nodiscard]] int GetVertexCount() const;

  private:
    int vertex_count;
    std::vector<std::vector<int> > graph;

};

void Graph::AddVertex() {
    ++vertex_count;
}

void Graph::AddEdge(int from, int to) {
    graph[from].push_back(to);
}

bool Graph::HasEdge(int from, int to) const {
    for (int i = 0; i < graph[from].size(); ++i) {
        if (graph[from][i] == to) {
            return true;
        }
    }
    return false;
}

const std::vector<int> &Graph::GetNextVertices(int from) const {
    return graph[from];
}

int Graph::GetVertexCount() const {
    return vertex_count;
}

Graph BuildGraph(int vertex_count, int edge_count, std::vector<int> &from, std::vector<int> &to) {
    Graph graph(vertex_count);
    for (int i = 0; i < edge_count; ++i) {
        graph.AddEdge(from[i], to[i]);
    }
    return graph;
}

bool TopologicalSort(Graph &graph, std::vector<int> &calls_order, std::vector<unsigned char> &colour, int &zeros_count,
                     int start) {
    std::stack<int> stack;
    std::vector<int> next_vertices;
    int vertex, index;

    stack.push(start);
    while (!stack.empty()) {
        vertex = stack.top();
        if (colour[vertex] == 0) {
            colour[vertex] = 1;
            next_vertices = graph.GetNextVertices(vertex);
            for (int i = next_vertices.size() - 1; i >= 0; --i) {
                if (colour[next_vertices[i]] == 0) {
                    stack.push(next_vertices[i]);
                } else if (colour[next_vertices[i]] == 1) {
                    return false;
                }
            }
        } else if (colour[vertex] == 1) {
            colour[vertex] = 2;
            calls_order[index++] = vertex;
            stack.pop();
        } else {
            stack.pop();
        }
    }

    return true;
}

bool IsOperationExecutable(Graph &graph, std::vector<int> &calls_order) {
    std::vector<unsigned char> colour(graph.GetVertexCount(), 0);
    int zeros_count = 0;
    bool result;
    for (int i = 0; i < graph.GetVertexCount(); ++i) {
        if (colour[i] == 0) {
            ++zeros_count;
            result = TopologicalSort(graph, calls_order, colour, zeros_count, i);
            if (!result) {
                return false;
            }
        }
    }

    if (zeros_count == 0) {
        return false;
    }
    return true;
}

void FindCallsOrder(int vertex_count, int edge_count, std::vector<int> &from, std::vector<int> &to) {
    Graph graph = BuildGraph(vertex_count, edge_count, from, to);
    std::vector<int> calls_order(vertex_count, 0);
    std::vector<int> right_order;

    if (IsOperationExecutable(graph, calls_order)) {
        std::cout << "YES" << std::endl;

        for (int i = 0; i < vertex_count; ++i) {
            std::cout << calls_order[vertex_count - 1 - i] << " ";
            right_order.push_back(calls_order[vertex_count - 1 - i]);
        }
        std::cout << '\n';
    } else {
        std::cout << "NO\n";
    }

}

int main() {
    int policeman_count;
    int channels_count;
    std::cin >> policeman_count >> channels_count;

    std::vector<int> from(channels_count);
    std::vector<int> to(channels_count);

    for (int i = 0; i < channels_count; ++i) {
        std::cin >> from[i] >> to[i];
    }

    FindCallsOrder(policeman_count, channels_count, from, to);
    return 0;
}