#include <iostream>
#include <stack>
#include <vector>

class Graph {
public:
    explicit Graph(int vert_count = 0);
    ~Graph();
    Graph(const Graph &) = delete;
    Graph(Graph &&) = delete;
    Graph &operator=(const Graph &) = delete;
    Graph &operator=(Graph &&) = delete;

    void AddVertex();
    void AddEdge(int from, int to);
    bool HasEdge(int from, int to) const;
    const std::vector<int> &GetNextVertices(int from) const;
    int GetVertexCount() const;
    void PrintAdjacencyMatrix() const;
    void TransitiveCloserDFS();
    void SetAdjacencyMatrixValue(int i, int j, int value);
    int GetAdjacencyMatrixValue(int i, int j);

  private:
    void DFS(int s, int v, bool first_time);

  private:
    std::vector<std::vector<int> > graph;
    std::vector<std::vector<int> > matrix;
    int vertex_count;

};

Graph::Graph(int vert_count) : vertex_count(vert_count) {
    std::vector<int> vec(vertex_count, 0);
    for (int i = 0; i < vert_count; ++i) {
        graph.emplace_back();
        matrix.push_back(vec);
    }
}

Graph::~Graph() {
    vertex_count = -1;
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

void Graph::PrintAdjacencyMatrix() const {
    for (int i = 0; i < vertex_count; ++i) {
        for (int j = 0; j < vertex_count; ++j)
            std::cout << matrix[i][j];
        std::cout << std::endl;
    }
}

void Graph::DFS(int s, int v, bool first_time) {
    if (!first_time) {
        matrix[s][v] = 1;
    }
    for (int i(0); i < graph[v].size(); ++i) {
        if (!matrix[s][graph[v][i]]) {
            DFS(s, graph[v][i], false);
        }
    }
}

void Graph::TransitiveCloserDFS() {
    for (int i(0); i < vertex_count; ++i) {
        DFS(i, i, true);
    }
}

void Graph::SetAdjacencyMatrixValue(int i, int j, int value) {
    matrix[i][j] = value;
}

int Graph::GetAdjacencyMatrixValue(int i, int j) {
    return matrix[i][j];
}

//// Solution

void SortVerticesByOut(Graph &graph, std::vector<int> &out) {
    std::vector<unsigned char> colour(graph.GetVertexCount(), 0);
    int vertex = 0, time = 0;
    std::stack<int> stack;
    std::vector<int> next_vertices;
    for (int i = 0; i < graph.GetVertexCount(); ++i) {
        if (colour[i] == 0) {
            stack.push(i);
            while (!stack.empty()) {
                vertex = stack.top();
                if (colour[vertex] == 0) {
                    colour[vertex] = 1;
                    next_vertices = graph.GetNextVertices(vertex);
                    for (int i = 0; i < next_vertices.size(); ++i) {
                        if (colour[next_vertices[i]] == 0) {
                            stack.push(next_vertices[i]);
                        }
                    }
                } else if (colour[vertex] == 1) {
                    colour[vertex] = 2;
                    out[time++] = vertex;
                    stack.pop();
                } else {
                    stack.pop();
                }
            }
        }
    }
}

void SetInitialValues(std::vector<int> &connected_components, std::vector<unsigned char> &colour, int n) {
    for (int i = 0; i < n; ++i) {
        connected_components[i] = 0;
        colour[i] = 0;
    }
    connected_components[n] = 0;
}

void BreakIntoConnectedComponents(Graph &graph, std::vector<int> &connected_components, Graph *inverted_graph) {
    int n = graph.GetVertexCount();
    std::vector<int> sorted_by_out_time_vertices(n);
    SortVerticesByOut(graph, sorted_by_out_time_vertices);
    std::vector<unsigned char> colour(n);
    SetInitialValues(connected_components, colour, n);

    int vertex = 0, number_of_component = -1;
    std::stack<int> stack;
    std::vector<int> next_vertices;
    for (int i = n - 1; i >= 0; --i) {
        if (colour[sorted_by_out_time_vertices[i]] == 0) {
            ++number_of_component;
            stack.push(sorted_by_out_time_vertices[i]);
            while (!stack.empty()) {
                vertex = stack.top();
                if (colour[vertex] == 0) {
                    colour[vertex] = 1;
                    next_vertices = inverted_graph->GetNextVertices(vertex);
                    for (int i = 0; i < next_vertices.size(); ++i) {
                        if (colour[next_vertices[i]] == 0) {
                            stack.push(next_vertices[i]);
                        }
                    }
                } else if (colour[vertex] == 1) {
                    colour[vertex] = 2;
                    connected_components[vertex] = number_of_component;
                    stack.pop();
                } else {
                    stack.pop();
                }
            }
        }
    }
    connected_components[graph.GetVertexCount()] = number_of_component + 1;
}

void BuildGraphConnectedComponentsGraph(const Graph &graph, Graph &connected_graph,
                                        std::vector<int> &connection_components) {
    int n = graph.GetVertexCount(), to;
    std::vector<int> next_vertices;
    for (int i = 0; i < n; ++i) {
        next_vertices = graph.GetNextVertices(i);
        for (int j = 0; j < next_vertices.size(); ++j) {
            to = next_vertices[j];
            if (connection_components[i] != connection_components[to]) {
                if (!connected_graph.HasEdge(connection_components[i], connection_components[to])) {
                    connected_graph.AddEdge(connection_components[i], connection_components[to]);
                } else if (i != to) {
                    if (!connected_graph.HasEdge(connection_components[i], connection_components[to])) {
                        connected_graph.AddEdge(connection_components[i], connection_components[to]);
                    }
                }
            }
        }
    }
}

void RecoverTransitionOfInitialGraph(Graph &graph, Graph &connected_graph, std::vector<int> components) {
    int n = graph.GetVertexCount();
    int h_comp;
    int j_comp;
    for (int h = 0; h < n; ++h) {
        for (int j = 0; j < n; ++j) {
            h_comp = components[h];
            j_comp = components[j];
            if (h_comp == j_comp) {
                if (connected_graph.GetAdjacencyMatrixValue(h_comp, j_comp) == 1) {
                    graph.SetAdjacencyMatrixValue(h, j, 1);
                } else if (h != j) {
                    graph.SetAdjacencyMatrixValue(h, j, 1);
                }
            } else if (connected_graph.GetAdjacencyMatrixValue(h_comp, j_comp) == 1) {
                graph.SetAdjacencyMatrixValue(h, j, 1);
            }
        }
    }
}

void MakeGraphTransitivelyClosed(Graph &graph, Graph &inverted_graph) {
    int n = graph.GetVertexCount();
    std::vector<int> connection_components(n + 1);

    BreakIntoConnectedComponents(graph, connection_components, inverted_graph);
    Graph connected_components_graph(connection_components[n]);

    BuildGraphConnectedComponentsGraph(graph, connected_components_graph, connection_components);
    connected_components_graph.TransitiveCloserDFS();

    RecoverTransitionOfInitialGraph(graph, connected_components_graph, connection_components);
    graph.PrintAdjacencyMatrix();
}

int main() {
    int n = 0;
    std::cin >> n;
    std::string str;
    Graph graph(n);
    Graph inverted_graph(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> str;
        for (int j = 0; j < n; ++j) {
            if (str[j] == '1') {
                graph.AddEdge(i, j);
                graph.SetAdjacencyMatrixValue(i, j, 1);
                inverted_graph.AddEdge(j, i);
            }
        }
    }

    MakeGraphTransitivelyClosed(graph, inverted_graph);
    return 0;
}