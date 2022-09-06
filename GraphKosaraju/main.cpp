#include <iostream>
#include <stack>
#include <vector>

class Graph {
  public:
    explicit Graph(int vert_count = 0) : vertex_count(vert_count) {
        for (int i = 0; i < vert_count; ++i) {
            graph.emplace_back();
        }
    }
    ~Graph();
    Graph(const Graph &) = default;
    Graph(Graph &&) = default;

    Graph &operator=(const Graph &) = delete;
    Graph &operator=(Graph &&) = delete;

    void AddVertex();
    void AddEdge(int from, int to);
    [[nosicard]] bool HasEdge(int from, int to) const;
    [[nodiscard]] const std::vector<int> &GetNextVertices(int from) const;
    int GetVertexCount() const;

  private:
    int vertex_count;
    std::vector<std::vector<int> > graph;

};

Graph::~Graph() {
    vertex_count = -1;
}

void Graph::AddVertex() {
    ++vertex_count;
    graph.emplace_back();
}

void Graph::AddEdge(int from, int to) {
    graph[from].push_back(to);
}

bool Graph::HasEdge(int from, int to) const {
    for (int x: graph[from]) {
        if (x == to) {
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

//// Solution

void BuildGraph(Graph &graph, int edge_count, std::vector<int> &from, std::vector<int> &to) {
    for (int i = 0; i < edge_count; ++i) {
        graph.AddEdge(from[i] - 1, to[i] - 1);
    }
}

void SortVerticesByOut(Graph &graph, std::vector<int> &out) {
    std::vector<unsigned char> colour(graph.GetVertexCount(), 0);
    int vertex = 0;
    int time = 0;
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

void BreakIntoConnectedComponents(Graph &graph, std::vector<int> &connected_components, Graph &inverted_graph) {
    int n = graph.GetVertexCount();
    std::vector<int> sorted_by_out_time_vertices(n);
    SortVerticesByOut(graph, sorted_by_out_time_vertices);
    std::vector<unsigned char> colour(n);
    SetInitialValues(connected_components, colour, n);

    int vertex = 0;
    int number_of_component = 0;

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
                    next_vertices = inverted_graph.GetNextVertices(vertex);
                    for (auto x: next_vertices) {
                        if (colour[x] == 0) {
                            stack.push(x);
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
    connected_components[graph.GetVertexCount()] = number_of_component;
}

void
FindStartsAndEnds(Graph &graph, std::vector<int> &connected_components, std::vector<std::pair<int, int>> &starts_ends) {
    for (int i = 0; i <= connected_components[graph.GetVertexCount()]; ++i) {
        starts_ends[0].first = 0;
        starts_ends[0].second = 0;
    }
    std::vector<int> next_vertices;
    for (int i = 0; i < graph.GetVertexCount(); ++i) {
        next_vertices = graph.GetNextVertices(i);
        for (int x: next_vertices) {
            if (connected_components[i] != connected_components[x]) {
                ++starts_ends[connected_components[i]].second;
                ++starts_ends[connected_components[x]].first;
            }
        }
    }
}

void FindInOutNodesCount(int &incoming, int &outcoming, std::vector<std::pair<int, int>> &starts_ends, int count) {
    for (int i = 1; i <= count; ++i) {
        if (starts_ends[i].first == 0) {
            ++incoming;
        }
        if (starts_ends[i].second == 0) {
            ++outcoming;
        }
    }
}

int FindNumberOfStreetsToAdd(int offices_count, int streets_count, std::vector<int> &from, std::vector<int> &to) {
    int ans = 0;

    Graph graph(offices_count);
    BuildGraph(graph, streets_count, from, to);
    Graph inverted_graph(offices_count);
    BuildGraph(inverted_graph, streets_count, to, from);

    std::vector<int> connected_components(offices_count + 1);
    BreakIntoConnectedComponents(graph, connected_components, inverted_graph);
    std::vector<std::pair<int, int>> starts_ends(connected_components[graph.GetVertexCount()] + 1);
    FindStartsAndEnds(graph, connected_components, starts_ends);

    if (connected_components[offices_count] == 1) {
        ans = 0;
    } else {
        int incoming = 0;
        int outcoming = 0;
        FindInOutNodesCount(incoming, outcoming, starts_ends, connected_components[offices_count]);
        ans = (incoming >= outcoming ? incoming : outcoming);
    }

    return ans;
}

int main() {
    int offices_count = 0;
    int streets_count = 0;
    std::cin >> offices_count >> streets_count;

    std::vector<int> from(streets_count);
    std::vector<int> to(streets_count);

    for (int i = 0; i < streets_count; ++i) {
        std::cin >> from[i] >> to[i];
    }

    std::cout << FindNumberOfStreetsToAdd(offices_count, streets_count, from, to) << std::endl;
    return 0;
}
