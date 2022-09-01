#include <iostream>
#include <queue>
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

    void AddVertex();
    void AddEdge(int from, int to);
    bool HasEdge(int from, int to) const;
    const std::vector<int> &GetNextVertices(int from) const;
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
        graph.AddEdge(to[i], from[i]);
    }
    return std::move(graph);
}

struct pair {
    int vertex;
    int dist;

    explicit pair(int v = -1, int d = -1) : vertex(v), dist(d) {}
};

void FindDistsBFS(int start, std::vector<int> &distances, Graph &graph) {
    std::queue<pair> queue;
    std::vector<bool> colour(graph.GetVertexCount(), false);
    pair top;
    colour[start] = true;
    queue.push(pair(start, 0));
    std::vector<int> neighbours;

    while (!queue.empty()) {
        top = queue.front();
        queue.pop();
        distances[top.vertex] = top.dist;
        neighbours = graph.GetNextVertices(top.vertex);
        for (int i = 0; i < neighbours.size(); ++i) {
            if (!colour[neighbours[i]]) {
                colour[neighbours[i]] = true;
                queue.push(pair(neighbours[i], top.dist + 1));
            }
        }
    }
}

int find_min_edges_count(int vertex_count, int edge_count, std::vector<int> &from, std::vector<int> &to, int leon,
                         int matilda, int milk) {
    int ans = 0;
    Graph graph = BuildGraph(vertex_count, edge_count, from, to);

    std::vector<int> leon_dist(vertex_count);
    std::vector<int> matilda_dist(vertex_count);
    std::vector<int> milk_dist(vertex_count);

    FindDistsBFS(leon, leon_dist, graph);
    FindDistsBFS(matilda, matilda_dist, graph);
    FindDistsBFS(milk, milk_dist, graph);


    for (int i = 1; i < vertex_count; ++i) {
        if (i == 1 || leon_dist[i] + matilda_dist[i] + milk_dist[i] <= ans) {
            ans = leon_dist[i] + matilda_dist[i] + milk_dist[i];
        }
    }
    return ans;
}

int main() {
    int vertex_count;
    int edge_count;
    int leon_pos;
    int matilda_pos;
    int milk_pos;
    std::cin >> vertex_count >> edge_count >> leon_pos >> matilda_pos >> milk_pos;
    ++vertex_count;

    std::vector<int> from(edge_count);
    std::vector<int> to(edge_count);

    for (int i = 0; i < edge_count; ++i) {
        std::cin >> from[i] >> to[i];
    }

    std::cout << find_min_edges_count(vertex_count, edge_count, from, to, leon_pos, matilda_pos, milk_pos) << '\n';
    return 0;
}