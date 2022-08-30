#include <iostream>
#include <vector>

struct Edge {
    int to;
    int capacity;
    int flow;
    int friend_edge_number;

    explicit Edge(int t = 0, int c = 0, int f = 0, int i = 0) : to(t), capacity(c), flow(f), friend_edge_number(i) {}
};

int start = 0;
int target = 0;

void BuildGraph(int n, std::vector<int> &from, std::vector<int> &to, std::vector<std::vector<Edge> > &graph,
                std::vector<std::vector<Edge> > &indirect_graph) {
    std::vector<Edge> a;
    int m = to.size();
    for (int i = 0; i <= n; ++i) {
        graph.push_back(a);
        indirect_graph.push_back(a);
    }
    for (int i = 0; i < m; ++i) {
        if (from[i] == to[i]) {
            continue;
        }
        int size = graph[from[i]].size();
        int indirect_size = indirect_graph[to[i]].size();
        graph[from[i]].push_back(Edge(to[i], 1, 0, indirect_size));
        indirect_graph[to[i]].push_back(Edge(from[i], 0, 0, size));
    }
}

int DFS(int vertex, int min_flow, std::vector<bool> &marked, std::vector<std::vector<Edge> > &graph,
        std::vector<std::vector<Edge> > &indirect_graph) {
    if (marked[vertex]) {
        return 0;
    }
    marked[vertex] = true;
    if (vertex == target) {
        return min_flow;
    }
    for (int i = 0; i < graph[vertex].size(); ++i) {
        Edge u = graph[vertex][i];
        if (u.flow < u.capacity) {
            int new_delta = DFS(u.to, std::min(min_flow, u.capacity - u.flow), marked, graph, indirect_graph);
            if (new_delta > 0) {
                graph[vertex][i].flow += new_delta;
                indirect_graph[u.to][u.friend_edge_number].flow -= new_delta;
                return new_delta;
            }
        }
    }
    for (int i = 0; i < indirect_graph[vertex].size(); ++i) {
        Edge u = indirect_graph[vertex][i];
        if (u.flow < u.capacity) {
            int new_delta = DFS(u.to, std::min(min_flow, u.capacity - u.flow), marked, graph, indirect_graph);
            if (new_delta > 0) {
                indirect_graph[vertex][i].flow += new_delta;
                graph[u.to][u.friend_edge_number].flow -= new_delta;
                return new_delta;
            }
        }
    }
    return 0;
}

int
DFSForMaxFlow(int vertex, std::vector<int> &path, std::vector<bool> &marked, std::vector<std::vector<Edge> > &graph) {
    if (marked[vertex]) {
        return 0;
    }
    marked[vertex] = true;
    if (vertex == target) {
        return target;
    }
    for (int i = 0; i < graph[vertex].size(); ++i) {
        Edge u = graph[vertex][i];
        if (u.capacity == u.flow) {
            int new_vertex = DFSForMaxFlow(u.to, path, marked, graph);
            if (new_vertex > 0) {
                path.push_back(new_vertex);
                graph[vertex][i].capacity += 1;
                return vertex;
            }
        }
    }
    return 0;
}

int FordFulkerson(std::vector<bool> &marked, std::vector<std::vector<Edge> > &graph,
                  std::vector<std::vector<Edge> > &indirect_graph) {
    int answer = 0;
    while (true) {
        marked.assign(marked.size(), false);
        int flow = DFS(start, 1000000, marked, graph, indirect_graph);
        if (flow > 0) {
            answer += flow;
        } else {
            break;
        }
    }
    return answer;
}

void PrintPath(std::vector<bool> &marked, std::vector<std::vector<Edge> > &graph) {
    std::vector<int> path;
    marked.assign(marked.size(), false);
    DFSForMaxFlow(start, path, marked, graph);
    path.push_back(start);
    for (int i = path.size() - 1; i >= 0; --i) {
        std::cout << path[i] << " ";
    }
    std::cout << std::endl;
}

void solve(int n, std::vector<int> &to, std::vector<int> &from) {
    std::vector<std::vector<Edge> > graph;
    std::vector<std::vector<Edge> > indirect_graph;
    BuildGraph(n, from, to, graph, indirect_graph);
    std::vector<bool> marked(n + 1, false);
    if (FordFulkerson(marked, graph, indirect_graph) < 2) {
        std::cout << "NO\n";
        return;
    }
    std::cout << "YES\n";
    PrintPath(marked, graph);
    PrintPath(marked, graph);
}


int main() {
    int n, m;
    std::cin >> n >> m >> start >> target;
    std::vector<int> from(m);
    std::vector<int> to(m);

    for (int i = 0; i < m; ++i) {
        std::cin >> from[i] >> to[i];
    }

    solve(n, to, from);
    return 0;
}