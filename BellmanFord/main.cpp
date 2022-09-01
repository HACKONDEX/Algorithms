#include <iostream>
#include <vector>

const long long MAX = 10e15;

struct Edge {
    int from;
    int to;
    int cost;
};

std::vector<Edge>
BuildEdgeGraph(int edge_count, std::vector<int> &from, std::vector<int> &to, std::vector<int> &cost) {
    std::vector<Edge> graph(edge_count);
    for (int i = 0; i < edge_count; ++i) {
        graph[i].from = from[i];
        graph[i].to = to[i];
        graph[i].cost = cost[i];
    }

    return std::move(graph);
}

std::vector<long long> DistWithInfinity(int ver_count) {
    std::vector<long long> dist(ver_count + 1, MAX);
    return std::move(dist);
}

void
BellmanFordModified(std::vector<long long> &dist, std::vector<Edge> &graph, int vert_count, int edge_count, int start,
                    int max_verts) {
    std::vector<bool> colour(vert_count + 1);

    int to;
    int from;
    int cost;
    dist[start] = 0;
    for (int i = 0; i < max_verts; ++i) {
        for (int j = 0; j <= vert_count; ++j) {
            colour[j] = false;
        }

        for (int j = 0; j < edge_count; ++j) {
            from = graph[j].from;
            if (dist[from] < MAX && !colour[from]) {
                to = graph[j].to;
                cost = graph[j].cost;
                if (dist[to] > dist[from] + cost) {
                    dist[to] = dist[from] + cost;
                    colour[to] = true;
                }
            }
        }
    }
}

long long
FindMinimalCost(int vert_count, int edge_count, int max_vert, int start, int finish, std::vector<int> &from,
                std::vector<int> &to,
                std::vector<int> &cost) {
    if (start == finish) {
        return 0;
    }

    std::vector<Edge> graph = BuildEdgeGraph(edge_count, from, to, cost);
    std::vector<long long> dist = DistWithInfinity(vert_count);

    BellmanFordModified(dist, graph, vert_count, edge_count, start, max_vert);

    long long ans = dist[finish] < MAX ? dist[finish] : -1;
    return ans;
}

int main() {
    int n, m, k, s, f;
    std::cin >> n >> m >> k >> s >> f;

    std::vector<int> from(m);
    std::vector<int> to(m);
    std::vector<int> cost(m);
    for (int i = 0; i < m; ++i) {
        std::cin >> from[i] >> to[i] >> cost[i];
    }

    std::cout << FindMinimalCost(n, m, k, s, f, from, to, cost) << std::endl;
    return 0;
}