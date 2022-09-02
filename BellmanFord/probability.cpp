#include <iostream>
#include <iomanip>
#include <vector>

struct Edge {
    int from;
    int to;
    double cost;

    Edge(int f, int t, double c) : from(f), to(t), cost(c) {}
};

static const double MIN = -2;

void BuildEdgeGraph(std::vector<Edge> &graph, int m, const std::vector<int> &from, const std::vector<int> &to,
                    const std::vector<double> &cost) {
    double c;
    for (int i = 0; i < m; ++i) {
        c = 100 - cost[i];
        c /= 100;
        graph.emplace_back(Edge(from[i], to[i], c));
        graph.emplace_back(Edge(to[i], from[i], c));
    }
}

void BellmanFord(std::vector<Edge> &graph, int vert_count, std::vector<double>& dist, int start, int finish) {
    int to, from;
    dist[start] = 1;
    for (int i = 0; i < vert_count; ++i) {
        for (auto &j: graph) {
            from = j.from;
            to = j.to;

            if (from == finish || to == start) {
                continue;
            }

            if (dist[from] != MIN && dist[to] < dist[from] * j.cost) {
                dist[to] = dist[from] * j.cost;
            }
        }
    }
}

double FindMinimalProbability(int vert_count, int m, int start, int finish, const std::vector<int> &from,
                              const std::vector<int> &to,
                              const std::vector<double> &cost) {
    std::vector<Edge> graph;
    BuildEdgeGraph(graph, m, from, to, cost);
    std::vector<double> prob(vert_count + 1, MIN);
    BellmanFord(graph, vert_count, prob, start, finish);
    return 1 - prob[finish];
}

int main() {
    std::setprecision(6);
    int n, m, s, f;
    std::cin >> n >> m >> s >> f;
    std::vector<int> from(m);
    std::vector<int> to(m);
    std::vector<double> cost(m);

    for (int i = 0; i < m; ++i) {
        std::cin >> from[i] >> to[i] >> cost[i];
    }

    std::cout << FindMinimalProbability(n, m, s, f, from, to, cost) << std::endl;
    return 0;
}