#include <iomanip>
#include <iostream>
#include <vector>

static const double MINIM = -10e6;

struct Edge {
    int from_vertex;
    int to_vertex;
    double tarif;
    double commis;

    explicit Edge(int from = 0, int ver = 0, double t = 0, double c = 0) : from_vertex(from), to_vertex(ver), tarif(t),
                                                                           commis(c) {}
};

template<typename T>
void Copy(const std::vector<T> &this_one, std::vector<T> &into_this, int size) {
    for (int i = 0; i < size; ++i)
        into_this[i] = this_one[i];
}

void BuildEdgeGraph(std::vector<Edge> &edges, int m, std::vector<int> &from, std::vector<int> &to,
                    std::vector<double> &a_b, std::vector<double> &b_a, std::vector<double> &c_a_b,
                    std::vector<double> &c_b_a) {
    int i = 0;
    for (int j = 0; j < m; ++j) {
        edges[i++] = Edge(from[j], to[j], a_b[j], c_a_b[j]);
        edges[i++] = Edge(to[j], from[j], b_a[j], c_b_a[j]);
    }
}

void BellmanFordNegativeCycle(std::vector<Edge> &edges, int n, int edges_count, std::vector<double> &d_last) {
    for (int j(0); j < n; ++j) {
        for (int i(0); i < edges_count; ++i) {
            int from = edges[i].from_vertex;
            if (d_last[from] == MINIM) {
                continue;
            }
            int to = edges[i].to_vertex;
            double tarif = edges[i].tarif;
            double commis = edges[i].commis;
            if (d_last[to] == MINIM || tarif * (d_last[from] - commis) > d_last[to]) {
                d_last[to] = tarif * (d_last[from] - commis);
            }
        }
    }
}

void GetAnswer(const std::vector<double> &d_next, std::vector<double> &d_last, int n) {
    for (int i(1); i <= n; ++i) {
        if (d_next[i] > d_last[i]) {
            std::cout << "YES" << std::endl;
            return;
        }
    }
    std::cout << "NO" << std::endl;
}


void CanIncrease(int n, int m, int s, double value, std::vector<int> &from, std::vector<int> &to,
                 std::vector<double> &from_tarif_to,
                 std::vector<double> &from_commission_to, std::vector<double> &to_tarif_from,
                 std::vector<double> &to_commission_from) {
    std::vector<Edge> edges(2 * m);
    BuildEdgeGraph(edges, m, from, to, from_tarif_to, to_tarif_from, from_commission_to, to_commission_from);
    std::vector<double> d_last(n + 1, MINIM);
    std::vector<double> d_next(n + 1);

    d_last[s] = value;
    BellmanFordNegativeCycle(edges, n, 2 * m, d_last);
    Copy(d_next, d_last, n + 1);
    BellmanFordNegativeCycle(edges, 1, 2 * m, d_next);
    GetAnswer(d_next, d_last, n);
}

int main() {
    int n, m, s;
    double value;
    std::setprecision(8);
    std::cin >> n >> m >> s >> value;
    std::vector<int> from(m);
    std::vector<int> to(m);
    std::vector<double> from_tarif_to(m);
    std::vector<double> from_commission_to(m);
    std::vector<double> to_tarif_from(m);
    std::vector<double> to_commission_from(m);
    for (int i(0); i < m; ++i) {
        std::cin >> from[i] >> to[i];
        std::cin >> from_tarif_to[i] >> from_commission_to[i];
        std::cin >> to_tarif_from[i] >> to_commission_from[i];

    }
    CanIncrease(n, m, s, value, from, to, from_tarif_to, from_commission_to, to_tarif_from,
                to_commission_from);
    return 0;
}