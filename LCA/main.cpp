#include <iostream>
#include <stack>
#include <vector>

struct Node {
    int vertex;
    int parent;

    Node(int v, int p) : vertex(v), parent(p) {}
};

void ReprocessDFS(std::vector<std::vector<int>> &graph_tree, std::vector<std::vector<int>> &jump,
                  std::vector<int> &in_time, std::vector<int> &out_time, int vertex, int parent) {
    std::stack<Node> stack;
    std::vector<bool> colour(graph_tree.size(), false);
    stack.push(Node(vertex, parent));
    int timer = 0;
    while (!stack.empty()) {
        vertex = stack.top().vertex;
        parent = stack.top().parent;
        if (!colour[vertex]) {
            colour[vertex] = true;
            in_time[vertex] = ++timer;
            jump[vertex][0] = parent;
            for (int i = 1; i < jump[vertex].size(); ++i) {
                jump[vertex][i] = jump[jump[vertex][i - 1]][i - 1];
            }
            for (int to: graph_tree[vertex]) {
                if (to != parent) {
                    stack.push(Node(to, vertex));
                }
            }
        } else {
            out_time[vertex] = ++timer;
            stack.pop();
        }
    }
}

int CalcLogarithm(int n) {
    int log = 0, value = 1;
    while (value <= n) {
        value <<= 1;
        ++log;
    }
    return log;

}

void SetJumpSize(int log, std::vector<std::vector<int>> &jump) {
    for (size_t i = 0; i < jump.size(); ++i) {
        jump[i].resize(log + 1, 0);
    }
}

int LCA(int u, int v, std::vector<int> &in_time, std::vector<int> &out_time,
        std::vector<std::vector<int>> &jump) {
    if (in_time[u] <= in_time[v] && out_time[u] >= out_time[v]) {
        return u;
    }
    if (in_time[v] <= in_time[u] && out_time[v] >= out_time[u]) {
        return v;
    }
    for (int i = jump[u].size() - 1; i >= 0; --i) {
        if (!(in_time[jump[u][i]] <= in_time[v] && out_time[jump[u][i]] >= out_time[v])) {
            u = jump[u][i];
        }
    }
    return jump[u][0];
}

long long ProcessQueries(int n, int m, std::vector<int> &in_time, std::vector<int> &out_time,
                         std::vector<std::vector<int> > &jump) {
    long long a1, a2, x, y, z;
    std::cin >> a1 >> a2 >> x >> y >> z;
    long long sum, v = sum = LCA(a1, a2, in_time, out_time, jump);
    while (m >= 2) {
        a1 = (x * a1 + y * a2 + z) % n;
        a2 = (x * a2 + y * a1 + z) % n;
        v = LCA((a1 + v) % n, a2, in_time, out_time, jump);
        sum += v;
        --m;
    }

    return sum;
}

long long ProcessLCAQueries(int n, int m, std::vector<std::vector<int> > &graph_tree) {
    std::vector<std::vector<int> > jump(n);
    int log_n = CalcLogarithm(n);
    SetJumpSize(log_n, jump);
    std::vector<int> in_time(n), out_time(n);
    ReprocessDFS(graph_tree, jump, in_time, out_time, 0, 0);
    return ProcessQueries(n, m, in_time, out_time, jump);
}

int main() {
    int n = 0;
    int m = 0;
    std::cin >> n >> m;
    std::vector<std::vector<int>> graph_tree(n);
    int parent;
    for (size_t i = 1; i < n; ++i) {
        std::cin >> parent;
        graph_tree[parent].emplace_back(i);
    }
    std::cout << ProcessLCAQueries(n, m, graph_tree) << '\n';
    return 0;
}