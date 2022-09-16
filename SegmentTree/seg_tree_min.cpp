#include <iostream>
#include <vector>

static const int INF = 10e8;

class SegmentTree {
private:
    struct Node {
        int min;
        int assign;

        explicit Node(int m = INF, int a = -1) : min(m), assign(a) {}
    };

public:
    SegmentTree(const std::vector<int>& data, int m);
    ~SegmentTree() = default;
    SegmentTree(const SegmentTree &) = delete;
    SegmentTree(SegmentTree &&) = delete;

    SegmentTree &operator=(const SegmentTree &) = delete;
    SegmentTree &operator=(SegmentTree &&) = delete;

    void SetOnSegment(int left, int right, int value);
    int MinOnSegment(int left, int right);

private:
    void Propagate(int index, int left, int right);
    int MinHelper(int left, int right, int index, int left_border, int right_border);
    void Update(int left, int right, int index, int value, int left_border, int right_border);

    static int MinDegreeOfTwo(int n);


private:
    int nodes_count;
    int vertices_count;
    int n;
    std::vector<Node> nodes;
};

SegmentTree::SegmentTree(const std::vector<int>& data, int m) {
    n = m;
    vertices_count = MinDegreeOfTwo(n);
    nodes_count = 2 * vertices_count - 1;
    nodes.assign(nodes_count, Node());
    for (int i = n - 1; i >= 0; --i) {
        nodes[nodes_count - vertices_count + i].min = data[i];
    }
    for (int i = nodes_count - vertices_count - 1; i >= 0; --i) {
        nodes[i].min = std::min(nodes[2 * i + 1].min, nodes[2 * i + 2].min);
    }
}

int SegmentTree::MinDegreeOfTwo(int n) {
    int value = 1;
    while (value < n) {
        value *= 2;
    }
    return value;
}

void SegmentTree::Propagate(int index, int left_border, int right_border) {
    if (right_border - left_border == 0) {
        return;
    }
    if (nodes[index].assign == -1) {
        return;
    }
    nodes[2 * index + 1].assign = nodes[2 * index + 2].assign = nodes[2 * index + 1].min = nodes[2 * index +
                                                                                                 2].min = nodes[index].assign;
    nodes[index].assign = -1;
}

int SegmentTree::MinHelper(int left, int right, int index, int left_border, int right_border) {
    Propagate(index, left_border, right_border);
    if (left_border > right || left > right_border) {
        return INF;
    }
    if (left_border >= left && right_border <= right) {
        return nodes[index].min;
    }
    int separation_index = (left_border + right_border) / 2;
    int left_sub_segment = MinHelper(left, right, 2 * index + 1, left_border, separation_index);
    int right_sub_segment = MinHelper(left, right, 2 * index + 2, separation_index + 1, right_border);
    return std::min(left_sub_segment, right_sub_segment);
}

void SegmentTree::Update(int left, int right, int index, int value, int left_border, int right_border) {
    Propagate(index, left_border, right_border);
    if (left_border > right || left > right_border) {
        return;
    }
    if (left_border >= left && right_border <= right) {
        nodes[index].min = nodes[index].assign = value;
        return;
    }
    int separation_index = (left_border + right_border) / 2;
    Update(left, right, 2 * index + 1, value, left_border, separation_index);
    Update(left, right, 2 * index + 2, value, separation_index + 1, right_border);
    nodes[index].min = std::min(nodes[2 * index + 1].min, nodes[2 * index + 2].min);
}

int SegmentTree::MinOnSegment(int left, int right) {
    return MinHelper(left, right, 0, 0, vertices_count - 1);
}

void SegmentTree::SetOnSegment(int left, int right, int value) {
    Update(left, right, 0, value, 0, vertices_count - 1);
}

void ProcessQueries(int n, std::vector<int>& numbers) {
    SegmentTree segment_tree(numbers, n);
    int k;
    int r;
    int g;
    int b;
    int change_request;
    int change_left;
    int change_right;
    int min_left;
    int min_right;
    std::cin >> k;
    std::vector<int> answers;
    for (int i = 0; i < k; ++i) {
        std::cin >> change_left >> change_right >> r >> g >> b >> min_left >> min_right;
        change_request = r + g + b;
        segment_tree.SetOnSegment(change_left, change_right, change_request);
        answers.push_back(segment_tree.MinOnSegment(min_left, min_right));
    }
    for (const auto& x : answers) {
        std::cout << x << " ";
    }
}

int main() {
    int n;
    std::cin >> n;
    std::vector<int> numbers(n);
    int r;
    int g;
    int b;
    for (int i = 0; i < n; ++i) {
        std::cin >> r >> g >> b;
        numbers[i] = r + g + b;
    }
    ProcessQueries(n, numbers);
    return 0;
}