#include <iostream>
#include <vector>

static const int MIN = -10e8;

class SegmentTree {
private:
    struct Node {
        int max;
        int add;

        explicit Node(int m = MIN, int a = 0) : max(m), add(a) {}
    };

public:
    SegmentTree(const std::vector<int>& data, int n);

    ~SegmentTree() = default;
    SegmentTree(const SegmentTree &) = delete;
    SegmentTree(SegmentTree &&) = delete;

    SegmentTree &operator=(const SegmentTree &) = delete;
    SegmentTree &operator=(SegmentTree &&) = delete;

    void AddOnSegment(int left, int right, int value);
    int MaxOnSegment(int left, int right);

private:
    void Propagate(int index, int left, int right);
    int MaxHelper(int left, int right, int index, int left_border, int right_border);
    void UpdateAdd(int left, int right, int index, int value, int left_border, int right_border);

    static int GetMinDegreeOfTwo(int n);

private:
    int nodes_count;
    int vertices_count;
    std::vector<Node> nodes;
};

SegmentTree::SegmentTree(const std::vector<int>& data, int n) {
    vertices_count = GetMinDegreeOfTwo(n);
    nodes_count = 2 * vertices_count - 1;
    nodes.assign(nodes_count, Node());
    for (int i = n - 1; i >= 0; --i) {
        nodes[nodes_count - vertices_count + i].max = data[i];
    }
    for (int i = nodes_count - vertices_count - 1; i >= 0; --i) {
        nodes[i].max = std::max(nodes[2 * i + 1].max, nodes[2 * i + 2].max);
    }
}

int SegmentTree::GetMinDegreeOfTwo(int n) {
    int value = 1;
    while (value < n)
        value *= 2;
    return value;
}

void SegmentTree::Propagate(int index, int left_border, int right_border) {
    if (right_border - left_border == 0) {
        return;
    }
    if (nodes[index].add == 0) {
        return;
    }
    nodes[2 * index + 1].add += nodes[index].add;
    nodes[2 * index + 2].add += nodes[index].add;
    nodes[2 * index + 1].max += nodes[index].add;
    nodes[2 * index + 2].max += nodes[index].add;
    nodes[index].add = 0;
}

int SegmentTree::MaxHelper(int left, int right, int index, int left_border, int right_border) {
    Propagate(index, left_border, right_border);
    if (left_border > right || left > right_border) {
        return MIN;
    }
    if (left_border >= left && right_border <= right) {
        return nodes[index].max;
    }
    int separation_index = (left_border + right_border) / 2;
    int left_sub_segment = MaxHelper(left, right, 2 * index + 1, left_border, separation_index);
    int right_sub_segment = MaxHelper(left, right, 2 * index + 2, separation_index + 1, right_border);
    return std::max(left_sub_segment, right_sub_segment);
}

void SegmentTree::UpdateAdd(int left, int right, int index, int value, int left_border, int right_border) {
    Propagate(index, left_border, right_border);
    if (left_border > right || left > right_border) {
        return;
    }
    if (left_border >= left && right_border <= right) {
        nodes[index].max += value;
        nodes[index].add += value;
        return;
    }
    int separation_index = (left_border + right_border) / 2;
    UpdateAdd(left, right, 2 * index + 1, value, left_border, separation_index);
    UpdateAdd(left, right, 2 * index + 2, value, separation_index + 1, right_border);
    nodes[index].max = std::max(nodes[2 * index + 1].max, nodes[2 * index + 2].max);
}

int SegmentTree::MaxOnSegment(int left, int right) {
    return MaxHelper(left, right, 0, 0, vertices_count - 1);
}

void SegmentTree::AddOnSegment(int left, int right, int value) {
    UpdateAdd(left, right, 0, value, 0, vertices_count - 1);
}

void
ProcessQueries(int n, int max_tickets_count, int requests_count, std::vector<int>& places, std::vector<int> &left,
               std::vector<int> &right, std::vector<int>& tickets) {
    SegmentTree segment_tree(places, n - 1);
    for (int i = 0; i < requests_count; ++i) {
        if (segment_tree.MaxOnSegment(left[i], right[i] - 1) + tickets[i] <= max_tickets_count) {
            segment_tree.AddOnSegment(left[i], right[i] - 1, tickets[i]);
        } else {
            std::cout << i << " ";
        }
    }
}

int main() {
    int n;
    int requests_count;
    int max_tickets_count;
    std::cin >> n;

    std::vector<int> places(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        std::cin >> places[i];
    }
    std::cin >> max_tickets_count >> requests_count;
    std::vector<int> left(requests_count);
    std::vector<int> right(requests_count);
    std::vector<int> tickets(requests_count);

    for (int i = 0; i < requests_count; ++i) {
        std::cin >> left[i] >> right[i] >> tickets[i];
    }

    ProcessQueries(n, max_tickets_count, requests_count, places, left, right, tickets);
    return 0;
}