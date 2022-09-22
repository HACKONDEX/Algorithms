#include <iostream>
#include <vector>

struct Point {
    Point() : x(0), y(0) {}
    Point(int64_t x, int64_t y) : x(x), y(y) {}

    int64_t x;
    int64_t y;
};

struct Segment {
    Segment() = default;
    Segment(const Point &a, const Point &b) : a(a), b(b) {}
    Segment(int64_t a_x, int64_t a_y, int64_t b_x, int64_t b_y) : a(Point(a_x, a_y)), b(Point(b_x, b_y)) {}

    Point a;
    Point b;
};

bool IfProjectionOnAxisIntersect(int64_t x_a, int64_t x_b, int64_t x_c, int64_t x_d) {
    if (x_a > x_b) {
        std::swap(x_a, x_b);
    }
    if (x_c > x_d) {
        std::swap(x_c, x_d);
    }
    return std::max(x_a, x_c) <= std::min(x_b, x_d);
}

int64_t OrientedArea(const Point &base, const Point &right, const Point &left) {
    return (right.x - base.x) * (left.y - base.y) - (right.y - base.y) * (left.x - base.x);
}

bool DoSegmentsIntersect(const Segment &first, const Segment &second) {
    if (!(IfProjectionOnAxisIntersect(first.a.x, first.b.x, second.a.x, second.b.x) &&
          IfProjectionOnAxisIntersect(first.a.y, first.b.y, second.a.y, second.b.y))) {
        return false;
    }
    bool on_different_sides_first =
            OrientedArea(first.a, first.b, second.a) * OrientedArea(first.a, first.b, second.b) <= 0;
    bool on_different_sides_second =
            OrientedArea(second.a, second.b, first.a) * OrientedArea(second.a, second.b, first.b) <= 0;
    return on_different_sides_first && on_different_sides_second;
}

int main() {
    Segment main_seg;
    std::cin >> main_seg.a.x >> main_seg.a.y >> main_seg.b.x >> main_seg.b.y;
    int64_t n;
    std::cin >> n;
    std::vector<Segment> segments;
    segments.reserve(n);
    Segment other;
    int64_t count = 0;
    for (int i = 0; i < n; ++i) {
        std::cin >> other.a.x >> other.a.y >> other.b.x >> other.b.y;
        if (DoSegmentsIntersect(main_seg, other)) {
            ++count;
        }
    }

    std::cout << count << '\n';
    return 0;
}
