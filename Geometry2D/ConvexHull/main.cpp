#include <cmath>
#include <iostream>
#include <vector>

struct Point {
    Point() : x(0), y(0) {}
    Point(int64_t x, int64_t y) : x(x), y(y) {}

    int64_t x;
    int64_t y;
};

struct GeomVector {

    GeomVector() : x(0), y(0) {}

    GeomVector(int64_t x, int64_t y) : x(x), y(y) {}

    GeomVector(const Point &a, const Point &b) : x(b.x - a.x), y(b.y - a.y) {}

    int64_t x;
    int64_t y;
};

int64_t VectorProduct(const GeomVector &a, const GeomVector &b) {
    return a.x * b.y - a.y * b.x;
}

long double DistOfPoints(const Point &a, const Point &b) {
    return std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
}

class ComparePointsByAngle {
public:

    ComparePointsByAngle(const Point &point_) : base_point(point_) {}

    bool operator()(const Point &first, const Point &second) {
        GeomVector first_vector(base_point, first);
        GeomVector second_vector(base_point, second);
        int64_t vector_product = VectorProduct(first_vector, second_vector);
        if (vector_product > 0) {
            return true;
        } else if (vector_product < 0) {
            return false;
        }
        return DistOfPoints(base_point, first) < DistOfPoints(base_point, second);
    }

private:
    Point base_point;
};

void GetMostLeftBottomPoint(Point &left_bottom_point, const std::vector<Point> &points_) {
    left_bottom_point = points_[0];
    for (const auto &current_point: points_) {
        if (current_point.y < left_bottom_point.y) {
            left_bottom_point = current_point;
        } else if (current_point.y == left_bottom_point.y && current_point.x < left_bottom_point.x) {
            left_bottom_point = current_point;
        }
    }
}

void GetConvexHull(std::vector<Point> &convex_hull, std::vector<Point> &points_) {
    Point compare_base;
    GetMostLeftBottomPoint(compare_base, points_);
    std::sort(points_.begin(), points_.end(), ComparePointsByAngle(compare_base));
    std::unique(points_.begin(), points_.end(), ComparePointsByAngle(compare_base));
    Point top;
    Point bottom;
    bool was_point_added = false;
    for (const auto &current_point: points_) {
        was_point_added = false;
        if (convex_hull.size() < 2) {
            convex_hull.push_back(current_point);
            continue;
        }
        ComparePointsByAngle compare_base(current_point);
        while (convex_hull.size() >= 2) {
            top = convex_hull.back();
            bottom = *(convex_hull.end() - 2);
            GeomVector current_top(current_point, top);
            GeomVector current_bottom(current_point, bottom);
            if (!compare_base.operator()(top, bottom)) {
                convex_hull.push_back(current_point);
                was_point_added = true;
                break;
            }
            convex_hull.pop_back();
        }
        if (!was_point_added) {
            convex_hull.push_back(current_point);
        }
    }
}

long double GetPerimeterOfConvexHull(std::vector<Point> &convex_hull) {
    long double perimeter = 0;
    if (convex_hull.size() == 1) {
        return 0;
    }
    if (convex_hull.size() > 2) {
        convex_hull.push_back(convex_hull.front());
    }
    for (size_t i = 0; i < convex_hull.size() - 1; ++i) {
        perimeter += DistOfPoints(convex_hull[i], convex_hull[i + 1]);
    }
    return perimeter;
}

long double GetMinimalFenceLength(const std::vector<Point> &points) {
    std::vector<Point> copied_points = points;
    std::vector<Point> convex_hull;
    GetConvexHull(convex_hull, copied_points);
    return GetPerimeterOfConvexHull(convex_hull);
}

int main() {
    std::cout.precision(10);
    int64_t points_count = 0;
    std::cin >> points_count;
    std::vector<Point> points(points_count);
    for (size_t i = 0; i < points_count; ++i) {
        std::cin >> points[i].x >> points[i].y;
    }
    std::cout << GetMinimalFenceLength(points) << std::endl;
    return 0;
}