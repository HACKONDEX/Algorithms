#include <geometry.hpp>

int main() {
    Rectangle x(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0));
    x.Rotate(Point(0.5, 0.5), 90);
    return 0;
}