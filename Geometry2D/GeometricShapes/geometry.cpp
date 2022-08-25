#include <algorithm>
#include <cmath>
#include <geometry.hpp>

const double PI = std::acos(-1);
const double EPS = 0.0000000001;

/// Point Methods

bool Point::operator==(const Point &other) const {
    return EqualWithEps(x, other.x) && EqualWithEps(y, other.y);
}

bool Point::operator!=(const Point &other) const {
    return !(*this == other);
}

Point &Point::operator-=(const Point &other) {
    x -= other.x;
    y -= other.y;
    return *this;
}

Point &Point::operator+=(const Point &other) {
    x += other.x;
    y += other.y;
    return *this;
}

Point &Point::operator*=(double k) {
    x *= k;
    y *= k;
    return *this;
}

void Point::Reflect(Point regarding_other) {
    x = 2 * regarding_other.x - x;
    y = 2 * regarding_other.y - y;
}

void Point::Rotate(Point centre, double angle) {
    double radian_angle = (angle * PI) / double(180.0);
    double sin = std::sin(radian_angle);
    double cos = std::cos(radian_angle);

    //// Change cartesian centre
    Point shifted_point = *this;
    shifted_point -= centre;

    //// Rotate cartesian system -Angle
    x = shifted_point.x * cos - shifted_point.y * sin;
    y = shifted_point.x * sin + shifted_point.y * cos;

    //// Recover centre
    *this += centre;
}

void Point::Scale(Point centre, double coefficient) {
    Point direction_vector = *this;
    direction_vector -= centre;
    --coefficient;
    direction_vector *= coefficient;
    *this += direction_vector;
}

void Point::Reflect(Line axis) {
    Point centre = FindNormalsFundament(*this, axis);
    Reflect(centre);
}

/// Line Methods

Line::Line(double k, double b) {
    point_1.x = 0;
    point_1.y = b;
    point_2.x = 1.0;
    point_2.y = k + b;
}

Line::Line(Point A, double k) {
    point_1 = A;
    double b = A.y - k * A.x;
    point_2.x = point_1.x + 1.0;
    point_2.y = k * point_2.x + b;
}

bool Line::IsHorizontal() const {
    return EqualWithEps(point_2.y, point_1.y);
}

bool Line::IsVertical() const {
    return EqualWithEps(point_1.x, point_2.x);
}

bool Line::operator==(const Line &other) const {
    return other.ContainPoint(point_1) && other.ContainPoint(point_2);
}

bool Line::operator!=(const Line &other) const {
    return !(*this == other);
}

void Line::Rotate(Point centre, double angle) {
    point_1.Rotate(centre, angle);
    point_2.Rotate(centre, angle);
}

void Line::GetNormalEquation(double &a, double &b, double &c) const {
    if (IsVertical()) {
        b = 0.0;
        a = 1.0;
        c = -point_1.x;
    } else if (IsHorizontal()) {
        a = 0.0;
        b = 1.0;
        c = -point_1.y;
    } else {
        a = (point_1.y - point_2.y) / (point_1.x - point_2.x);
        b = 1.0;
        c = -(point_1.y - point_1.x * a);
        a *= -1.0;
    }
}

bool Line::ContainPoint(Point point) const {
    double a, b, c, d;
    GetNormalEquation(a, b, c);
    d = a * point.x + b * point.y + c;
    return EqualWithEps(d, 0.0);
}

bool Line::IsAbove(Point point) const {
    if (IsVertical()) {
        return point.x >= point_1.x;
    }
    double aa, bb, cc;
    GetNormalEquation(aa, bb, cc);
    double put = aa * point.x + bb * point.y + cc;
    return put >= 0;
}

bool Line::IsBelow(Point point) const {
    if (IsVertical()) {
        return point.x <= point_1.x;
    }
    double aa, bb, cc;
    GetNormalEquation(aa, bb, cc);
    double put = aa * point.x + bb * point.y + cc;
    return put <= 0;
}

std::pair<Point, Point> Line::GetTwoPoints() const {
    return std::make_pair(point_1, point_2);
}

/// GeomVector Methods

GeomVector::GeomVector(Point A, Point B) {
    x = B.x - A.x;
    y = B.y - A.y;
}

double GeomVector::Angle(GeomVector other) const {
    double module = sqrt(x * x + y * y);
    double module_other = sqrt(other.x * other.x + other.y * other.y);
    double scalar = other.x * x + other.y * y;
    scalar /= (module * module_other);
    return acos(scalar);
}

double GeomVector::CosineAngle(GeomVector other) const {
    double module = x * x + y * y;
    double module_other = other.x * other.x + other.y * other.y;
    double scalar = other.x * x + other.y * y;
    return (scalar / (module * module_other)) * scalar;
}

double GeomVector::SinePseudoScalar(GeomVector other) const {
    double module = sqrt(x * x + y * y);
    double module_other = sqrt(other.x * other.x + other.y * other.y);
    double scalar = x * other.y - other.x * y;
    return scalar / (module * module_other);
}

/// Circle Methods

Point Circle::GetCentre() const {
    return focus_1;
}

Circle::Circle(Point centre, double radius) {
    focus_2 = focus_1 = centre;
    b = a = radius;
}

double Circle::Radius() const {
    return a;
}

double Circle::Perimeter() const {
    return (2 * PI * a);
}

double Circle::Area() const {
    return (PI * a * a);
}

void Circle::Scale(Point centre, double coefficient) {
    focus_1.Scale(centre, coefficient);
    focus_2.Scale(centre, coefficient);
    a *= fabs(coefficient);
    b *= fabs(coefficient);
}

void Circle::Rotate(Point centre, double angle) {
    focus_1.Rotate(centre, angle);
    focus_2.Rotate(centre, angle);
}

void Circle::Reflect(Line axis) {
    focus_1.Reflect(axis);
    focus_2.Reflect(axis);
}

void Circle::Reflect(Point centre) {
    focus_1.Reflect(centre);
    focus_2.Reflect(centre);
}

bool Circle::ContainsPoint(Point point) const {
    double d = Dist(point, focus_1);
    return d <= a + EPS;
}

bool Circle::IsSimilarTo(const Shape &another) const {
    if (typeid(another) == typeid(Circle)) {
        return true;
    } else if (typeid(another) == typeid(Ellipse)) {
        const Ellipse &other = dynamic_cast<const Ellipse &>(another);
        return other.IsCircle();
    }
    return false;
}

bool Circle::IsCongruentTo(const Shape &another) const {
    if (typeid(another) == typeid(Circle)) {
        const Circle &other = dynamic_cast<const Circle & >( another );
        return EqualWithEps(other.a, a);
    } else if (typeid(another) == typeid(Ellipse)) {
        const Ellipse &other = dynamic_cast<const Ellipse &>(another);
        return other.IsCircle() && EqualWithEps(a, other.GetA());
    }
    return false;
}

bool Circle::operator==(const Shape &another) const {
    if (typeid(another) == typeid(Circle)) {
        const Circle &other = dynamic_cast<const Circle & >( another );
        return EqualWithEps(other.a, a) && (other.focus_1 == focus_1);
    } else if (typeid(another) == typeid(Ellipse)) {
        const Ellipse &other = dynamic_cast<const Ellipse &>(another);
        return other.IsCircle() && EqualWithEps(other.Area(), this->Area()) && (other.GetFocuses().first == focus_1);
    }
    return false;
}

bool Circle::operator!=(const Shape &another) const {
    return !(*this == another);
}

/// Ellipse Methods

double Ellipse::GetA() const {
    return a;
}

Ellipse::Ellipse(Point A, Point B, double a_2) {
    focus_1 = A;
    focus_2 = B;
    a = a_2 / 2.0;
    CalcB();
}

void Ellipse::CalcB() {
    double aa = GetFocalDist();
    b = sqrt((a * a - aa * aa));
}

double Ellipse::GetFocalDist() const {
    return (Dist(focus_1, focus_2) / 2);
}

void Ellipse::Rotate(Point centre, double angle) {
    focus_1.Rotate(centre, angle);
    focus_2.Rotate(centre, angle);
}

void Ellipse::Scale(Point centre, double coefficient) {
    focus_1.Scale(centre, coefficient);
    focus_2.Scale(centre, coefficient);
    a *= fabs(coefficient);
    CalcB();
}

void Ellipse::Reflect(Line axis) {
    focus_1.Reflect(axis);
    focus_2.Reflect(axis);
}

void Ellipse::Reflect(Point centre) {
    focus_1.Reflect(centre);
    focus_2.Reflect(centre);
}

double Ellipse::Perimeter() const {

    double armat = sqrt((3 * a + b) * (a + 3 * b));
    double amboxj = 3 * (a + b) - armat;
    return PI * amboxj;
}

bool Ellipse::IsSimilarTo(const Shape &another) const {
    if (typeid(another) == typeid(Ellipse)) {
        const Ellipse &other = dynamic_cast<const Ellipse &>( another );
        return EqualWithEps(a * other.b, b * other.a);
    } else if (typeid(another) == typeid(Circle)) {
        return IsCircle();
    }
    return false;
}

bool Ellipse::IsCongruentTo(const Shape &another) const {
    if (typeid(another) == typeid(Ellipse)) {
        const Ellipse &other = dynamic_cast<const Ellipse &>( another);
        return EqualWithEps(a, other.a) && EqualWithEps(b, other.b);
    } else if (typeid(another) == typeid(Circle)) {
        const Circle &other = dynamic_cast<const Circle &>( another );
        return (focus_1 == focus_2) && EqualWithEps(a, other.a);
    } else {
        return false;
    }
}

std::pair<Point, Point> Ellipse::GetFocuses() const {
    return std::make_pair(focus_1, focus_2);
}

double Ellipse::Eccentricity() const {
    return (GetFocalDist() / a);
}

Point Ellipse::GetCentre() const {
    return CutCentre(focus_1, focus_2);
}

std::pair<Line, Line> Ellipse::Headmistresses() const {
    double aa, bb, cc;
    Line line(focus_1, focus_2);
    line.GetNormalEquation(aa, bb, cc);
    Point solution;
    double r, d, e;
    e = Eccentricity();
    r = a - GetFocalDist();
    d = r / e;
    double dist = r + d;
    Line first_directrice = GetHeadmistresses(focus_1, aa, bb, cc, dist);
    Line second_directrice = GetHeadmistresses(focus_2, aa, bb, cc, dist);
    return std::make_pair(first_directrice, second_directrice);
}

bool Ellipse::operator==(const Shape &another) const {
    if (typeid(another) == typeid(Ellipse)) {
        const Ellipse &other = dynamic_cast<const Ellipse &>(another);
        if ((focus_1 == other.focus_1 && focus_2 == other.focus_2) ||
            (focus_1 == other.focus_2 && focus_2 == other.focus_1)) {
            return EqualWithEps(a, other.a);
        }
        return false;
    } else if (typeid(another) == typeid(Circle)) {
        const Circle &other = dynamic_cast<const Circle &>( another );
        return (focus_1 == focus_2) && (focus_1 == other.focus_1) && EqualWithEps(a, other.a);
    }
    return false;
}

double Ellipse::Area() const {
    return (PI * a * b);
}

bool Ellipse::ContainsPoint(Point point) const {
    double dist_from_focs = Dist(point, focus_1) + Dist(point, focus_2);
    return (dist_from_focs < (2 * a) + EPS);
}

Line Ellipse::GetHeadmistresses(Point focus, double aa, double bb, double cc, double dist) const {
    Point solution;
    if (EqualWithEps(aa, 0)) {
        solution.y = -cc / bb;
        double tmp = sqrt(dist * dist - (solution.y - focus.y) * (solution.y - focus.y));
        solution.x = tmp + focus.x;
        if (ContainsPoint(solution))
            solution.x = -tmp + focus.x;
    } else if (EqualWithEps(bb, 0)) {
        solution.x = -cc / bb;
        double tmp = sqrt(dist * dist - (solution.x - focus.x) * (solution.x - focus.x));
        solution.y = tmp + focus.y;
        if (ContainsPoint(solution))
            solution.y = -tmp + focus.y;
    } else {
        double aaa = (bb / aa) * (bb / aa) + 1;
        double bbb = 2 * ((bb / aa) * (cc / aa + focus.x) - focus.y);
        double ccc = focus.y * focus.y + (cc / aa + focus.x) * (cc / aa + focus.x) - dist * dist;
        std::pair<double, double> resh = SquareEquation(aaa, bbb, ccc);
        solution.y = resh.first;
        solution.x = -(cc / aa) - (bb / aa) * solution.y;
        if (ContainsPoint(solution)) {
            solution.y = resh.second;
            solution.x = -(cc / aa) - (bb / aa) * solution.y;
        }
    }

    Point second_point(solution.x + aa, solution.y + bb);
    Line line(solution, second_point);
    return line;
}

bool Ellipse::IsCircle() const {
    return (focus_1 == focus_2) && (EqualWithEps(a, b));
}

bool Ellipse::operator!=(const Shape &another) const {
    return !(*this == another);
}

/// Polygon Methods

Polygon::Polygon(const std::vector<Point> &p) {
    points = p;
}

int Polygon::NextIndex(int i) const {
    if (i == points.size() - 1) {
        return 0;
    }
    return i + 1;
}

int Polygon::PreviousIndex(int i) const {
    if (i == 0) {
        return points.size() - 1;
    }
    return i - 1;
}

std::vector<Line> Polygon::GetBorder() const {
    std::vector<Line> lines;
    for (int i = 0; i < points.size(); ++i) {
        lines.emplace_back(Line(points[i], points[NextIndex(i)]));
    }
    return lines;
}

int Polygon::VerticesCount() const {
    return points.size();
}

std::vector<Point> Polygon::GetVertices() const {
    return points;
}

bool Polygon::IsConvex() const {
    std::vector<Line> lines = GetBorder();
    std::vector<bool> upper_lower;
    int count = points.size();
    for (int i = 0; i < lines.size(); ++i) {
        Line current_line = lines[i];
        upper_lower.clear();
        for (int j = 0; j < count; ++j) {
            if (!current_line.ContainPoint(points[j])) {
                upper_lower.emplace_back(current_line.IsAbove(points[j]));
            }
        }
        bool value;
        if (!upper_lower.empty()) {
            value = upper_lower[0];
        } else {
            value = true;
        }
        for (int j = 0; j < upper_lower.size(); ++j) {
            if (upper_lower[j] != value) {
                return false;
            }
        }
    }
    return true;
}

bool Polygon::ContainsPoint(Point point) const {
    for (int i = 0; i < points.size(); ++i) {
        if (CBelongsACutB(point, points[i], points[NextIndex(i)])) {
            return true;
        }
    }

    std::vector<Line> lines = GetBorder();
    Line axis(point, Point(point.x, point.y - 1.0));
    std::vector<bool> answers;
    answers.push_back(ShootRay(point, axis, 0.0, lines));
    answers.push_back(ShootRay(point, axis, 3.0, lines));
    answers.push_back(ShootRay(point, axis, 6.0, lines));
    answers.push_back(ShootRay(point, axis, 13.0, lines));
    answers.push_back(ShootRay(point, axis, 19.0, lines));
    answers.push_back(ShootRay(point, axis, 29.0, lines));
    answers.push_back(ShootRay(point, axis, 67.0, lines));
    answers.push_back(ShootRay(point, axis, 71.0, lines));
    answers.push_back(ShootRay(point, axis, 41.0, lines));
    int false_count = 0;
    int true_count = 0;
    for (int j = 0; j < answers.size(); ++j) {
        if (answers[j]) {
            ++true_count;
        } else {
            ++false_count;
        }
    }
    return true_count >= false_count;
}

bool Polygon::ShootRay(Point centre, Line ray, double angle, std::vector<Line> &lines) const {
    int intersections_count = 0;
    ray.Rotate(centre, angle);
    Line horizontal(centre, Point(centre.x + 1.0, centre.y));
    std::vector<bool> colour(points.size(), false);
    for (int j = 0; j < lines.size(); ++j) {
        Point H;
        if (LineIntersection(ray, lines[j], H)) {
            std::pair<Point, Point> AB = lines[j].GetTwoPoints();
            if (CBelongsACutB(H, AB.first, AB.second) && horizontal.IsAbove(H)) {
                bool any_coincidence = false;
                for (int i = 0; i < points.size(); ++i) {
                    if (H == points[i]) {
                        if (colour[i]) {
                            any_coincidence = true;
                        } else {
                            colour[i] = true;
                        }
                        break;
                    }
                }
                if (!any_coincidence) {
                    ++intersections_count;
                }
            }
        }
    }
    return (intersections_count % 2 == 1);
}

double Polygon::Area() const {
    double area = 0.0, minimum = 0.0;
    Point a;
    Point b;
    for (int j = 0; j < points.size(); ++j) {
        if (points[j].y < minimum) {
            minimum = points[j].y;
        }
    }
    for (int j = 0; j < points.size(); ++j) {
        a = points[NextIndex(j)];
        b = points[j];
        area += (a.x - b.x) * ((a.y - minimum + b.y - minimum) / 2.0);
    }
    return fabs(area);
}

double Polygon::Perimeter() const {
    double perimeter = 0.0;
    for (int j = 0; j < points.size(); ++j) {
        perimeter += Dist(points[j], points[NextIndex(j)]);
    }
    return perimeter;
}

void Polygon::Scale(Point centre, double coefficient) {
    for (int j = 0; j < points.size(); ++j) {
        points[j].Scale(centre, coefficient);
    }
}

void Polygon::Rotate(Point centre, double angle) {
    for (int j = 0; j < points.size(); ++j) {
        points[j].Rotate(centre, angle);
    }
}

void Polygon::Reflect(Point centre) {
    for (int j = 0; j < points.size(); ++j) {
        points[j].Reflect(centre);
    }
}

void Polygon::Reflect(Line axis) {
    for (int j = 0; j < points.size(); ++j) {
        points[j].Reflect(axis);
    }
}

bool Polygon::operator!=(const Shape &another) const {
    return !(*this == another);
}

double Polygon::FindLargestEdge() const {
    double max = 0, tmp;
    for (int i = 0; i < points.size(); ++i) {
        tmp = Dist(points[i], points[NextIndex(i)]);
        if (tmp > max) {
            max = tmp;
        }
    }
    return max;
}

bool Polygon::AllAnglesEqual(const Polygon &other) const {
    std::vector<double> this_angles;
    std::vector<double> other_angles;
    std::vector<double> reversed_other_angles;
    std::vector<Point> reversed_other = other.points;
    std::reverse(reversed_other.begin(), reversed_other.end());
    for (int i = 0; i < points.size(); ++i) {
        this_angles.emplace_back(SineA(points[i], points[PreviousIndex(i)], points[NextIndex(i)]));
        other_angles.emplace_back(SineA(other.points[i], other.points[PreviousIndex(i)], other.points[NextIndex(i)]));
        reversed_other_angles.emplace_back(
                SineA(reversed_other[i], reversed_other[PreviousIndex(i)], reversed_other[NextIndex(i)]));
    }
    bool reversed, straight;
    for (int shift = 0; shift < points.size(); ++shift) {
        reversed = false;
        straight = false;
        for (int i = 0; i < this_angles.size(); ++i) {
            if (!EqualWithEps(this_angles[i], other_angles[NextShiftedIndex(i, shift, this_angles.size())])) {
                straight = true;
            }
            if (!EqualWithEps(this_angles[i],
                              reversed_other_angles[NextShiftedIndex(i, shift, this_angles.size())])) {
                reversed = true;
            }
        }
        if (!reversed || !straight) {
            return true;
        }
    }
    return false;
}

double Polygon::AreSimilarWithCoefficient(const Polygon &other, double similarity_coefficient) const {
    if (!AllAnglesEqual(other)) {
        return false;
    }
    std::vector<double> this_edges;
    std::vector<double> other_edges;
    std::vector<double> reversed_other_edges;
    std::vector<Point> other_points = other.points;
    std::reverse(other_points.begin(), other_points.end());
    for (int i = 0; i < points.size(); ++i) {
        this_edges.emplace_back(Dist(points[i], points[NextIndex(i)]));
        other_edges.emplace_back(Dist(other.points[i], other.points[NextIndex(i)]));
        reversed_other_edges.emplace_back(Dist(other_points[i], other_points[NextIndex(i)]));
    }
    bool straight, reversed;
    for (int shift = 0; shift < this_edges.size(); ++shift) {
        straight = false;
        reversed = false;
        for (int i = 0; i < this_edges.size(); ++i) {
            if (!EqualWithEps((this_edges[i] / other_edges[NextShiftedIndex(i, shift, this_edges.size())]),
                              similarity_coefficient)) {
                straight = true;
            }
            if (!EqualWithEps(
                    (this_edges[i] / reversed_other_edges[NextShiftedIndex(i, shift, this_edges.size())]),
                    similarity_coefficient)) {
                reversed = true;
            }
        }
        if (!straight || !reversed) {
            return true;
        }
    }
    return false;
}

bool Polygon::IsSimilarTo(const Shape &another) const {
    try {
        const Polygon &other = dynamic_cast<const Polygon &>( another );
        if (other.VerticesCount() == VerticesCount()) {
            double this_max_edge = FindLargestEdge();
            double other_max_edge = other.FindLargestEdge();
            return AreSimilarWithCoefficient(other, this_max_edge / other_max_edge);
        }
        return false;
    }
    catch (...) {
        return false;
    }
}

bool Polygon::IsCongruentTo(const Shape &another) const {
    try {
        const Polygon &other = dynamic_cast<const Polygon &>(another);
        if (other.VerticesCount() == VerticesCount())
            return AreSimilarWithCoefficient(other, 1.0);
        return false;
    }
    catch (...) {
        return false;
    }
}


bool Polygon::IsAnyShiftWhenPointsCoincided(const std::vector<Point> &other_points) const {
    bool any_dismatch = false;
    int count = points.size();
    for (int i = 0; i < count; ++i) {
        any_dismatch = false;
        for (int j = 0; j < count; ++j) {
            Point tmp = other_points[NextShiftedIndex(j, i, count)];
            if (tmp != points[j]) {
                any_dismatch = true;
                break;
            }
        }
        if (!any_dismatch) {
            return true;
        }
    }
    return false;
}

bool Polygon::operator==(const Shape &another) const {
    try {
        const Polygon &other = dynamic_cast<const Polygon &>(another);
        if (other.VerticesCount() == VerticesCount()) {
            std::vector<Point> other_points = other.GetVertices();
            bool first_wave = IsAnyShiftWhenPointsCoincided(other_points);
            std::reverse(other_points.begin(), other_points.end());
            bool second_wave = IsAnyShiftWhenPointsCoincided(other_points);
            return first_wave || second_wave;
        }
        return false;
    }
    catch (...) {
        return false;
    }
}

/// Rectangle Methods

Rectangle::Rectangle(Point first, Point second, double relation) : Polygon(first) {
    if (relation < 1.0)
        relation = 1.0 / relation;
    double diameter = Dist(first, second);
    double shorter_edge = diameter / sqrt(1.0 + relation * relation);
    double m = shorter_edge * shorter_edge / diameter;
    Line diagonal(first, second);
    std::pair<Point, Point> f = FindTwoPointsOnLineWithDFromA(first, diagonal, m);
    Point F = CBelongsACutB(f.first, first, second) ? f.first : f.second;
    double a, b, c;
    diagonal.GetNormalEquation(a, b, c);
    Point another = F;
    another += Point(a, b);
    Line normal_from_F_to_diagonal(F, another);
    double normal = sqrt(shorter_edge * shorter_edge - m * m);
    std::pair<Point, Point> h = FindTwoPointsOnLineWithDFromA(F, normal_from_F_to_diagonal, normal);
    Point H;
    H = (h.first.x < first.x || EqualWithEps(h.first.x, first.x)) ? h.first : h.second;
    Point L = CutCentre(first, second);
    Point J = RecoverSecondVertexOfCut(H, L);
    points.emplace_back(H);
    points.emplace_back(second);
    points.emplace_back(J);
}

std::pair<double, double> Rectangle::GetEdges() const {
    return std::make_pair(Dist(points[0], points[1]), Dist(points[1], points[2]));
}

int Rectangle::NextIndex(int i) const {
    return (i == 3) ? 0 : (i + 1);
}

Point Rectangle::Centre() const {
    Point L = CutCentre(points[0], points[2]);
    return L;
}

std::pair<Line, Line> Rectangle::Diagonals() const {
    Line first_diagonal(points[0], points[2]);
    Line second_diagonal(points[1], points[3]);
    return std::make_pair(first_diagonal, second_diagonal);
}

double Rectangle::Area() const {
    std::pair<double, double> edges = GetEdges();
    return edges.first * edges.second;
}

double Rectangle::Perimeter() const {
    std::pair<double, double> edges = GetEdges();
    return 2 * (edges.second + edges.first);
}

bool Rectangle::ContainsPoint(Point point) const {
    for (int j = 0; j < 4; ++j)
        if (CBelongsACutB(point, points[j], points[NextIndex(j)]))
            return true;
    std::vector<GeomVector> vectors;
    for (int i = 0; i < 4; ++i)
        vectors.push_back(GeomVector(point, points[i]));
    double angle_sum = 0;
    for (int j = 0; j < 4; ++j)
        angle_sum += vectors[j].Angle(vectors[NextIndex(j)]);
    return EqualWithEps(2 * PI, angle_sum);
}

/// Square Methods

int Square::NextIndex(int i) const {
    if (i == 3) {
        return 0;
    }
    return i + 1;
}

double Square::GetEdgeLen() const {
    return Dist(points[0], points[1]);
}

Circle Square::CircumscribedCircle() const {
    double radius = (GetEdgeLen() / sqrt(2));
    Point centre = CutCentre(points[0], points[2]);
    return Circle(centre, radius);
}

Circle Square::InscribedCircle() const {
    double radius = (GetEdgeLen() / 2);
    Point centre = CutCentre(points[0], points[2]);
    return Circle(centre, radius);
}

double Square::Perimeter() const {
    double edge = GetEdgeLen();
    return (4 * edge);
}

double Square::Area() const {
    double edge = GetEdgeLen();
    return (edge * edge);
}

bool Square::ContainsPoint(Point point) const {
    for (int j = 0; j < 4; ++j) {
        if (CBelongsACutB(point, points[j], points[NextIndex(j)])) {
            return true;
        }
    }
    std::vector<GeomVector> vectors;
    for (int i = 0; i < 4; ++i) {
        vectors.push_back(GeomVector(point, points[i]));
    }
    double angle_sum = 0;
    for (int j = 0; j < 4; ++j) {
        angle_sum += vectors[j].Angle(vectors[NextIndex(j)]);
    }
    return EqualWithEps(2 * PI, angle_sum);
}

/// Triangle Methods

int Triangle::NextIndex(int i) const {
    if (i == 2) {
        return 0;
    }
    return i + 1;
}

std::vector<double> Triangle::GetEdges() const {
    std::vector<double> edges(3);
    edges[0] = Dist(points[0], points[1]);
    edges[1] = Dist(points[1], points[2]);
    edges[2] = Dist(points[2], points[0]);
    return edges;
}

double Triangle::CircumscribedRadius() const {
    std::vector<double> edges = GetEdges();
    double area_s = Area();
    double radius = (edges[0] * edges[1] * edges[2]) / (4.0 * area_s);
    return radius;
}

double Triangle::InscribedRadius() const {
    double p = Perimeter();
    p /= 2.0;
    double s = Area();
    double radius = s / p;
    return radius;
}

Circle Triangle::CircumscribedCircle() const {
    Point A = points[0], B = points[1], C = points[2];
    Point M = CutCentre(A, C);
    Point N = CutCentre(A, B);
    Line A_line_C(A, C);
    Line A_line_B(A, B);
    double a, b, c;
    Point M_1 = M, N_1 = N;
    A_line_C.GetNormalEquation(a, b, c);
    M_1 += Point(a, b);
    A_line_B.GetNormalEquation(a, b, c);
    N_1 += Point(a, b);
    Line mid_perpend_M(M, M_1);
    Line mid_perpend_N(N, N_1);
    Point O_centre;
    bool has = LineIntersection(mid_perpend_M, mid_perpend_N, O_centre);
    double radius = CircumscribedRadius();
    return Circle(O_centre, radius);
}

Circle Triangle::InscribedCircle() const {
    Point A = points[0], B = points[1], C = points[2];
    double p = Perimeter();
    p /= 2.0;
    double x = p - Dist(B, C);
    Line A_line_B(A, B);
    std::pair<Point, Point> two_points = FindTwoPointsOnLineWithDFromA(A, A_line_B, x);
    Point K;
    K = (ContainsPoint(two_points.first)) ? two_points.first : two_points.second;
    double a, b, c;
    A_line_B.GetNormalEquation(a, b, c);
    Point F = K;
    F += Point(a, b);
    Line K_line_I(K, F);
    double radius = InscribedRadius();
    two_points = FindTwoPointsOnLineWithDFromA(K, K_line_I, radius);
    Point centre = (ContainsPoint(two_points.first)) ? two_points.first : two_points.second;
    return Circle(centre, radius);
}

Point Triangle::Centroid() const {
    Point A = points[0], B = points[1], C = points[2];
    Point P = CutCentre(A, C);
    Point Q = CutCentre(B, C);
    Line A_line_Q(A, Q);
    Line B_line_p(B, P);
    Point M;
    bool has = LineIntersection(A_line_Q, B_line_p, M);
    return M;
}

Point Triangle::Orthocentre() const {
    Point A = points[0], B = points[1], C = points[2];
    Line A_line_C(A, C);
    Line B_line_C(B, C);
    Point D = FindNormalsFundament(B, A_line_C);
    Point K = FindNormalsFundament(A, B_line_C);
    Line D_line_B(D, B);
    Line K_line_A(K, A);
    Point H;
    bool has = LineIntersection(D_line_B, K_line_A, H);
    return H;
}

Line Triangle::EulerLine() const {
    Point M = Centroid();
    Point H = Orthocentre();
    return Line(M, H);
}

Circle Triangle::NinePointsCircle() const {
    Point A = points[0], B = points[1], C = points[2];
    Point A_1 = CutCentre(B, C);
    Point B_1 = CutCentre(A, C);
    Point C_1 = CutCentre(A, B);
    Triangle medians_triangle(A_1, B_1, C_1);
    return medians_triangle.CircumscribedCircle();
}

bool Triangle::ContainsPoint(Point point) const {
    for (int j = 0; j < 3; ++j) {
        if (CBelongsACutB(point, points[j], points[NextIndex(j)])) {
            return true;
        }
    }
    std::vector<GeomVector> vectors;
    for (int i = 0; i < 3; ++i) {
        vectors.push_back(GeomVector(point, points[i]));
    }
    double angle_sum = 0;
    for (int j = 0; j < 3; ++j) {
        angle_sum += vectors[j].Angle(vectors[NextIndex(j)]);
    }
    return EqualWithEps(2 * PI, angle_sum);
}

double Triangle::Perimeter() const {
    std::vector<double> edges = GetEdges();
    double perimeter = edges[0] + edges[1] + edges[2];
    return perimeter;
}

double Triangle::Area() const {
    std::vector<double> edges = GetEdges();
    double p = edges[0] + edges[1] + edges[2];
    p /= 2;
    double s_quadrat = p * (p - edges[0]) * (p - edges[1]) * (p - edges[2]);
    return sqrt(s_quadrat);
}

/// Global Functions

double Dist(Point a, Point b) {
    double dist = std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    return dist;
}

Point CutCentre(Point A, Point B) {
    Point c;
    c.x = (A.x + B.x) / 2;
    c.y = (A.y + B.y) / 2;
    return c;
}

std::pair<double, double> SquareEquation(double a, double b, double c) {
    std::pair<double, double> answer;
    double D = b * b - 4 * a * c;
    if (D < 0) {
        answer.first = 0;
        answer.second = 0;
        return answer;
    }
    answer.first = (-b + sqrt(D)) / (2 * a);
    answer.second = (-b - sqrt(D)) / (2 * a);
    return answer;
}

Point FindNormalsFundament(Point point, Line line) {
    double a, b, c, x, y;
    line.GetNormalEquation(a, b, c);
    if (EqualWithEps(a, 0)) {
        y = -c / b;
        x = ((y - point.y) * a) / b + point.x;
    } else if (EqualWithEps(b, 0)) {
        x = -c / a;
        y = ((x - point.x) * b) / a + point.y;
    } else {
        y = (point.y * a - b * ((c / a) + point.x)) / (a + (b * b) / a);
        x = -(c / a) - ((b / a) * y);
    }
    return Point(x, y);
}

bool EqualWithEps(double a, double b) {
    return std::fabs(a - b) <= EPS;
}

bool CBelongsACutB(Point C, Point A, Point B) {
    double dist_A = Dist(C, A);
    double dist_B = Dist(C, B);
    double AB = Dist(A, B);
    double sum = dist_A + dist_B;
    double a, b, c;
    Line line(A, B);
    line.GetNormalEquation(a, b, c);
    double k = a * C.x + b * C.y + c;
    return EqualWithEps(sum, AB) && EqualWithEps(k, 0);
}

bool LineIntersection(Line first, Line second, Point &intersection_point) {
    double a_1, b_1, c_1, a_2, b_2, c_2;
    first.GetNormalEquation(a_1, b_1, c_1);
    second.GetNormalEquation(a_2, b_2, c_2);
    if (EqualWithEps(a_1 * b_2, a_2 * b_1)) {
        return false;
    }
    if (EqualWithEps(a_1, 0.0)) {
        intersection_point.y = -c_1 / b_1;
        intersection_point.x = (-c_2 - b_2 * intersection_point.y) / a_2;
    } else if (EqualWithEps(b_1, 0)) {
        intersection_point.x = -c_1 / a_1;
        intersection_point.y = (-c_2 - a_2 * intersection_point.x) / b_2;
    } else if (EqualWithEps(a_2, 0.0)) {
        intersection_point.y = -c_2 / b_2;
        intersection_point.x = (-c_1 - b_1 * intersection_point.y) / a_1;
    } else if (EqualWithEps(b_2, 0.0)) {
        intersection_point.x = -c_2 / a_2;
        intersection_point.y = (-c_1 - a_1 * intersection_point.x) / b_1;
    } else {
        intersection_point.y = (((c_1 * a_2) / a_1) - c_2) / (b_2 - ((b_1 * a_2) / a_1));
        intersection_point.x = (-c_1 - b_1 * intersection_point.y) / a_1;
    }
    return true;
}

std::pair<Point, Point> FindTwoPointsOnLineWithDFromA(Point A, Line line, double d) {
    double a, b, c;
    line.GetNormalEquation(a, b, c);
    return SystemSolver(a, b, c, A.x, A.y, d);
}

std::pair<Point, Point> SystemSolver(double a, double b, double c, double x_0, double y_0, double d) {
    Point first;
    Point second;
    if (EqualWithEps(a, 0)) {
        first.y = -c / b;
        second.y = -c / b;
        double s = sqrt((d * d) - (y_0 + (c / b)) * (y_0 + (c / b)));
        first.x = x_0 - s;
        second.x = x_0 + s;
    } else if (EqualWithEps(b, 0)) {
        first.x = -c / a;
        second.x = -c / a;
        double s = sqrt((d * d) - (x_0 - first.x) * (x_0 - first.x));
        first.y = y_0 - s;
        second.y = y_0 + s;
    } else {
        double aa, bb, cc;
        aa = 1 + (b / a) * (b / a);
        bb = 2 * ((b / a) * (x_0 + (c / a)) - y_0);
        cc = (x_0 + (c / a)) * (x_0 + (c / a)) + y_0 * y_0 - d * d;
        std::pair<double, double> y = SquareEquation(aa, bb, cc);
        first.y = y.first;
        second.y = y.second;
        first.x = (-b / a) * first.y - (c / a);
        second.x = (-b / a) * second.y - (c / a);
    }
    return std::make_pair(first, second);
}

Point RecoverSecondVertexOfCut(Point A, Point centre) {
    Point B;
    B.x = 2 * centre.x - A.x;
    B.y = 2 * centre.y - A.y;
    return B;
}

int NextShiftedIndex(int i, int shift, int whole_variants_count) {
    return (i + shift) % whole_variants_count;
}

double SineA(Point A, Point B, Point C) {
    Triangle tri(A, B, C);
    double s = tri.Area();
    double a = Dist(A, B);
    double b = Dist(A, C);
    double sin = (2.0 * s) / (a * b);
    return sin;
}

//void Triangle::Rotate(Point centre, double angle) {
//    for (int j = 0; j < 3; ++j)
//        points[j].Rotate(centre, angle);
//}
//
//void Triangle::Reflect(Line axis) {
//    for (int j = 0; j < 3; ++j)
//        points[j].Reflect(axis);
//}
//
//void Triangle::Reflect(Point centre) {
//    for (int j = 0; j < 3; ++j)
//        points[j].Reflect(centre);
//}
//
//void Triangle::Scale(Point centre, double coefficient) {
//    for (int j = 0; j < 3; ++j)
//        points[j].Scale(centre, coefficient);
//}

//void Square::Rotate(Point centre, double angle) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Rotate(centre, angle);
//}
//
//void Square::Reflect(Line axis) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Reflect(axis);
//}
//
//void Square::Reflect(Point centre) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Reflect(centre);
//}
//
//void Square::Scale(Point centre, double coefficient) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Scale(centre, coefficient);
//}

//void Rectangle::Rotate(Point centre, double angle) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Rotate(centre, angle);
//}
//
//void Rectangle::Reflect(Line axis) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Reflect(axis);
//}
//
//void Rectangle::Reflect(Point centre) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Reflect(centre);
//}
//
//void Rectangle::Scale(Point centre, double coefficient) {
//    for (int j = 0; j < 4; ++j)
//        points[j].Scale(centre, coefficient);
//}

