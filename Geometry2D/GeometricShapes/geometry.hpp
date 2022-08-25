#pragma once

#include <typeinfo>
#include <vector>

struct Point;
struct GeomVector;
class Line;
class Shape;
class Ellipse;
class Circle;
class Polygon;
class Rectangle;
class Square;
class Triangle;

struct Coordinates {
    double x;
    double y;

    explicit Coordinates(double a = 0.0, double b = 0.0) : x(a), y(b) {}
};

struct Point : public Coordinates{
    explicit Point(double a = 0.0, double b = 0.0) : Coordinates(a, b) {}
    bool operator==(const Point &other) const;
    bool operator!=(const Point &other) const;
    Point &operator-=(const Point &other);
    Point &operator+=(const Point &other);
    Point &operator*=(double k);

    void Reflect(Point regarding_other);
    void Reflect(Line axis);
    void Rotate(Point centre, double angle);
    void Scale(Point centre, double coefficient); // Homothety
};

struct GeomVector : public Coordinates{
    explicit GeomVector(double a = 0.0, double b = 0.0) : Coordinates(a, b) {}
    GeomVector(Point A, Point B);

    [[nodiscard]] double Angle(GeomVector other) const;
    [[nodiscard]] double CosineAngle(GeomVector other) const;
    [[nodiscard]] double SinePseudoScalar(GeomVector other) const;
};

class Line {
  public:
    Line(Point A, Point B) : point_1(A), point_2(B) {}
    Line(double k, double b); // kx + b
    Line(Point A, double k);

    Line() = default;
    ~Line() = default;

    bool operator==(const Line &other) const;
    bool operator!=(const Line &other) const;

    void Rotate(Point centre, double angle);
    void GetNormalEquation(double &a, double &b, double &c) const;
    [[nodiscard]] bool ContainPoint(Point point) const;
    [[nodiscard]] bool IsAbove(Point point) const;
    [[nodiscard]] bool IsBelow(Point point) const;
    [[nodiscard]] std::pair<Point, Point> GetTwoPoints() const;

  private:
    [[nodiscard]] bool IsVertical() const;
    [[nodiscard]] bool IsHorizontal() const;

  private:
    Point point_1;
    Point point_2;
};

class Shape {
public:
    virtual ~Shape() = default;

    virtual bool operator==(const Shape &another) const = 0;
    virtual bool operator!=(const Shape &another) const = 0;

    virtual void Rotate(Point centre, double angle) = 0;
    virtual void Reflect(Point centre) = 0;
    virtual void Reflect(Line axis) = 0;
    virtual void Scale(Point centre, double coefficient) = 0;

    [[nodiscard]] virtual double Perimeter() const = 0;
    [[nodiscard]] virtual double Area() const = 0;


    [[nodiscard]] virtual bool IsCongruentTo(const Shape &another) const = 0;
    [[nodiscard]] virtual bool IsSimilarTo(const Shape &another) const = 0;
    [[nodiscard]] virtual bool ContainsPoint(Point point) const = 0;
};

class Ellipse : public Shape {
  public:
    Ellipse(Point A, Point B, double a_2);
    Ellipse() = default;

    bool operator==(const Shape &other) const override;
    bool operator!=(const Shape &other) const override;

    void Rotate(Point centre, double angle) override;
    void Reflect(Point centre) override;
    void Reflect(Line axis) override;
    void Scale(Point centre, double coefficient) override;

    [[nodiscard]] double Perimeter() const override;
    [[nodiscard]] double Area() const override;
    [[nodiscard]] bool IsCongruentTo(const Shape &another) const override;
    [[nodiscard]] bool IsSimilarTo(const Shape &another) const override;
    [[nodiscard]] bool ContainsPoint(Point point) const override;


    [[nodiscard]] bool IsCircle() const;
    [[nodiscard]] double GetA() const;
    [[nodiscard]] std::pair<Point, Point> GetFocuses() const;
    [[nodiscard]] virtual Point GetCentre() const;
    [[nodiscard]] std::pair<Line, Line> Headmistresses() const;
    [[nodiscard]] double Eccentricity() const;

  private:
    void CalcB();
    [[nodiscard]] Line GetHeadmistresses(Point focus, double aa, double bb, double cc, double dist) const;
    [[nodiscard]] double GetFocalDist() const;

  protected:
    double a;
    double b;
    Point focus_1;
    Point focus_2;
};

class Circle : public Ellipse {
public:
    Circle(Point centre, double radius);
    Circle() = default;

    bool operator==(const Shape &another) const override;
    bool operator!=(const Shape &another) const override;

    void Rotate(Point centre, double angle) override;
    void Reflect(Point centre) override;
    void Reflect(Line axis) override;
    void Scale(Point centre, double coefficient) override;
    [[nodiscard]] Point GetCentre() const override;
    [[nodiscard]] double Perimeter() const override;
    [[nodiscard]] double Area() const override;
    [[nodiscard]] bool IsCongruentTo(const Shape &another) const override;
    [[nodiscard]] bool IsSimilarTo(const Shape &another) const override;
    [[nodiscard]] bool ContainsPoint(Point point) const override;
    [[nodiscard]] double Radius() const;

};

class Polygon : public Shape {
  public:
    explicit Polygon(const std::vector<Point> &p);

    template<typename T, typename... Args>
    explicit Polygon(T first, Args... args) {
        AddPoints(first, args...);
    }

    bool operator==(const Shape &another) const override;
    bool operator!=(const Shape &another) const override;

    void Rotate(Point centre, double angle) override;
    void Reflect(Point centre) override;
    void Reflect(Line axis) override;
    void Scale(Point centre, double coefficient) override;
    [[nodiscard]] double Perimeter() const override;
    [[nodiscard]] double Area() const override;
    [[nodiscard]] bool IsCongruentTo(const Shape &another) const override;
    [[nodiscard]] bool IsSimilarTo(const Shape &another) const override;
    [[nodiscard]] bool ContainsPoint(Point point) const override;

    [[nodiscard]] int VerticesCount() const;
    [[nodiscard]] std::vector<Point> GetVertices() const;
    [[nodiscard]] bool IsConvex() const;

  private:
    template<typename T>
    void AddPoints(T __x) {
        points.emplace_back(__x);
    }

    template<typename T, typename... Args>
    void AddPoints(const T &first, const Args &... args) {
        AddPoints(first);
        AddPoints(args...);
    }

    [[nodiscard]] int NextIndex(int i) const;
    [[nodiscard]] int PreviousIndex(int i) const;
    [[nodiscard]] std::vector<Line> GetBorder() const;
    [[nodiscard]] bool ShootRay(Point centre, Line ray, double angle, std::vector<Line> &lines) const;
    [[nodiscard]] double FindLargestEdge() const;
    [[nodiscard]] double AreSimilarWithCoefficient(const Polygon &other, double similarity_coefficient) const;
    [[nodiscard]] bool IsAnyShiftWhenPointsCoincided(const std::vector<Point> &other_points) const;
    [[nodiscard]] bool AllAnglesEqual(const Polygon &other) const;

  protected:
    std::vector<Point> points;
};

class Rectangle : public Polygon {
  public:
    Rectangle(Point p_1, Point p_2, Point p_3, Point p_4) : Polygon(p_1, p_2, p_3, p_4) {}
    Rectangle(Point first, Point second, double relation);

    [[nodiscard]] double Perimeter() const override;
    [[nodiscard]] double Area() const override;
    [[nodiscard]] bool ContainsPoint(Point point) const override;

    [[nodiscard]] std::pair<double, double> GetEdges() const;
    [[nodiscard]] Point Centre() const;
    [[nodiscard]] std::pair<Line, Line> Diagonals() const;

  private:
    [[noidscard]] int NextIndex(int i) const;
};

class Square : public Rectangle {
  public:
    Square(Point p_1, Point p_2, Point p_3, Point p_4) : Rectangle(p_1, p_2, p_3, p_4) {}
    Square(Point p_1, Point p_2) : Rectangle(p_1, p_2, 1.0) {}

    [[nodiscard]] double Perimeter() const override;
    [[nodiscard]] double Area() const override;
    [[nodiscard]] bool ContainsPoint(Point point) const override;

    [[nodiscard]] Circle CircumscribedCircle() const;
    [[nodiscard]] Circle InscribedCircle() const;

  private:
    [[nodiscard]] double GetEdgeLen() const;
    [[nodiscard]] int NextIndex(int i) const;

};

//// Triangle

class Triangle : public Polygon {
  public:
    Triangle(Point p_1, Point p_2, Point p_3) : Polygon(p_1, p_2, p_3) {}

    [[nodiscard]] double Perimeter() const override;
    [[nodiscard]] double Area() const override;
    [[nodiscard]] bool ContainsPoint(Point point) const override;
    [[nodiscard]] Circle CircumscribedCircle() const;
    [[nodiscard]] Circle InscribedCircle() const;
    [[nodiscard]] Point Centroid() const;
    [[nodiscard]] Point Orthocentre() const;
    [[nodiscard]] Line EulerLine() const;
    [[nodiscard]] Circle NinePointsCircle() const;
    [[nodiscard]] std::vector<double> GetEdges() const;
    [[nodiscard]] double CircumscribedRadius() const;
    [[nodiscard]] double InscribedRadius() const;

  private:
    [[nodiscard]] int NextIndex(int i) const;

};

std::pair<double, double> SquareEquation(double a, double b, double c);
double Dist(Point, Point);
Point CutCentre(Point A, Point B);
Point FindNormalsFundament(Point point, Line line);
bool EqualWithEps(double a, double b);
bool CBelongsACutB(Point C, Point A, Point B);
bool LineIntersection(Line first, Line second, Point &intersection_point);
std::pair<Point, Point> FindTwoPointsOnLineWithDFromA(Point A, Line line, double d);
std::pair<Point, Point> SystemSolver(double a, double b, double c, double x_0, double y_0, double d);
Point RecoverSecondVertexOfCut(Point A, Point centre);
int NextShiftedIndex(int i, int shift, int whole_variants_count);
double GetSimilarityCoefficient(const Polygon &first, const Polygon &second);
double SineA(Point A, Point B, Point C);

