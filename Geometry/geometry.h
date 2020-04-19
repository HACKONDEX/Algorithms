#ifndef GEOMETRY_GEOMETRY_CPP
#define GEOMETRY_GEOMETRY_CPP

#include <cmath>
#include <iomanip>
#include <vector>
#include <typeinfo>

using std::sqrt;
using std::fabs;

const double PI = std::acos(-1);
const double Eps = 0.0000001;
const double EPS = 0.0000000001;

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

std::pair< double, double > square_equation( double a, double b, double c );
double dist( Point, Point );
Point cut_centre( Point A, Point B );
Point find_normals_fundament( Point point, Line line );
bool equale_with_eps( double a, double b );
bool C_belongs_A_cut_B( Point C, Point A, Point B );
bool line_intersection( Line first, Line second, Point& intersection_point );
std::pair< Point, Point > find_two_points_on_a_line_with_d_from_A( Point A, Line line, double d );
std::pair< Point, Point > specific_system_solver( double a, double b, double c, double x_0, double y_0, double d);
Point recover_second_vertex_of_cut( Point A, Point centre );
int next_shifted_index( int i, int shift, int whole_variants_count);
template < typename T >
void reverse_vector( std::vector<T>& vec)
{
    //// Я тупой балван, слобак!!
    int count = vec.size();
    for( int i = 0; i < count / 2; ++i )
        std::swap(vec[i],vec[count - 1 - i] );
}
double get_similarity_coefficient( const Polygon& first, const Polygon& second );
double sin_A( Point A, Point B, Point C );


//// Point

struct Point
{
    double x;
    double y;

    Point( double a = 0, double b = 0 ) : x( a ), y( b ) {}

    bool operator == ( const Point& other ) const;
    bool operator != ( const Point& other ) const;
    Point& operator -= ( const Point& other );
    Point& operator += ( const Point& other );
    Point& operator *= ( double k );

    void reflex( Point regarding_other );
    void reflex( Line axis );
    void rotate( Point centre, double angle );
    void homothety( Point centre, double coefficient );

};

bool Point::operator==(const Point &other) const
{
    return equale_with_eps(x, other.x) && equale_with_eps(y, other.y);
}

bool Point::operator != ( const Point& other ) const
{
    return !( *this == other );
}

Point& Point::operator-=(const Point &other)
{
    x -= other.x;
    y -= other.y;
    return *this;
}

Point& Point::operator+=(const Point &other)
{
    x += other.x;
    y += other.y;
    return *this;
}

Point& Point::operator*=(double k)
{
    x *= k;
    y *= k;
    return *this;
}

void Point::reflex( Point regarding_other )
{
    x = 2 * regarding_other.x - x;
    y = 2 * regarding_other.y - y;
}

void Point::rotate(Point centre, double angle)
{
    double radian_angle = ( angle * PI ) / double(180.0);
    double sin = std::sin( radian_angle );
    double cos = std::cos( radian_angle );
    //// Change cartesian centre
    Point shifted_point = *this;
    shifted_point -= centre;
    //// Rotate cartesian system -angle
    x = shifted_point.x * cos - shifted_point.y * sin;
    y = shifted_point.x * sin + shifted_point.y * cos;
    //// Recover centre
    *this += centre;
}

void Point::homothety(Point centre, double coefficient)
{
    Point direction_vector = *this;
    direction_vector -= centre;
    --coefficient;
    direction_vector *= coefficient;
    *this += direction_vector;
}

//// Line

class Line
{
public:

    Line( Point A, Point B ) : point_1( A ), point_2( B ) {}
    Line( double k, double b );
    Line( Point A, double k );
    Line() = default;
    ~Line() = default;

    bool operator == ( const Line& other ) const;
    bool operator != ( const Line& other ) const;

    void rotate( Point centre, double angle );
    void get_normal_equation( double& a, double& b, double& c ) const;

    bool contain_point( Point point ) const;
    bool is_above( Point point ) const;
    bool is_below( Point point ) const;
    std::pair< Point, Point > get_two_points() const;

private:

    Point point_1;
    Point point_2;

    bool is_vertical() const;
    bool is_horizontal() const;
};

Line::Line( double k, double b )
{
    point_1.x = 0;
    point_1.y = b;
    point_2.x = 1.0;
    point_2.y = k + b;
}

Line::Line(Point A, double k)
{
    point_1 = A;
    double b = A.y - k * A.x;
    point_2.x = point_1.x + 1.0;
    point_2.y = k * point_2.x + b;
}

bool Line::is_horizontal() const
{
    return equale_with_eps(point_2.y, point_1.y );
}

bool Line::is_vertical() const
{
    return equale_with_eps( point_1.x, point_2.x );
}

bool Line::operator==(const Line &other) const
{
    return other.contain_point(point_1) && other.contain_point(point_2);
}

bool Line::operator!=(const Line &other) const
{
    return !( *this == other );
}

void Line::rotate(Point centre, double angle)
{
    point_1.rotate( centre, angle );
    point_2.rotate( centre, angle );
}

void Line::get_normal_equation(double &a, double &b, double &c) const
{
    if( is_vertical() )
    {
        b = 0;
        a = 1;
        c = -point_1.x;
    }
    else if( is_horizontal() )
    {
        a = 0.0;
        b = 1.0;
        c = -point_1.y;
    }
    else
    {
        a = ( point_1.y - point_2.y ) / ( point_1.x - point_2.x );
        b = 1.0;
        c = -( point_1.y - point_1.x * a );
        a *= -1.0;
    }
}

bool Line::contain_point(Point point) const
{
    double a, b, c, d;
    get_normal_equation(a,b,c);
    d = a * point.x + b * point.y + c;
    return equale_with_eps(d,0);
}

bool Line::is_above(Point point) const
{
    if( is_vertical() )
        return point.x >= point_1.x;
    double aa, bb, cc;
    get_normal_equation(aa, bb, cc);
    double put = aa * point.x + bb * point.y + cc;
    return put >= 0;
}

bool Line::is_below(Point point) const
{
    if( is_vertical() )
        return point.x <= point_1.x;
    double aa, bb, cc;
    get_normal_equation(aa, bb, cc);
    double put = aa * point.x + bb * point.y + cc;
    return put <= 0;
}

std::pair< Point, Point > Line::get_two_points() const
{
    return std::make_pair(point_1, point_2);
}

//// GeomVector

struct GeomVector
{
    double x;
    double y;

    double angle( GeomVector other ) const;
    double cosinusn_of_angle( GeomVector other ) const;
    double sinus_psevdo_scalyar( GeomVector other ) const;

    GeomVector( double a = 0, double b = 0) : x(a), y(b){}
    GeomVector( Point A, Point B);
};

double GeomVector::sinus_psevdo_scalyar(GeomVector other) const
{
    double module = sqrt( x*x + y*y );
    double module_other = sqrt( other.x * other.x + other.y * other.y );
    double scalyar = x * other.y - other.x * y;
    return scalyar / ( module * module_other );
}

double GeomVector::cosinusn_of_angle(GeomVector other) const
{
    double module =  x*x + y*y ;
    double module_other =  other.x * other.x + other.y * other.y ;
    double scalyar = other.x * x + other.y * y;
    return (scalyar / ( module * module_other ) ) * scalyar;
}

GeomVector::GeomVector(Point A, Point B)
{
    x = B.x - A.x;
    y = B.y - A.y;
}

double GeomVector::angle(GeomVector other) const
{
    double module = sqrt( x*x + y*y );
    double module_other = sqrt( other.x * other.x + other.y * other.y );
    double scalyar = other.x * x + other.y * y;
    scalyar /= ( module * module_other );
    return acos( scalyar );
}

//// Point methods

void Point::reflex(Line axis)
{
    Point centre = find_normals_fundament( *this, axis );
    reflex( centre );
}

//// Shape

class Shape {
public:

    virtual ~Shape() = default;

    virtual void rotate( Point centre, double angle ) = 0; /// Поворот вокруг точки
    virtual void reflex( Point centre ) = 0; /// Симметричное отражение относительно точки
    virtual void reflex( Line axis ) = 0; /// Симметричное отражение относительно прямой
    virtual void scale( Point centre, double coefficient ) = 0;/// Homothety

    virtual double perimeter() const = 0; /// Перимитр
    virtual double area() const = 0; /// Площадь
    virtual bool operator ==  ( const Shape& another ) const = 0; /// Совподает ли
    virtual bool operator != ( const Shape& another ) const = 0;
    virtual bool isCongruentTo( const Shape& another) const = 0; /// Равен ли
    virtual bool isSimilarTo( const Shape& another ) const = 0; /// Подобие
    virtual bool containsPoint( Point point ) const = 0; /// Содержит ли точку
};

class Ellipse : public Shape
{
public:
    Ellipse( Point A, Point B, double a_2 );
    Ellipse() = default;
//    ~Ellipse() = default;

    std::pair<Point,Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    virtual Point center() const;

    void rotate( Point centre, double angle ) override;
    void reflex( Point centre ) override;
    void reflex( Line axis ) override;
    void scale( Point centre, double coefficient ) override;

    double perimeter() const override;
    double area() const override;
    bool operator ==  ( const Shape& another ) const override;
    bool operator != (const Shape& another ) const override;
    bool isCongruentTo( const Shape& another) const override;
    bool isSimilarTo( const Shape& another ) const override;
    bool containsPoint( Point point ) const override;


    bool is_circle() const;
    double return_a() const;

protected:
    double a;
    double b;
    Point focus_1;
    Point focus_2;
private:
    void calc_b();
    Line get_directrice( Point focus, double aa, double bb, double cc, double dist ) const;
    double get_focal_dist() const;
};

class Circle : public Ellipse
{
public:

    Circle( Point centre, double radius );
    Circle() = default;
//    ~Circle() = default;

    double radius() const;

    void rotate( Point centre, double angle ) override;
    void reflex( Point centre ) override;
    void reflex( Line axis ) override;
    void scale( Point centre, double coefficient ) override;
    Point center() const override;

    double perimeter() const override;
    double area() const override;
    bool operator ==  ( const Shape& another ) const override;
    bool operator != (const Shape& another ) const override;
    bool isCongruentTo( const Shape& another) const override;
    bool isSimilarTo( const Shape& another ) const override;
    bool containsPoint( Point point ) const override;

private:
};

Point Circle::center() const
{
    return  focus_1;
}

double Ellipse::return_a() const
{
    return a;
}

class Polygon : public Shape
{
public:
    template < typename T, typename... Args >
    explicit Polygon(T first, Args... args);

    explicit Polygon( const std::vector<Point>& p );

    int verticesCount() const;
    std::vector<Point> getVertices() const;
    bool isConvex() const;

    void rotate( Point centre, double angle ) override;
    void reflex( Point centre ) override;
    void reflex( Line axis ) override;
    void scale( Point centre, double coefficient ) override;

    double perimeter() const override;
    double area() const override;
    bool operator ==  ( const Shape& another ) const override;
    bool operator != (const Shape& another ) const override;
    bool isCongruentTo( const Shape& another) const override;
    bool isSimilarTo( const Shape& another ) const override;
    bool containsPoint( Point point ) const override;

protected:
    std::vector<Point> points;


private:
    template < typename T >
    void add_points( T __x );

    template < typename T, typename... Args >
    void add_points( const T& first, const Args&... args );

    int next_index( int i ) const;
    int previous_index( int i ) const;
    std::vector<Line> get_boarder() const;
    bool contain_point_shoot_ray( Point centre, Line ray, double angle, std::vector<Line>& lines ) const;
    double find_the_largest_edge() const;
    double are_similar_with_coefficient( const Polygon& other, double similarity_coefficient ) const;
    bool is_there_any_shift_when_all_points_coincided(const std::vector<Point> &other_points) const;
    bool if_all_angles_are_equal(const Polygon& other) const;
};

//// Polygon

template < typename T >
void Polygon::add_points(T __x)
{
    points.emplace_back( __x );
}

template < typename T, typename... Args >
void Polygon::add_points( const T& first, const Args&... args )
{
    add_points( first );
    add_points(args...);
}

template < typename T, typename... Args >
Polygon::Polygon(T first, Args... args)
{
    add_points( first, args...);
}

Polygon::Polygon(const std::vector<Point> &p)
{
    points = p;
}

int Polygon::next_index(int i) const
{
    if( i == points.size() - 1 )
        return 0;
    return i + 1;
}

int Polygon::previous_index(int i) const
{
    if( i == 0 )
        return points.size() - 1;
    return i - 1;
}

std::vector<Line> Polygon::get_boarder() const
{
    std::vector<Line> lines;
    for( int i = 0; i < points.size(); ++i )
        lines.emplace_back( Line(points[i], points[next_index(i)]));
    return lines;
}

int Polygon::verticesCount() const
{
    return points.size();
}

std::vector<Point> Polygon::getVertices() const
{
    return points;
}

bool Polygon::isConvex() const
{
    std::vector<Line> lines = get_boarder();
    std::vector<bool> upper_lower;
    int count = points.size();
    for( int i = 0; i < lines.size(); ++i )
    {
        Line current_line = lines[i];
        upper_lower.clear();
        for (int j = 0; j < count; ++j)
        {
            if (!current_line.contain_point(points[j]))
                upper_lower.emplace_back(current_line.is_above(points[j]));
        }
        bool value;
        if( !upper_lower.empty() )
            value = upper_lower[0];
        else
            value = true;
        for (int j = 0; j < upper_lower.size(); ++j)
            if (upper_lower[j] != value)
                return false;
    }
    return true;
}

bool Polygon::containsPoint(Point point) const {
    for (int i = 0; i < points.size(); ++i)
        if (C_belongs_A_cut_B(point, points[i], points[next_index(i)]))
            return true;

    std::vector<Line> lines = get_boarder();
    Line axis(point, Point(point.x, point.y - 1.0));
    //// Слабая теория вероятности.
    std::vector<bool> answers;
    answers.push_back(contain_point_shoot_ray(point, axis, 0.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 3.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 6.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 13.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 19.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 29.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 67.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 71.0, lines));
    answers.push_back(contain_point_shoot_ray(point, axis, 41.0, lines));
    int false_count = 0;
    int true_count = 0;
    for( int j = 0; j < answers.size(); ++j )
    {
        if(answers[j])
            ++true_count;
        else
            ++false_count;
    }
    return true_count >= false_count;
}

bool Polygon::contain_point_shoot_ray(Point centre, Line ray, double angle, std::vector<Line>& lines) const
{
    int intersections_count = 0;
    ray.rotate(centre, angle);
    Line horizontal(centre, Point(centre.x + 1.0, centre.y));
    std::vector<bool> colour(points.size(), false);
    for( int j = 0; j < lines.size(); ++j )
    {
        Point H;
        if( line_intersection(ray, lines[j], H))
        {
            std::pair< Point, Point > AB = lines[j].get_two_points();
            if( C_belongs_A_cut_B( H, AB.first, AB.second ) && horizontal.is_above(H) )
            {
                bool any_coincidence = false;
                for( int i = 0; i < points.size(); ++i )
                {
                    if( H == points[i] )
                    {
                        if( colour[i] )
                            any_coincidence = true;
                        else
                            colour[i] = true;
                        break;
                    }
                }
                if( !any_coincidence )
                    ++intersections_count;
            }
        }
    }
    return (intersections_count % 2 == 1);
}

double Polygon::area() const
{
    double area = 0.0, minimum = 0.0;
    Point a;
    Point b;
    for( int j = 0; j < points.size(); ++j )
    {
        if( points[j].y < minimum )
        {
            minimum = points[j].y;
        }
    }
    for( int j = 0; j < points.size(); ++j )
    {
        a = points[next_index(j)];
        b = points[j];
        area += ( a.x - b.x ) * ( ( a.y - minimum + b.y - minimum ) / 2.0 );
    }
    return fabs(area);
}

double Polygon::perimeter() const
{
    double perimeter = 0.0;
    for( int j = 0; j < points.size(); ++j )
        perimeter += dist( points[j], points[next_index(j)]);
    return perimeter;
}

void Polygon::scale(Point centre, double coefficient)
{
    for( int j = 0; j < points.size(); ++j )
        points[j].homothety(centre,coefficient);
}

void Polygon::rotate(Point centre, double angle)
{
    for( int j = 0; j < points.size(); ++j )
        points[j].rotate(centre,angle);
}

void Polygon::reflex(Point centre)
{
    for( int j = 0; j < points.size(); ++j )
        points[j].reflex(centre);
}

void Polygon::reflex(Line axis)
{
    for( int j = 0; j < points.size(); ++j )
        points[j].reflex(axis);
}

//// Rectangle

class Rectangle : public Polygon
{
public:

    Rectangle( Point p_1, Point p_2, Point p_3, Point p_4 ) : Polygon(p_1, p_2, p_3, p_4){}
    Rectangle( Point first, Point second, double relation );

    Point center() const;
    std::pair<Line, Line> diagonals() const;


    void rotate( Point centre, double angle ) override;
    void reflex( Point centre ) override;
    void reflex( Line axis ) override;
    void scale( Point centre, double coefficient ) override;

    double perimeter() const override;
    double area() const override;
    bool containsPoint( Point point ) const override;

    std::pair< double, double > get_edges() const;

private:


    int next_index( int i ) const;
};

Rectangle::Rectangle(Point first, Point second, double relation) : Polygon( first )
{
    if( relation < 1.0 )
        relation = 1.0 / relation;
    double diameter = dist( first, second );
    double shorter_edge = diameter / sqrt( 1.0 + relation * relation );
    double m = shorter_edge * shorter_edge / diameter;
    Line diagonal(first, second);
    std::pair< Point, Point > f = find_two_points_on_a_line_with_d_from_A( first,diagonal, m );
    Point F = C_belongs_A_cut_B( f.first, first, second ) ? f.first : f.second;
    double a,b,c;
    diagonal.get_normal_equation(a,b,c);
    Point another = F;
    another += Point( a, b);
    Line normal_from_F_to_diagonal( F, another );
    double normal = sqrt( shorter_edge * shorter_edge - m * m );
    std::pair< Point, Point > h = find_two_points_on_a_line_with_d_from_A( F, normal_from_F_to_diagonal, normal);
    Point H;
    H = ( h.first.x < first.x || equale_with_eps(h.first.x, first.x ) ) ? h.first : h.second;
    Point L = cut_centre( first, second );
    Point J = recover_second_vertex_of_cut( H, L );
    points.emplace_back( H );
    points.emplace_back( second );
    points.emplace_back( J );
}

std::pair< double, double > Rectangle::get_edges() const
{
    return std::make_pair( dist( points[0], points[1]), dist( points[1], points[2]) );
}

int Rectangle::next_index(int i) const
{
    return ( i == 3 ) ? 0 : (i + 1) ;
}

Point Rectangle::center() const
{
    Point L = cut_centre(points[0], points[2]);
    return L;
}

std::pair<Line, Line> Rectangle::diagonals() const
{
    Line first_diagonal(points[0], points[2]);
    Line second_diagonal( points[1], points[3]);
    return std::make_pair( first_diagonal, second_diagonal );
}

double Rectangle::area() const
{
    std::pair< double , double > edges = get_edges();
    return edges.first * edges.second;
}

double Rectangle::perimeter() const
{
    std::pair< double , double > edges = get_edges();
    return 2 * ( edges.second + edges.first );
}

bool Rectangle::containsPoint(Point point) const
{
    for( int j = 0; j < 4; ++j )
        if( C_belongs_A_cut_B( point, points[j],points[next_index(j)]))
            return true;
    std::vector<GeomVector> vectors;
    for( int i = 0; i < 4; ++i )
        vectors.push_back( GeomVector( point, points[i]) );
    double angle_sum = 0;
    for( int j = 0; j < 4; ++j )
        angle_sum += vectors[j].angle(vectors[next_index(j)]);
    return equale_with_eps( 2 * PI, angle_sum );
}

void Rectangle::rotate(Point centre, double angle)
{
    for( int j = 0; j < 4; ++j)
        points[j].rotate(centre, angle);
}

void Rectangle::reflex(Line axis)
{
    for( int j = 0; j < 4; ++j)
        points[j].reflex(axis);
}

void Rectangle::reflex(Point centre)
{
    for( int j = 0; j < 4; ++j)
        points[j].reflex(centre);
}

void Rectangle::scale(Point centre, double coefficient)
{
    for( int j = 0; j < 4; ++j)
        points[j].homothety(centre, coefficient);
}

//// Square

class Square : public Rectangle
{
public:
    Square( Point p_1, Point p_2, Point p_3, Point p_4 ) : Rectangle( p_1, p_2, p_3, p_4 ) {}
    Square( Point p_1, Point p_2 ) : Rectangle( p_1, p_2, 1.0 ) {}

    Circle circumscribedCircle() const;/// Описанная окружность
    Circle inscribedCircle() const;/// Вписанная окружнотсь

    void rotate( Point centre, double angle ) override;
    void reflex( Point centre ) override;
    void reflex( Line axis ) override;
    void scale( Point centre, double coefficient ) override;

    double perimeter() const override;
    double area() const override;
    bool containsPoint( Point point ) const override;

private:

    double get_edge() const;
    int next_index(int i) const;
};

int Square::next_index(int i) const
{
    if( i == 3)
        return 0;
    return i + 1;
}

double Square::get_edge() const
{
    return dist(points[0], points[1]);
}

Circle Square::circumscribedCircle() const
{
    double radius = ( get_edge() / sqrt(2) );
    Point centre = cut_centre(points[0],points[2]);
    return Circle(centre, radius);
}

Circle Square::inscribedCircle() const
{
    double radius = ( get_edge() / 2 );
    Point centre = cut_centre( points[0], points[2]);
    return Circle( centre, radius );
}

double Square::perimeter() const
{
    double edge = get_edge();
    return ( 4 * edge );
}

double Square::area() const
{
    double edge = get_edge();
    return  ( edge * edge );
}

bool Square::containsPoint(Point point) const
{
    for( int j = 0; j < 4; ++j )
        if( C_belongs_A_cut_B(point, points[j], points[next_index(j)]))
            return true;
    std::vector<GeomVector> vectors;
    for( int i = 0; i < 4; ++i )
        vectors.push_back( GeomVector( point, points[i]) );
    double angle_sum = 0;
    for( int j = 0; j < 4; ++j )
        angle_sum += vectors[j].angle(vectors[next_index(j)]);
    return equale_with_eps( 2 * PI, angle_sum );
}

void Square::rotate(Point centre, double angle)
{
    for( int j = 0; j < 4; ++j )
        points[j].rotate(centre, angle);
}

void Square::reflex(Line axis)
{
    for( int j = 0; j < 4; ++j )
        points[j].reflex(axis);
}

void Square::reflex(Point centre)
{
    for( int j = 0; j < 4; ++j )
        points[j].reflex(centre);
}

void Square::scale(Point centre, double coefficient)
{
    for( int j = 0; j < 4; ++j )
        points[j].homothety(centre, coefficient);
}

//// Triangle

class Triangle : public Polygon
{
public:

    Triangle( Point p_1, Point p_2, Point p_3 ) : Polygon( p_1, p_2, p_3 ) {}

    Circle circumscribedCircle() const; /// Описанная окружнсть
    Circle inscribedCircle() const; /// Вписанная окружность
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;

    void rotate( Point centre, double angle ) override;
    void reflex( Point centre ) override;
    void reflex( Line axis ) override;
    void scale( Point centre, double coefficient ) override;

    double perimeter() const override;
    double area() const override;
    bool containsPoint( Point point ) const override;

    std::vector< double > get_edges() const;
    double circumscribed_radius() const;
    double inscribed_radius() const;

private:

    int next_index( int i ) const;
};

int Triangle::next_index(int i) const
{
    if( i == 2 )
        return 0;
    return i + 1;
}

std::vector< double > Triangle::get_edges() const
{
    std::vector< double > edges(3);
    edges[0] = dist(points[0], points[1]);
    edges[1] = dist(points[1], points[2]);
    edges[2] = dist(points[2], points[0]);
    return edges;
}

double Triangle::circumscribed_radius() const
{
    std::vector< double > edges = get_edges();
    double area_s = area();
    double radius = ( edges[0] * edges[1] * edges[2] ) / ( 4.0 * area_s );
    return radius;
}

double Triangle::inscribed_radius() const
{
    double p = perimeter();
    p /= 2.0;
    double s = area();
    double radius = s / p;
    return radius;
}

Circle Triangle::circumscribedCircle() const
{
    Point A = points[0], B = points[1], C = points[2];
    Point M = cut_centre(A, C);
    Point N = cut_centre(A, B);
    Line A_line_C( A, C );
    Line A_line_B( A, B );
    double a, b, c;
    Point M_1 = M, N_1 = N;
    A_line_C.get_normal_equation( a, b, c );
    M_1 += Point( a, b );
    A_line_B.get_normal_equation( a, b, c );
    N_1 += Point( a, b );
    Line mid_perpend_M( M, M_1 );
    Line mid_perpend_N( N, N_1 );
    Point O_centre;
    bool has = line_intersection(mid_perpend_M, mid_perpend_N, O_centre);
    double radius = circumscribed_radius();
    return Circle(O_centre, radius);
}

Circle Triangle::inscribedCircle() const
{
    Point A = points[0], B = points[1], C = points[2];
    double p = perimeter();
    p /= 2.0;
    double x = p - dist(B,C);
    Line A_line_B( A, B );
    std::pair< Point, Point > two_points = find_two_points_on_a_line_with_d_from_A(A, A_line_B, x);
    Point K;
    K = ( containsPoint( two_points.first ) ) ? two_points.first : two_points.second;
    double a, b, c;
    A_line_B.get_normal_equation( a, b, c );
    Point F = K;
    F += Point(a, b);
    Line K_line_I( K, F );
    double radius = inscribed_radius();
    two_points = find_two_points_on_a_line_with_d_from_A(K, K_line_I, radius);
    Point centre = ( containsPoint(two_points.first) ) ? two_points.first : two_points.second;
    return Circle(centre, radius);
}

Point Triangle::centroid() const
{
    Point A = points[0], B = points[1], C = points[2];
    Point P = cut_centre( A, C );
    Point Q = cut_centre( B, C );
    Line A_line_Q( A, Q );
    Line B_line_p( B, P );
    Point M;
    bool has = line_intersection(A_line_Q, B_line_p, M);
    return M;
}

Point Triangle::orthocenter() const
{
    Point A = points[0], B = points[1], C = points[2];
    Line A_line_C( A, C );
    Line B_line_C( B, C );
    Point D = find_normals_fundament(B, A_line_C);
    Point K = find_normals_fundament(A, B_line_C);
    Line D_line_B( D, B );
    Line K_line_A( K, A );
    Point H;
    bool has = line_intersection( D_line_B, K_line_A, H );
    return H;
}

Line Triangle::EulerLine() const
{
    Point M = centroid();
    Point H = orthocenter();
    return Line( M, H);
}

Circle Triangle::ninePointsCircle() const
{
    Point A = points[0], B = points[1], C = points[2];
    Point A_1 = cut_centre( B, C );
    Point B_1 = cut_centre( A, C );
    Point C_1 = cut_centre( A, B );
    Triangle medians_triangle( A_1, B_1, C_1 );
    return medians_triangle.circumscribedCircle();
}

bool Triangle::containsPoint(Point point) const
{
    for( int j = 0; j < 3; ++j )
        if( C_belongs_A_cut_B(point, points[j], points[next_index(j)]))
            return true;
    std::vector<GeomVector> vectors;
    for( int i = 0; i < 3; ++i )
        vectors.push_back( GeomVector( point, points[i]) );
    double angle_sum = 0;
    for( int j = 0; j < 3; ++j )
        angle_sum += vectors[j].angle(vectors[next_index(j)]);
    return equale_with_eps( 2 * PI, angle_sum );
}

double Triangle::perimeter() const
{
    std::vector< double > edges = get_edges();
    double perimeter = edges[0] + edges[1] + edges[2];
    return perimeter;
}

double Triangle::area() const
{
    std::vector< double > edges = get_edges();
    double p = edges[0] + edges[1] + edges[2];
    p /= 2;
    double s_quadrat = p * ( p - edges[0] ) * ( p - edges[1] ) * ( p - edges[2] );
    return sqrt(s_quadrat);
}

void Triangle::rotate(Point centre, double angle)
{
    for( int j = 0; j < 3; ++j )
        points[j].rotate(centre, angle);
}

void Triangle::reflex(Line axis)
{
    for( int j = 0; j < 3; ++j )
        points[j].reflex(axis);
}

void Triangle::reflex(Point centre)
{
    for( int j = 0; j < 3; ++j )
        points[j].reflex(centre);
}

void Triangle::scale(Point centre, double coefficient)
{
    for( int j = 0; j < 3; ++j )
        points[j].homothety(centre, coefficient);
}

//// Ellipse

Ellipse::Ellipse( Point A, Point B, double a_2 )
{
    focus_1 = A;
    focus_2 = B;
    a = a_2 / 2.0;
    calc_b();
}

void Ellipse::calc_b()
{
    double aa = get_focal_dist();
    b = sqrt( ( a * a - aa * aa ) );
}

double Ellipse::get_focal_dist() const
{
    return ( dist( focus_1, focus_2 ) / 2 );
}

void Ellipse::rotate(Point centre, double angle)
{
    focus_1.rotate( centre, angle );
    focus_2.rotate( centre, angle );
}

void Ellipse::scale(Point centre, double coefficient)
{
    focus_1.homothety( centre, coefficient );
    focus_2.homothety( centre, coefficient );
    a *= fabs( coefficient );
    calc_b();
}

void Ellipse::reflex(Line axis)
{
    focus_1.reflex( axis );
    focus_2.reflex( axis );
}

void Ellipse::reflex(Point centre)
{
    focus_1.reflex( centre );
    focus_2.reflex( centre );
}

double Ellipse::perimeter() const
{

    double armat = sqrt( ( 3 * a + b ) * ( a + 3 * b ) );
    double amboxj = 3 * ( a + b ) - armat;
    return PI * amboxj;
}

bool Ellipse::isSimilarTo(const Shape &another) const //// Подобие
{
    if( typeid(another) == typeid(Ellipse) )
    {
        const Ellipse& other = dynamic_cast<const Ellipse&>( another );
        return equale_with_eps( a * other.b, b * other.a );
    }
    else if( typeid(another) == typeid(Circle) )
        return is_circle();
    return false;
}

bool Ellipse::isCongruentTo(const Shape &another) const //// Геометрическое равенство
{
    if( typeid(another) == typeid(Ellipse) )
    {
        const Ellipse& other = dynamic_cast<const Ellipse&>( another);
        return equale_with_eps(a, other.a ) &&  equale_with_eps(b, other.b );
    }
    else if(typeid(another) == typeid(Circle) )
    {
        const Circle& other = dynamic_cast<const Circle&>( another );
        return ( focus_1 == focus_2 ) && equale_with_eps(a, other.a);
    }
    else
        return false;
}

std::pair<Point,Point> Ellipse::focuses() const
{
    return std::make_pair( focus_1, focus_2 );
}

double Ellipse::eccentricity() const
{
    return ( get_focal_dist() / a );
}

Point Ellipse::center() const
{
    return cut_centre( focus_1, focus_2 );
}

std::pair<Line, Line> Ellipse::directrices() const
{
    double aa, bb, cc;
    Line line( focus_1, focus_2 );
    line.get_normal_equation( aa, bb, cc );
    Point solution;
    double r, d, e;
    e = eccentricity();
    r = a - get_focal_dist();
    d = r / e;
    double dist = r + d;
    Line first_directrice = get_directrice( focus_1, aa, bb, cc, dist );
    Line second_directrice = get_directrice(focus_2, aa, bb, cc, dist );
    return std::make_pair( first_directrice, second_directrice );
}

bool Ellipse::operator==(const Shape &another) const
{
    if( typeid(another) == typeid(Ellipse) )
    {
        const Ellipse& other = dynamic_cast<const Ellipse&>(another);
        if( (focus_1 == other.focus_1 && focus_2 == other.focus_2) || (focus_1 == other.focus_2 && focus_2 == other.focus_1) )
            return equale_with_eps(a, other.a);
        return false;
    }
    else if( typeid(another) == typeid(Circle) )
    {
        const Circle& other = dynamic_cast<const Circle&>( another );
        return ( focus_1 == focus_2) && ( focus_1 == other.focus_1 ) && equale_with_eps( a, other.a );
    }
    return false;
}

double Ellipse::area() const
{
    return ( PI * a * b );
}

bool Ellipse::containsPoint(Point point) const
{
    double dist_from_focs = dist( point, focus_1 ) + dist( point, focus_2 );
    return ( dist_from_focs < ( 2 * a ) + Eps );
}

Line Ellipse::get_directrice(Point focus, double aa, double bb, double cc, double dist) const
{
    Point solution;
    if( equale_with_eps( aa, 0 ) )
    {
        solution.y = -cc / bb;
        double tmp = sqrt( dist * dist - ( solution.y - focus.y ) * ( solution.y - focus.y ) );
        solution.x =  tmp + focus.x;
        if( containsPoint( solution ) )
            solution.x = - tmp + focus.x;
    }
    else if( equale_with_eps( bb, 0 ) )
    {
        solution.x = -cc / bb;
        double tmp = sqrt( dist * dist - ( solution.x - focus.x ) * ( solution.x - focus.x ) );
        solution.y = tmp + focus.y;
        if( containsPoint( solution ) )
            solution.y = -tmp + focus.y;
    }
    else
    {
        double aaa = ( bb / aa ) * ( bb / aa ) + 1;
        double bbb = 2 * ( ( bb / aa ) * ( cc / aa + focus.x ) - focus.y );
        double ccc = focus.y * focus.y + ( cc / aa + focus.x ) * ( cc / aa + focus.x ) - dist * dist;
        std::pair< double, double > resh = square_equation( aaa, bbb, ccc );
        solution.y = resh.first;
        solution.x = - ( cc / aa ) - ( bb / aa ) * solution.y;
        if( containsPoint( solution ) )
        {
            solution.y = resh.second;
            solution.x = - ( cc / aa ) - ( bb / aa ) * solution.y;
        }
    }

    Point second_point( solution.x + aa, solution.y + bb );
    Line line( solution, second_point );
    return line;
}

bool Ellipse::is_circle() const
{
    return (focus_1 == focus_2) && ( equale_with_eps( a, b)) ;
}

//// Circle

Circle::Circle( Point centre, double radius )
{
    focus_2 = focus_1 = centre;
    b = a = radius;
}

double Circle::radius() const
{
    return a;
}

double Circle::perimeter() const
{
    return ( 2 * PI * a );
}

double Circle::area() const
{
    return ( PI * a * a );
}

void Circle::scale(Point centre, double coefficient)
{
    focus_1.homothety(centre, coefficient);
    focus_2.homothety( centre, coefficient );
    a *= fabs(coefficient);
    b *= fabs(coefficient);
}

void Circle::rotate(Point centre, double angle)
{
    focus_1.rotate( centre, angle );
    focus_2.rotate(centre, angle);
}

void Circle::reflex(Line axis)
{
    focus_1.reflex( axis );
    focus_2.reflex( axis );
}

void Circle::reflex(Point centre)
{
    focus_1.reflex( centre );
    focus_2.reflex( centre );
}

bool Circle::containsPoint(Point point) const
{
    double d = dist( point, focus_1 );
    return d <= a + Eps;
}

bool Circle::isSimilarTo(const Shape &another) const ////Подобие
{
    if(typeid(another) == typeid(Circle))
        return true;
    else if(typeid(another) == typeid(Ellipse) )
    {
        const Ellipse& other = dynamic_cast<const Ellipse&>(another);
        return other.is_circle();
    }
    return false;
}

bool Circle::isCongruentTo(const Shape &another) const
{
    if(typeid(another) == typeid(Circle))
    {
        const Circle& other = dynamic_cast<const  Circle& >( another );
        return equale_with_eps( other.a, a );
    }
    else if(typeid(another) == typeid(Ellipse))
    {
        const Ellipse& other = dynamic_cast<const Ellipse&>(another);
        return other.is_circle() && equale_with_eps( a, other.return_a() );
    }
    return false;
}

bool Circle::operator==(const Shape &another) const
{
    if(typeid(another) == typeid(Circle))
    {
        const Circle& other = dynamic_cast<const  Circle& >( another );
        return equale_with_eps( other.a, a ) && ( other.focus_1 == focus_1 );
    }
    else if(typeid(another) == typeid(Ellipse))
    {
        const Ellipse& other = dynamic_cast<const Ellipse&>(another);
        return other.is_circle() && equale_with_eps( other.area(), this->area() ) && (other.focuses().first == focus_1);
    }
    return false;
}

//// Similarity of Polygon

double Polygon::find_the_largest_edge() const
{
    double max = 0, tmp;
    for( int i = 0; i < points.size(); ++i )
    {
        tmp = dist(points[i], points[next_index(i)]);
        if( tmp > max )
            max = tmp;
    }
    return max;
}

bool Polygon::if_all_angles_are_equal(const Polygon &other) const
{
    std::vector< double > this_angles;
    std::vector< double > other_angles;
    std::vector< double > reversed_other_angles;
    std::vector< Point > reversed_other = other.points;
    reverse_vector(reversed_other);
    for( int i = 0; i < points.size(); ++i )
    {
        this_angles.emplace_back( sin_A(points[i],points[previous_index(i)],points[next_index(i)]));
        other_angles.emplace_back(sin_A( other.points[i], other.points[previous_index(i)], other.points[next_index(i)]));
        reversed_other_angles.emplace_back(sin_A(reversed_other[i],reversed_other[previous_index(i)], reversed_other[next_index(i)]));
    }
    bool reversed, straight;
    for( int shift = 0; shift < points.size(); ++shift )
    {
        reversed = false;
        straight = false;
        for( int i = 0; i < this_angles.size(); ++i )
        {
            if( !equale_with_eps(this_angles[i], other_angles[next_shifted_index(i,shift,this_angles.size())]) )
                straight = true;
            if( !equale_with_eps(this_angles[i], reversed_other_angles[next_shifted_index(i,shift,this_angles.size())]))
                reversed = true;
        }
        if( !reversed || !straight )
            return true;
    }
    return false;
}

double Polygon::are_similar_with_coefficient(const Polygon &other, double similarity_coefficient) const
{
    if( !if_all_angles_are_equal(other))
        return false;
    std::vector< double > this_edges;
    std::vector< double > other_edges;
    std::vector< double > reversed_other_edges;
    std::vector< Point > other_points = other.points;
    reverse_vector( other_points );
    for( int i = 0; i < points.size(); ++i )
    {
        this_edges.emplace_back( dist( points[i], points[next_index(i)] ) );
        other_edges.emplace_back( dist( other.points[i], other.points[next_index(i)] ) );
        reversed_other_edges.emplace_back( dist( other_points[i], other_points[next_index(i)] ) );
    }
    bool straight, reversed;
    for( int shift = 0; shift < this_edges.size(); ++shift )
    {
        straight = false;
        reversed = false;
        for( int i = 0; i < this_edges.size(); ++i )
        {
            if( !equale_with_eps(  ( this_edges[i] / other_edges[next_shifted_index(i,shift,this_edges.size())]), similarity_coefficient ) )
                straight = true;
            if( ! equale_with_eps( (this_edges[i] / reversed_other_edges[next_shifted_index(i,shift,this_edges.size())]), similarity_coefficient))
                reversed = true;
        }
        if( !straight || !reversed )
            return true;
    }
    return false;
}

bool Polygon::isSimilarTo(const Shape &another) const
{
    try
    {
        const Polygon& other = dynamic_cast<const Polygon&>( another );
        if( other.verticesCount() == verticesCount() )
        {
            double this_max_edge = find_the_largest_edge();
            double other_max_edge = other.find_the_largest_edge();
            return are_similar_with_coefficient( other, this_max_edge / other_max_edge );
        }
        return false;
    }
    catch(...)
    {
        return false;
    }
}

//// Conguengerency of Polygon

bool Polygon::isCongruentTo(const Shape &another) const
{
    try
    {
        const Polygon& other = dynamic_cast<const Polygon&>(another);
        if( other.verticesCount() == verticesCount() )
            return are_similar_with_coefficient(other, 1.0);
        return false;
    }
    catch(...)
    {
        return false;
    }
}

//// == of Polygon

bool Polygon::is_there_any_shift_when_all_points_coincided(const std::vector<Point> &other_points) const
{
    bool any_dismatch = false;
    int count = points.size();
    for( int i = 0; i < count; ++ i )
    {
        any_dismatch = false;
        for( int j = 0; j < count; ++j )
        {
            Point tmp = other_points[next_shifted_index(j,i,count)];
            if( tmp != points[j])
            {
                any_dismatch = true;
                break;
            }
        }
        if( !any_dismatch )
            return true;
    }
    return false;
}

bool Polygon::operator==(const Shape &another) const
{
    try
    {
        const Polygon& other = dynamic_cast<const Polygon&>(another);
        if( other.verticesCount() == verticesCount() )
        {
            std::vector< Point > other_points = other.getVertices();
            bool first_wave = is_there_any_shift_when_all_points_coincided( other_points );
            reverse_vector(other_points);
            bool second_wave = is_there_any_shift_when_all_points_coincided( other_points );
            return first_wave || second_wave;
        }
        return false;
    }
    catch(...)
    {
        return false;
    }
}

//// != for Polygon, Rectangle, Square, Triangle, Ellipse and Circle

bool Polygon::operator!=(const Shape &another) const
{
    return !( *this == another );
}

bool Ellipse::operator!=(const Shape &another) const
{
    return !( *this == another );
}

bool Circle::operator!=(const Shape &another) const
{
    return !( *this == another );
}

//// Other functions

double dist( Point a, Point b )
{
    double dist = std::sqrt( ( a.x - b.x ) * ( a.x - b.x ) + ( a.y - b.y ) * ( a.y - b.y ) );
    return dist;
}

Point cut_centre( Point A, Point B )
{
    Point c;
    c.x = ( A.x + B.x ) / 2;
    c.y = ( A.y + B.y ) / 2;
    return c;
}

std::pair< double, double > square_equation( double a, double b, double c )
{
    std::pair< double, double > answer;
    double D = b * b - 4 * a * c;
    if( D < 0 )
    {
        answer.first = 0;
        answer.second = 0;
        return answer;
    }
    answer.first = ( -b + sqrt(D) ) / ( 2 * a );
    answer.second = ( -b - sqrt(D) ) / ( 2 * a );
    return answer;
}

Point find_normals_fundament( Point point, Line line )
{
    double a, b, c, x, y;
    line.get_normal_equation( a, b, c );
    if( equale_with_eps( a, 0 ) )
    {
        y = -c / b;
        x = ( ( y - point.y ) * a ) / b + point.x;
    }
    else if( equale_with_eps( b, 0 ) )
    {
        x = -c / a;
        y = ( ( x - point.x ) * b ) / a + point.y;
    }
    else
    {
        y = ( point.y * a - b *( ( c / a ) + point.x ) ) / ( a + ( b * b ) / a );
        x = - ( c / a) - ( ( b / a ) * y );
    }
    return Point( x, y );
}

bool equale_with_eps( double a, double b )
{
    double s = std::fabs(a - b);
    return  s <= Eps;
}

bool C_belongs_A_cut_B( Point C, Point A, Point B )
{
    double dist_A = dist( C, A);
    double dist_B = dist( C, B);
    double AB = dist( A, B);
    double sum = dist_A + dist_B;
    double a,b,c;
    Line line( A, B);
    line.get_normal_equation(a,b,c);
    double k = a * C.x + b * C.y + c;
    return equale_with_eps( sum, AB ) && equale_with_eps( k, 0);
}

bool line_intersection( Line first, Line second, Point& intersection_point )
{
    double a_1,b_1,c_1,a_2,b_2,c_2;
    first.get_normal_equation( a_1, b_1, c_1 );
    second.get_normal_equation( a_2, b_2, c_2);
    if( equale_with_eps(a_1 * b_2, a_2 * b_1 ) )
        return false;
    if( equale_with_eps( a_1, 0.0 ) )
    {
        intersection_point.y = -c_1 / b_1;
        intersection_point.x = ( -c_2 - b_2 * intersection_point.y ) / a_2;
    }
    else if( equale_with_eps( b_1, 0 ) )
    {
        intersection_point.x = -c_1 / a_1;
        intersection_point.y = ( -c_2 - a_2 * intersection_point.x ) / b_2;
    }
    else if( equale_with_eps( a_2, 0.0 ) )
    {
        intersection_point.y = -c_2 / b_2;
        intersection_point.x = ( -c_1 - b_1 * intersection_point.y ) / a_1;
    }
    else if( equale_with_eps( b_2, 0.0) )
    {
        intersection_point.x = -c_2 / a_2;
        intersection_point.y = ( -c_1 - a_1 * intersection_point.x ) / b_1;
    }
    else
    {
        intersection_point.y = ( ( ( c_1 * a_2 ) / a_1 ) - c_2 ) / ( b_2 - ( ( b_1 * a_2 ) / a_1 ) );
        intersection_point.x = ( -c_1 - b_1 * intersection_point.y ) / a_1;
    }
    return true;
}

std::pair< Point, Point > find_two_points_on_a_line_with_d_from_A( Point A, Line line, double d )
{
    double a, b, c;
    line.get_normal_equation(a,b,c);
    return specific_system_solver( a, b, c, A.x, A.y, d );
}

std::pair< Point, Point > specific_system_solver( double a, double b, double c, double x_0, double y_0, double d)
{
    Point first;
    Point second;
    if( equale_with_eps(a, 0) )
    {
        first.y = -c / b;
        second.y = -c / b;
        double s = sqrt((d * d) - (y_0 + ( c / b )) * (y_0 + ( c / b )));
        first.x = x_0 - s;
        second.x = x_0 + s;
    }
    else if( equale_with_eps(b, 0) )
    {
        first.x = - c / a;
        second.x = - c / a;
        double s = sqrt((d * d) - (x_0 - first.x) * (x_0 - first.x));
        first.y = y_0 - s;
        second.y = y_0 + s;
    }
    else
    {
        double aa, bb, cc;
        aa = 1 + ( b / a ) * ( b / a );
        bb = 2 *( ( b / a ) * ( x_0 + ( c / a ) ) - y_0 );
        cc = ( x_0 + ( c / a ) ) * ( x_0 + ( c / a ) ) + y_0 * y_0 - d * d;
        std::pair< double, double > y = square_equation( aa, bb, cc );
        first.y = y.first;
        second.y = y.second;
        first.x = ( -b / a ) * first.y - ( c / a );
        second.x = ( -b / a ) * second.y - ( c / a );
    }
    return std::make_pair( first, second );
}

Point recover_second_vertex_of_cut( Point A, Point centre )
{
    Point B;
    B.x = 2 * centre.x - A.x;
    B.y = 2 * centre.y - A.y;
    return B;
}

int next_shifted_index( int i, int shift, int whole_variants_count)
{
    return ( i + shift ) % whole_variants_count;
}

double sin_A( Point A, Point B, Point C )
{
    Triangle tri(A,B,C);
    double s = tri.area();
    double a = dist( A, B);
    double b = dist( A, C);
    double sin = (2.0 * s) / ( a * b );
    return sin;
}

#endif //GEOMETRY_GEOMETRY_CPP