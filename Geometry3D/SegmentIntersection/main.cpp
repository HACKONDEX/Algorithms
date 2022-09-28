#include <iostream>
#include <vector>
#include <cmath>

const double EPSILION = 0.000000001;

struct Point
{
  Point() : x(0), y(0), z(0) {}

  double x;
  double y;
  double z;

  Point& operator += (const Point& other)
  {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  Point& operator -= (const Point& other)
  {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  bool operator == (const Point& other) const
  {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct GeomVector
{
  GeomVector() : x(0), y(0), z(0) {}
  GeomVector(const Point& point) : x(point.x), y(point.y), z(point.z) {}
  GeomVector(const Point& a, const Point& b) : x(b.x -a.x), y(b.y - a.y), z(b.z - a.z) {}

  double x;
  double y;
  double z;

  GeomVector& operator += (const GeomVector& other)
  {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  GeomVector& operator -= (const GeomVector& other)
  {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  double scalar_square() const
  {
    return x * x + y * y + z * z;
  }

  bool is_parallel_to(const GeomVector& other) const
  {
    return x * other.y == y * other.x && x * other.z == z * other.x && y * other.z == z * other.y;
  }

  double get_length() const
  {
    return std::sqrt(x * x + y * y + z * z);
  }
};

struct Line
{
  Line(const Point& point, const GeomVector& vector) : l(point), v(vector){}

  Point l;
  GeomVector v;
};

GeomVector operator + (const GeomVector& first, const GeomVector& second)
{
  GeomVector sum = first;
  sum += second;
  return sum;
}

GeomVector operator - (const GeomVector& first, const GeomVector& second)
{
  GeomVector sum = first;
  sum -= second;
  return sum;
}

double dist(const Point& a, const Point& b)
{
  GeomVector A(a), B(b);
  GeomVector AB = B - A;
  return AB.get_length();
}

double scalar_multiplication(const GeomVector& a, const GeomVector& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

void get_two_variable_function(std::vector<double>& coefs, const GeomVector& a, const GeomVector& b, const GeomVector& c )
{
  coefs[0] = a.scalar_square();
  coefs[1] = b.scalar_square();
  coefs[2] = c.scalar_square();
  coefs[3] = 2 * scalar_multiplication(a, b);
  coefs[4] = 2 * scalar_multiplication(a, c);
  coefs[5] = 2 * scalar_multiplication(b, c);
}

double function_value_on_points(const std::vector<double>& coefs, double k, double t)
{
  return coefs[0] + coefs[1] * t * t + coefs[2] * k * k + coefs[3] * t + coefs[4] * k + coefs[5] * k * t;
}

bool find_extremum_points_of_the_function(const std::vector<double>& coefs, double& k, double& t)
{
  if (coefs[1] == 0)
  {
    if (coefs[5] == 0)
      return false;
    k = -coefs[3] / coefs[5];
    t = (-coefs[4] - 2 * coefs[2] * k) / coefs[5];
    return true;
  }
  if (coefs[2] == 0)
  {
    if(coefs[5] == 0)
      return false;
    t = -coefs[4] / coefs[5];
    k = (-coefs[3] - 2 * coefs[1] * t) / coefs[5];
    return true;
  }
  if(4 * coefs[2] * coefs[1] == coefs[5] * coefs[5])
    return false;
  k = (coefs[5] * coefs[3] - 2 * coefs[1] * coefs[4]) / (4 * coefs[2] * coefs[1] - coefs[5] * coefs[5]);
  t = (-coefs[3] - coefs[5] * k) / (2 * coefs[1]);
  return true;
}

double herons_method(const Point& a, const Point& b, const Point& c)
{
  if(a == c || b == c)
    return 0;
  GeomVector A(a), B(b), C(c);
  GeomVector AB = B - A, BC = C - B, CA = A - C;
  double half_perimetr = (AB.get_length() + BC.get_length() + CA.get_length()) / 2;
  double area = std::sqrt(half_perimetr * (half_perimetr - AB.get_length()) * (half_perimetr - BC.get_length()) * (half_perimetr - CA.get_length()));
  if(area == 0)
    return 0;
  return (2 * area) / AB.get_length();
}

bool are_collinear(const Point& a, const Point& b, const Point& c)
{
  GeomVector B(b), A(a);
  B -= A;
  return ((c.x - a.x) * B.y) == ((c.y - a.y) * B.x) && ((c.x - a.x) * B.z) == ((c.z - a.z) * B.x);
}

void get_the_base_of_the_perpendicular(const Point& point, const Line& line, Point& base)
{
  double px = point.x, py = point.y, pz = point.z, vx = line.v.x, vy = line.v.y, vz = line.v.z, lx = line.l.x,
      ly = line.l.y, lz = line.l.z;
  if(line.v.x == 0)/// Vx == 0
  {
    if(line.v.y == 0)
    {
      base.x = line.l.x;
      base.y = line.l.y;
      base.z = point.z;
    }
    else if(line.v.z == 0)
    {
      base.x = line.l.x;
      base.y = point.y;
      base.z = line.l.z;
    }
    else
    {
      base.x = line.l.x;
      base.z = (pz * vz * vz + py * vy * vz - ly * vy * vz + lz * vy * vy) / (vy * vy + vz * vz);
      base.y = (base.z * vy - lz * vy +ly * vz) / vz;
    }
  }
  else/// Vx != 0
  {
    if(vy == 0 && vz == 0)
    {
      base.x = px;
      base.y = ly;
      base.z = lz;
    }
    else if(vy == 0 && vz != 0)
    {
      base.y = ly;
      base.z = (lz * vx * vx - lx * vx * vz + pz * vz * vx + pz * vz *vz) / (vz * vz + vx * vx);
      base.x = (base.z * vx - lz * vx + lx * vz) / vz;
    }
    else if(vy != 0 && vz == 0)
    {
      base.z = lz;
      base.x = (px * vx * vx + lx * vy * vy - ly * vx * vy +py *vy * vx) / (vx * vx + vy * vy);
      base.y = (base.x * vy - lx * vy + ly * vx) / vx;
    }
    else
    {
      base.x = (px * vy * vx * vx + lx * vy * vy * vy - ly * vx * vy * vy + py * vy * vy * vx + lx * vy * vz * vz - ly * vx * vz * vz + ly * vz * vz * vx - lz * vx * vy * vz + pz * vz * vy * vx) / (vy * (vx * vx + vy * vy + vz * vz));
      base.y = (base.x * vy - lx * vy + ly * vx) / vx;
      base.z = (base.y * vz - ly * vz + lz * vy) / vy;
    }
  }
}

bool point_belongs_to_segment_ab(const Point& point, const Point& a, const Point& b)
{
  GeomVector A(a), B(b), P(point);
  GeomVector AB = B - A, AP = P - A, BP = P - B;
  return AB.get_length() - (AP.get_length() + BP.get_length()) < EPSILION;
}

double Dot(const GeomVector& a, const GeomVector& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

void VectorProduct(const GeomVector& first, const GeomVector& second, GeomVector& product)
{
  product.x = first.y * second.z - second.y * first.z;
  product.y = -first.x * second.z + second.x * first.z;
  product.z = first.x * second.y - second.x * first.y;
}

double OrientedAreaModulo(const GeomVector& a, const GeomVector& b) {

  GeomVector vector_product;
  VectorProduct(a, b, vector_product);
  return vector_product.get_length();
}

bool IsProjectionOfPointOnSegment(const Point& point, const Point& a, const Point& b) {
  GeomVector a_vector_point(a, point);
  GeomVector b_vector_point(b, point);
  GeomVector a_vector_b(a, b);
  GeomVector b_vector_a(b, a);
  return Dot(a_vector_b, a_vector_point) * Dot(b_vector_a, b_vector_point) >= 0;
}

double dist_from_point_to_segment_ab(const Point& x, const Point& a, const Point& b) {
  GeomVector a_vector_b(a, b);
//    GeomVector b_vector_a(b, a);
  GeomVector b_vector_x(b, x);
  GeomVector a_vector_x(a, x);
//    if(Dot(b_vector_a, b_vector_x) * Dot(a_vector_b, a_vector_x) >= 0) {
//        return  std::abs(OrientedAreaModulo(a_vector_x, a_vector_b)) / a_vector_b.Length();
//    }
  if(IsProjectionOfPointOnSegment(x, a, b)) {
    return std::abs(OrientedAreaModulo(b_vector_x, a_vector_x)) / a_vector_b.get_length();
  }
  return std::min(a_vector_x.get_length(), b_vector_x.get_length());
}

double if_segments_projections_doesnt_intersect(const Point& a, const Point& b, const Point& c, const Point& d)
{
  double a_cd = dist_from_point_to_segment_ab(a, c, d);
  double b_cd = dist_from_point_to_segment_ab(b, c, d);
  double c_ab = dist_from_point_to_segment_ab(c, a, b);
  double d_ab = dist_from_point_to_segment_ab(d, a, b);
  return std::min(std::min(a_cd, b_cd), std::min(c_ab, d_ab));
}

double find_dist_between_segments(const Point& a, const Point& b, const Point& c, const Point& d)
{
  GeomVector A(a), B(b), C(c), D(d);
  //// A + (B - A) * t == x, если точка х пренодлежит AB
  //// C + (D - C) * k == x, если точка х пренодлежит CD
  GeomVector A_C = A - C, B_A = B - A, D_C = D - C, C_D = C - D;
  GeomVector AB = B - A, AC = C - A, BC = C - B, AD = D - A, BD = D - B, DC = C - D;
  double ax = std::min(std::min(AC.get_length(), AD.get_length()), std::min(BC.get_length(), BD.get_length()));
  /// g(t, k) == C_0 + C_1 * t^2 + C_2 * k^2 + C_3 * t + C_4 * k + C_5 * kt
  if(B_A.is_parallel_to(D_C))/// Если параллельны, "Метог Герона"
  {

    if(are_collinear(a, b, c))
    {
      if(AB.get_length() == AC.get_length() + BC.get_length())
        return 0;
      if(AD.get_length() + BD.get_length() == AB.get_length())
        return 0;
      if(DC.get_length() == AC.get_length() + AD.get_length())
        return 0;
      if(DC.get_length() == BC.get_length() + BD.get_length())
        return 0;
      return std::min(std::min(AC.get_length(), AD.get_length()), std::min(BC.get_length(), BD.get_length()));
    }
    else
    {
      Point perpendicular_base_A, perpendicular_base_B;
      Line line(d, DC);
      get_the_base_of_the_perpendicular(a, line, perpendicular_base_A);
      get_the_base_of_the_perpendicular(b, line, perpendicular_base_B);
      GeomVector H_A(perpendicular_base_A), H_B(perpendicular_base_B);
      GeomVector H_A_C = C - H_A, H_A_D = D - H_A, H_B_C = C - H_B, H_B_D = D - H_B, H_A_H_B = H_B - H_A;
      if(DC.get_length() == H_A_C.get_length() + H_A_D.get_length())
      {
        return herons_method(a, b, c);
      }
      if(DC.get_length() == H_B_C.get_length() + H_B_D.get_length())
      {
        return herons_method(a, b, c);
      }
      if(H_A_H_B.get_length() == H_A_C.get_length() + H_B_C.get_length())
      {
        return herons_method(a, b, c);
      }
      if(H_A_H_B.get_length() == H_A_D.get_length() + H_B_D.get_length())
      {
        return herons_method(a, b, c);
      }
      return std::min(std::min(AC.get_length(), AD.get_length()), std::min(BC.get_length(), BD.get_length()));
    }
  }
  std::vector<double> function_coefs(6);
  get_two_variable_function(function_coefs, A_C, B_A, C_D);
  double k = 0, t = 0;
  double possible_answer = std::min(std::min(function_value_on_points(function_coefs, 0, 0),
                                             function_value_on_points(function_coefs, 1, 1)),
                                    std::min(function_value_on_points(function_coefs, 1, 0),
                                             function_value_on_points(function_coefs, 0, 1)));
  if(find_extremum_points_of_the_function(function_coefs, k, t)) {
//        std::cout << k << t << std::endl;
    if (k <= 1 && t <= 1 && k >= 0 && t >= 0)
      return std::min(sqrt(std::min(possible_answer, function_value_on_points(function_coefs, k, t))), ax);
  }
  return std::min(std::min(sqrt(possible_answer),ax), if_segments_projections_doesnt_intersect(a, b, c, d));
}

int main()
{
  std::cout.precision(10);
  Point a, b;//// Первый отрезок
  std::cin >> a.x  >> a.y >> a.z;
  std::cin >> b.x >> b.y >> b.z;
  Point c, d;/// Второй
  std::cin >> c.x  >> c.y >> c.z;
  std::cin >> d.x >> d.y >> d.z;
  std::cout << find_dist_between_segments(a, b, c, d) << std::endl;
  return 0;
}