#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <unordered_map>
#include <vector>

struct Coordinates;
struct Point;
struct GeomVector;
struct Face;
struct Edge;
struct ConvexHullBuilder;

struct Coordinates {
    Coordinates() : x(0), y(0), z(0) {}
    Coordinates(double x, double y, double z) : x(x), y(y), z(z) {}

    double x;
    double y;
    double z;
};

struct Point : Coordinates {
    Point() : Coordinates() {}
    Point(double x, double y, double z) : Coordinates(x, y, z) {}
    explicit Point(const GeomVector &vector);

    Point &operator=(const GeomVector &vector);
    bool operator==(const Point &other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct GeomVector : Coordinates {
    GeomVector() : Coordinates() {}
    GeomVector(double x, double y, double z) : Coordinates(x, y, z) {}
    explicit GeomVector(const Point &point) : Coordinates(point.x, point.y, point.z) {}
    GeomVector(const Point &a, const Point &b) : Coordinates(b.x - a.x, b.y - a.y, b.z - a.z) {}
    GeomVector(const GeomVector &other) : Coordinates(other.x, other.y, other.z) {}

    GeomVector &operator*=(double coefficient);
    GeomVector &operator+=(const GeomVector &other);
    GeomVector &operator-=(const GeomVector &other);

    bool Normalize();
};

Point::Point(const GeomVector &vector) {
    x = vector.x;
    y = vector.y;
    z = vector.z;
}

Point &Point::operator=(const GeomVector &vector) {
    x = vector.x;
    y = vector.y;
    z = vector.z;
    return *this;
}

GeomVector &GeomVector::operator*=(double coefficient) {
    x *= coefficient;
    y *= coefficient;
    z *= coefficient;
    return *this;
}

GeomVector &GeomVector::operator+=(const GeomVector &other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
}

GeomVector &GeomVector::operator-=(const GeomVector &other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
}

bool GeomVector::Normalize() {
    double length = std::sqrt(x * x + y * y + z * z);
    if (length == 0) {
        return false;
    }
    x /= length;
    y /= length;
    z /= length;
    return true;
}

/// Face - 3 points, and their indexes
struct Face {
    Face() = default;
    Face(const Point &a, int64_t index_a,
         const Point &b, int64_t index_b,
         const Point &c, int64_t index_c) : points(std::array<Point, 3>()), indexes(std::array<int64_t, 3>()) {
        points[0] = a;
        points[1] = b;
        points[2] = c;
        indexes[0] = index_a;
        indexes[1] = index_b;
        indexes[2] = index_c;
    }

    void GetProjectionOnFace(Point &projection, const Point &point64_t) const;
    void GetNormalVector(GeomVector &normal) const;
    double Dist(const Point &point);

    std::array<Point, 3> points;
    std::array<int64_t, 3> indexes;
};

struct Edge {
    Edge() = default;
    Edge(int64_t a_index, int64_t b_index, int64_t first_face_index, int64_t second_face_index)
            : point_index(std::array<int64_t, 2>()), face_index(std::array<int64_t, 2>()) {
        point_index[0] = a_index;
        point_index[1] = b_index;
        face_index[0] = first_face_index;
        face_index[1] = second_face_index;
    }

    std::array<int64_t, 2> point_index;
    std::array<int64_t, 2> face_index;
};

double Dot(const GeomVector &first, const GeomVector &second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

void CrossProduct(GeomVector &product, const GeomVector &first, const GeomVector &second) {
    product.x = first.y * second.z - first.z * second.y;
    product.y = first.z * second.x - first.x * second.z;
    product.z = first.x * second.y - first.y * second.x;
}

double OrientedVolume(const GeomVector &vector_1, const GeomVector &vector_2, const GeomVector &vector_3) {
    return vector_1.x * (vector_2.y * vector_3.z - vector_2.z * vector_3.y) -
           vector_1.y * (vector_2.x * vector_3.z - vector_2.z * vector_3.x) +
           vector_1.z * (vector_2.x * vector_3.y - vector_2.y * vector_3.x);
}

void Face::GetNormalVector(GeomVector &normal) const {
    CrossProduct(normal, GeomVector(points[0], points[1]), GeomVector(points[0], points[2]));
}

struct ConvexHullBuilder {
    explicit ConvexHullBuilder(const std::vector<Point> &points) : points(points),
                                                                   point64_t_to_faces_graph(std::move(
                                                                           std::vector<std::set<int64_t
                                                                           >>(points.size(), std::set<int64_t>()))),
                                                                   edge_point64_ts_matrix(std::move(
                                                                           std::vector<std::vector<
                                                                                   int64_t >>(points.size(),
                                                                                              std::vector<int64_t>(
                                                                                                      points.size(),
                                                                                                      -1)))),
                                                                   maximal_face_index(0), maximal_edge_index(0),
                                                                   current_point64_t_index(0) {
        horizon.reserve(points.size());
    }


    void FindAndSetInnerPoint64_t();
    int64_t BuildFace(Face &new_face, int64_t index_0, int64_t index_1, int64_t index_2);
    void BuildInitialTetrahedron();
    void BuildEdge(int64_t face_key_1, int64_t face_key_2, int64_t point64_t_index_1,
                   int64_t point64_t_index_2);
    bool IsFaceVisible(const Face &face, int64_t point64_t_index);
    bool IsFaceVisible(int64_t face_key, int64_t point64_t_index);
    void BuildInitialGraph();
    Edge &GetEdge(int64_t point64_t_index_1, int64_t point64_t_index_2);
    int64_t GetEdgeKey(int64_t point64_t_index_1, int64_t point64_t_index_2);
    bool AddOutlineEdgeToHorizon(int64_t point64_t_index_1, int64_t point64_t_index_2);
    void GetHorizon(const std::set<int64_t> &visible_faces);
    void UpdateEdgeList(Edge &current_edge);
    void AddNewFace(Edge &near_edge);
    void UpdateGraph(const std::set<int64_t> &first_point64_ts_set, const std::set<int64_t> &second_point64_ts_set);
    void AddNewFacesFromHorizon();
    void RemoveOldFaces();
    void BuildConvexHull();

    void PrintGraph();
    void PrintHorizon();
    void PrintEdgeList();
    void PrintConvexHull();
    void PrintVisibleFaces(const std::set<int64_t> &visible_faces);
    void PrintFace(int64_t i);

    const std::vector<Point> &points;
    Point inner_point64_t;
    std::unordered_map<int64_t, Face> face_list;
    std::vector<std::set<int64_t>> point64_t_to_faces_graph;
    std::unordered_map<int64_t, std::set<int64_t>> face_to_point64_ts_graph;
    std::unordered_map<int64_t, Edge> edge_list;
    std::vector<std::vector<int64_t>> edge_point64_ts_matrix;
    std::set<int64_t> convex_hull; /// Hull face indices :)
    int64_t maximal_face_index;
    int64_t maximal_edge_index;
    int64_t current_point64_t_index;
    std::vector<int64_t> horizon;
    std::set<int64_t> faces_to_delete;
};

void ConvexHullBuilder::PrintFace(int64_t i) {
    Face &y = face_list[i];
    std::cout << y.indexes[0] << " " << y.indexes[1] << " " << y.indexes[2] << std::endl;
}

void ConvexHullBuilder::PrintVisibleFaces(const std::set<int64_t> &visible_faces) {
    std::cout << "visible faces from point64_t \n";
    for (int64_t x: visible_faces) {
        Face &y = face_list[x];
        std::cout << y.indexes[0] << " " << y.indexes[1] << " " << y.indexes[2] << std::endl;
    }
}

void ConvexHullBuilder::PrintGraph() {
    std::cout << "Point64_ts to faces \n";
    for (int64_t i = current_point64_t_index; i < point64_t_to_faces_graph.size(); ++i) {
        std::cout << i << " ---> ";
        for (const auto &x: point64_t_to_faces_graph[i]) {
            PrintFace(x);
        }
        std::cout << '\n';
    }

    std::cout << "Faces to point64_ts\n";
    for (const auto &x: face_to_point64_ts_graph) {
        Face &z = face_list[x.first];
        std::cout << z.indexes[0] << " " << z.indexes[1] << " " << z.indexes[2] << " ---> ";
        for (int64_t y: x.second) {
            std::cout << y << " ";
        }
        std::cout << '\n';
    }
}

void ConvexHullBuilder::PrintEdgeList() {
    std::cout << "Edge List \n";
    for (auto &x: edge_list) {
        std::cout << x.first << " ---> " << x.second.point_index[0] << " " << x.second.point_index[1] << '\n';
    }
}

void ConvexHullBuilder::PrintHorizon() {
    std::cout << "Horizon" << '\n';
    std::cout << horizon.size() << '\n';
    for (int64_t x: horizon) {
        Edge &edge_ = edge_list[x];
        std::cout << x << " ---> " << edge_.point_index[0] << " " << edge_.point_index[1] << '\n';
    }
}

void ConvexHullBuilder::PrintConvexHull() {
    std::cout << "Convex hull \n";
    std::cout << convex_hull.size() << std::endl;
    for (int64_t x: convex_hull) {
        Face &y = face_list[x];
        std::cout << y.indexes[0] << " " << y.indexes[1] << " " << y.indexes[2] << '\n';
    }
}

void ConvexHullBuilder::FindAndSetInnerPoint64_t() {
    GeomVector vector_a(points[0]);
    vector_a += GeomVector(points[1]);
    vector_a *= 0.5;
    vector_a -= GeomVector(points[2]);
    vector_a *= 2.0 / 3.0;
    vector_a += GeomVector(points[2]);
    vector_a += GeomVector(points[3]);
    vector_a *= 0.5;
    inner_point64_t = vector_a;
}

int64_t ConvexHullBuilder::BuildFace(Face &new_face, int64_t index_0, int64_t index_1, int64_t index_2) {
    new_face.points[0] = points[index_0];
    new_face.points[1] = points[index_1];
    new_face.points[2] = points[index_2];
    new_face.indexes[0] = index_0;
    new_face.indexes[1] = index_1;
    new_face.indexes[2] = index_2;
    face_list.insert({maximal_face_index, new_face});
    convex_hull.insert(maximal_face_index);
    ++maximal_face_index;
    return maximal_face_index - 1;
}

void ConvexHullBuilder::BuildEdge(int64_t face_key_1, int64_t face_key_2, int64_t point64_t_index_1,
                                  int64_t point64_t_index_2) {
    Edge new_edge{};
    new_edge.face_index[0] = face_key_1;
    new_edge.face_index[1] = face_key_2;
    new_edge.point_index[0] = point64_t_index_1;
    new_edge.point_index[1] = point64_t_index_2;
    edge_list.insert({maximal_edge_index, new_edge});
    edge_point64_ts_matrix[point64_t_index_1][point64_t_index_2] = edge_point64_ts_matrix[point64_t_index_2][point64_t_index_1] = maximal_edge_index;
    ++maximal_edge_index;
}

void ConvexHullBuilder::BuildInitialTetrahedron() {
    Face face_0_1_2;
    Face face_1_0_3;
    Face face_2_3_0;
    Face face_3_2_1;

    int64_t key_0_1_2 = BuildFace(face_0_1_2, 0, 1, 2);
    int64_t key_1_0_3 = BuildFace(face_1_0_3, 0, 1, 3);
    int64_t key_2_3_0 = BuildFace(face_2_3_0, 0, 2, 3);
    int64_t key_3_2_1 = BuildFace(face_3_2_1, 1, 2, 3);

    BuildEdge(key_0_1_2, key_1_0_3, 0, 1);
    BuildEdge(key_0_1_2, key_2_3_0, 0, 2);
    BuildEdge(key_0_1_2, key_3_2_1, 1, 2);
    BuildEdge(key_1_0_3, key_3_2_1, 1, 3);
    BuildEdge(key_3_2_1, key_2_3_0, 2, 3);
    BuildEdge(key_1_0_3, key_2_3_0, 0, 3);

    current_point64_t_index = 4;
}

bool ConvexHullBuilder::IsFaceVisible(const Face &face, int64_t point64_t_index) {
    GeomVector vector_0_point64_t(face.points[0], points[point64_t_index]);
    GeomVector vector_0_1(face.points[0], face.points[1]);
    GeomVector vector_0_2(face.points[0], face.points[2]);
    GeomVector vector_0_inner_point64_t(face.points[0], inner_point64_t);
    return OrientedVolume(vector_0_1, vector_0_2, vector_0_point64_t) *
           OrientedVolume(vector_0_1, vector_0_2, vector_0_inner_point64_t) < 0.0;
}

bool ConvexHullBuilder::IsFaceVisible(int64_t face_key, int64_t point64_t_index) {
    const Face &face_ = face_list[face_key];
    return IsFaceVisible(face_, point64_t_index);
}

void ConvexHullBuilder::BuildInitialGraph() {
    auto point64_ts_count = static_cast<int64_t>(points.size());
    for (int64_t i = 0; i < 4; ++i) {
        face_to_point64_ts_graph.insert({i, std::set<int64_t>()});
    }
    for (int64_t i = 4; i < point64_ts_count; ++i) {
        for (const auto &face_element: face_list) {
            if (IsFaceVisible(face_element.second, i)) {
                point64_t_to_faces_graph[i].insert(face_element.first);
                face_to_point64_ts_graph[face_element.first].insert(i);
            }
        }
    }
}

Edge &ConvexHullBuilder::GetEdge(int64_t point64_t_index_1, int64_t point64_t_index_2) {
    return edge_list[edge_point64_ts_matrix[point64_t_index_1][point64_t_index_2]];
}

int64_t ConvexHullBuilder::GetEdgeKey(int64_t point64_t_index_1, int64_t point64_t_index_2) {
    return edge_point64_ts_matrix[point64_t_index_1][point64_t_index_2];
}

bool ConvexHullBuilder::AddOutlineEdgeToHorizon(int64_t point64_t_index_1, int64_t point64_t_index_2) {
    int64_t edge_index = edge_point64_ts_matrix[point64_t_index_1][point64_t_index_2];
    if (edge_index == -1) {
        return false;
    }
    Edge &edge_ = edge_list[edge_index];
    bool first = IsFaceVisible(edge_.face_index[0], current_point64_t_index);
    bool second = IsFaceVisible(edge_.face_index[1], current_point64_t_index);
    if (first ^ second) {
        horizon.push_back(edge_index);
        if (second) {
            std::swap(edge_.face_index[0], edge_.face_index[1]);
        }
        return true;
    }
    edge_list.erase(edge_index);
    edge_point64_ts_matrix[point64_t_index_1][point64_t_index_2] = edge_point64_ts_matrix[point64_t_index_2][point64_t_index_1] = -1;
    return false;
}

void ConvexHullBuilder::GetHorizon(const std::set<int64_t> &visible_faces) {
    for (int64_t face_key: visible_faces) {
        Face &current_face = face_list[face_key];
        bool is_connected_to_horizon_1 = AddOutlineEdgeToHorizon(current_face.indexes[0], current_face.indexes[1]);
        bool is_connected_to_horizon_2 = AddOutlineEdgeToHorizon(current_face.indexes[1], current_face.indexes[2]);
        bool is_connected_to_horizon_3 = AddOutlineEdgeToHorizon(current_face.indexes[2], current_face.indexes[0]);
        if (!(is_connected_to_horizon_1 || is_connected_to_horizon_2 || is_connected_to_horizon_3)) {
            faces_to_delete.insert(face_key);
        }
    }
}

void ConvexHullBuilder::UpdateEdgeList(Edge &current_edge) {
    for (int64_t index = 0; index <= 1; ++index) {
        if (edge_point64_ts_matrix[current_edge.point_index[index]][current_point64_t_index] == -1) {
            edge_point64_ts_matrix[current_edge.point_index[index]][current_point64_t_index] =
            edge_point64_ts_matrix[current_point64_t_index][current_edge.point_index[index]] = maximal_edge_index;
            Edge new_edge(current_point64_t_index, current_edge.point_index[index], maximal_face_index, -1);
            edge_list.insert({maximal_edge_index, new_edge});
            ++maximal_edge_index;
        } else {
            Edge &edge_ = GetEdge(current_point64_t_index, current_edge.point_index[index]);
            if (edge_.face_index[0] == -1) {
                edge_.face_index[0] = maximal_face_index;
            } else {
                edge_.face_index[1] = maximal_face_index;
            }
        }
    }
}

void ConvexHullBuilder::UpdateGraph(const std::set<int64_t> &first_point64_ts_set,
                                    const std::set<int64_t> &second_point64_ts_set) {
    std::set<int64_t> &new_face_to_point64_ts = face_to_point64_ts_graph[maximal_face_index];
    for (int64_t point64_t_index: first_point64_ts_set) {
        if (point64_t_index != current_point64_t_index && IsFaceVisible(maximal_face_index, point64_t_index)) {
            point64_t_to_faces_graph[point64_t_index].insert(maximal_face_index);
            new_face_to_point64_ts.insert(point64_t_index);
        }
    }
    for (int64_t point64_t_index: second_point64_ts_set) {
        if (point64_t_index != current_point64_t_index && IsFaceVisible(maximal_face_index, point64_t_index)) {
            point64_t_to_faces_graph[point64_t_index].insert(maximal_face_index);
            new_face_to_point64_ts.insert(point64_t_index);
        }
    }
}

void ConvexHullBuilder::AddNewFace(Edge &near_edge) {
    Face new_face(points[near_edge.point_index[0]], near_edge.point_index[0],
                  points[near_edge.point_index[1]], near_edge.point_index[1],
                  points[current_point64_t_index], current_point64_t_index);
    face_list.insert({maximal_face_index, new_face});
    convex_hull.insert(maximal_face_index);
    face_to_point64_ts_graph.insert({maximal_face_index, std::set<int64_t>()});

}

void ConvexHullBuilder::AddNewFacesFromHorizon() {
    for (int64_t current_edge_key: horizon) {
        Edge &current_edge = edge_list[current_edge_key];
        std::set<int64_t> &visible_point64_ts_from_deletable_face = face_to_point64_ts_graph[current_edge.face_index[0]];
        std::set<int64_t> &visible_point64_ts = face_to_point64_ts_graph[current_edge.face_index[1]];
        visible_point64_ts_from_deletable_face.erase(current_point64_t_index);
        faces_to_delete.insert(current_edge.face_index[0]);
        AddNewFace(current_edge);
        UpdateEdgeList(current_edge);
        UpdateGraph(visible_point64_ts_from_deletable_face, visible_point64_ts);
        current_edge.face_index[0] = maximal_face_index;
        ++maximal_face_index;
    }
}

void ConvexHullBuilder::RemoveOldFaces() {
    for (int64_t face_key: faces_to_delete) {
        std::set<int64_t> &visible_point64_ts = face_to_point64_ts_graph[face_key];
        for (int64_t point64_t_index: visible_point64_ts) {
            point64_t_to_faces_graph[point64_t_index].erase(face_key);
        }
        face_to_point64_ts_graph.erase(face_key);
        convex_hull.erase(face_key);
        face_list.erase(face_key);
    }
}

void ConvexHullBuilder::BuildConvexHull() {
    FindAndSetInnerPoint64_t();
    BuildInitialTetrahedron();
    BuildInitialGraph();
    auto point64_ts_count = static_cast<int64_t>(points.size());
    for (current_point64_t_index = 4; current_point64_t_index < point64_ts_count; ++current_point64_t_index) {
        std::set<int64_t> &visible_faces = point64_t_to_faces_graph[current_point64_t_index];
        GetHorizon(visible_faces);
        visible_faces.clear();
        AddNewFacesFromHorizon();
        RemoveOldFaces();
        horizon.resize(0);
        faces_to_delete.clear();
    }
}

void Face::GetProjectionOnFace(Point &projection, const Point &point64_t) const {
    GeomVector normal;
    GetNormalVector(normal);
    GeomVector normalized_normal(normal);
    normalized_normal.Normalize();
    double plane_coefficient_d = -points[0].x * normal.x + points[0].y * normal.y + points[0].z * normal.z;
    double dist_to_projection =
            -(normal.z * point64_t.x + normal.y * point64_t.y + normal.x * point64_t.z + plane_coefficient_d)
            / (normal.x * normalized_normal.x + normal.y * normalized_normal.y + normal.z * normalized_normal.z);
    normalized_normal *= dist_to_projection;
    projection.x = point64_t.x + normalized_normal.x;
    projection.y = point64_t.y + normalized_normal.y;
    projection.z = point64_t.z + normalized_normal.z;
}

double DistOfPoint64_ts(const Point &first, const Point &second) {
    return std::sqrt((second.x - first.x) * (second.x - first.x) + (second.y - first.y) * (second.y - first.y) +
                     (second.z - first.z) * (second.z - first.z));
}

double DistFromPoint64_tToFace(const Face &face_, const Point &point64_t_) {
    Point projection;
    face_.GetProjectionOnFace(projection, point64_t_);
    return DistOfPoint64_ts(projection, point64_t_);
}

double Face::Dist(const Point &point) {
    GeomVector a(points[0], points[1]);
    GeomVector b(points[0], points[2]);
    GeomVector normal;
    CrossProduct(normal, a, b);
    double d = -(normal.x * points[0].x + normal.y * points[0].y + normal.z * points[0].z);
    double numerator = d + normal.x * point.x + normal.y * point.y + normal.z * point.z;
    double denominator = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    return std::abs(numerator) / denominator;
}

int main() {
    int64_t point64_ts_count = 0;
    std::cin >> point64_ts_count;
    std::vector<Point> point64_ts(point64_ts_count);
    for (int64_t i = 0; i < point64_ts_count; ++i) {
        std::cin >> point64_ts[i].x >> point64_ts[i].y >> point64_ts[i].z;
    }
    std::random_device random_device;
    std::mt19937 g(random_device());
    std::shuffle(point64_ts.begin(), point64_ts.end(), g);
    ConvexHullBuilder builder(point64_ts);
    builder.BuildConvexHull();
    std::vector<Face> faces;
    faces.reserve(builder.convex_hull.size());
    for (int64_t face_key: builder.convex_hull) {
        faces.push_back(builder.face_list[face_key]);
    }
}