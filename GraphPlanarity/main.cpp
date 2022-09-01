#include <algorithm>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_set>
#include <vector>

class Graph {
public:
    explicit Graph(int vert_count = 0) : vertex_count(vert_count), edge_count(0) {
        for (int i = 0; i < vert_count; ++i) {
            graph.emplace_back();
        }

    }

    ~Graph() = default;

    Graph(const Graph &) = delete;

    Graph(Graph &&) = delete;

    Graph &operator=(const Graph &) = delete;

    Graph &operator=(Graph &&) = delete;

    void AddVertex();

    void AddEdge(int from, int to);

    void RemoveEdge(int from, int to);

    bool HasEdge(int from, int to) const;

    const std::vector<int> &GetNextVertices(int from) const;

    int GetVertexCount();

    int GetEdgeCount();

private:
    int vertex_count;
    int edge_count;
    std::vector<std::vector<int>> graph;

};

void Graph::AddVertex() {
    ++vertex_count;
    graph.emplace_back();
}

int Graph::GetEdgeCount() {
    return edge_count;
}

void Graph::AddEdge(int from, int to) {
    graph[from].push_back(to);
    ++edge_count;
}

void Graph::RemoveEdge(int from, int to) {
    std::vector<int> perm;
    for (int i = 0; i < graph[from].size(); ++i) {
        if (graph[from][i] != to) {
            perm.push_back(graph[from][i]);
        }
    }
    --edge_count;
    graph[from] = perm;
}

bool Graph::HasEdge(int from, int to) const {
    for (int i = 0; i < graph[from].size(); ++i) {
        if (graph[from][i] == to) {
            return true;
        }
    }
    return false;
}

const std::vector<int> &Graph::GetNextVertices(int from) const {
    return graph[from];
}

int Graph::GetVertexCount() {
    return vertex_count;
}

//// Solution

void BuildGraph(Graph *graph, int edge, std::vector<int>& from, std::vector<int>& to, int **matrix_graph) {
    for (int i = 0; i < edge; ++i) {
        if (from[i] != to[i]) {
            graph->AddEdge(from[i], to[i]);
            graph->AddEdge(to[i], from[i]);
            matrix_graph[from[i]][to[i]] = 1;
            matrix_graph[to[i]][from[i]] = 1;
        }
    }
}

struct Node {
    int number;
    int parent;

    Node(int v = -1, int p = -1) : number(v), parent(p) {}
};

void IfNotColoured(std::stack<Node> &stack, std::vector<int> &next_vertices, Node &current_node,
                   unsigned char *colour, int *in_time, int *low) {
    int to;
    for (int i = 0; i < next_vertices.size(); ++i) {
        to = next_vertices[i];
        if (to == current_node.parent) {
            continue;
        } else if (colour[to] == 0) {
            stack.push(Node(to, current_node.number));
        } else {
            low[current_node.number] = std::min(low[current_node.number], in_time[to]);
        }
    }
}

void IfColoured(std::vector<int> &next_vertices, Node &current_node, int *in_time, int *low, int **matrix_graph) {
    int to;
    for (int i = 0; i < next_vertices.size(); ++i) {
        to = next_vertices[i];
        if (current_node.parent != to && low[to] < low[current_node.number])
            low[current_node.number] = low[to];
    }
    for (int i = 0; i < next_vertices.size(); ++i) {
        to = next_vertices[i];
        if (current_node.parent != to && low[to] > in_time[current_node.number])
            if (matrix_graph[to][current_node.number] == 1) {
                matrix_graph[to][current_node.number] = -1;
                matrix_graph[current_node.number][to] = -1;
            }
    }
}

void FindBridges(Graph *graph, int **matrix_graph) {
    int vertices = graph->GetVertexCount();
    auto *colour = new unsigned char[vertices];
    int *low = new int[vertices];
    auto *in_time = new int[vertices];
    for (int i = 0; i < vertices; ++i) {
        colour[i] = 0;
    }

    std::stack<Node> stack;
    std::vector<int> next_vertices;
    int timer = 0;
    Node current_node(0, 0);

    for (int i = 0; i < vertices; ++i) {
        if (colour[i] == 0) {
            stack.push(Node(i, -1));

            while (!stack.empty()) {
                current_node = stack.top();
                next_vertices = graph->GetNextVertices(current_node.number);
                if (colour[current_node.number] == 0) {
                    colour[current_node.number] = 1;
                    low[current_node.number] = in_time[current_node.number] = timer++;
                    IfNotColoured(stack, next_vertices, current_node, colour, in_time, low);
                } else {
                    stack.pop();
                    colour[current_node.number] = 2;
                    IfColoured(next_vertices, current_node, in_time, low, matrix_graph);
                }
            }
        }
    }

    delete[] in_time;
    delete[] low;
    delete[] colour;
}

void RemoveBridgesFromGraph(Graph *graph, int **matrix_graph) {
    int vertices = graph->GetVertexCount();
    FindBridges(graph, matrix_graph);

    for (int i = 0; i < vertices; ++i) {
        for (int j = i; j < vertices; ++j) {
            if (matrix_graph[i][j] == -1) {
                matrix_graph[i][j] = 0;
                matrix_graph[j][i] = 0;
                graph->RemoveEdge(i, j);
                graph->RemoveEdge(j, i);
            }
        }
    }
}

void BFSMarkingConnectionComponentToWhichVertexBelongs(int **matrix_graph, int start, std::vector<bool>& checked,
                                                       int vertices, std::vector<int> &component) {
    std::queue<int> q;
    checked[start] = true;
    q.push(start);
    int vertex = 0;
    while (!q.empty()) {
        vertex = q.front();
        component.push_back(vertex);
        q.pop();
        for (int i = 0; i < vertices; ++i) {
            if (matrix_graph[vertex][i] == 1 && !checked[i]) {
                checked[i] = true;
                q.push(i);
            }
        }
    }
}

void RecoverCycle(std::vector<int> &cycle, std::stack<Node> &stack, int cycle_end, std::vector<unsigned char> &colour) {
    cycle.push_back(cycle_end);
    Node v;

    while (!stack.empty()) {
        v = stack.top();

        stack.pop();
        if (v.number == cycle_end) {
            break;
        } else {
            if (colour[v.number] == 1) {
                cycle.push_back(v.number);
                colour[v.number] = 2;
            }
        }
    }
}

bool DFSAnyCycle(Graph *graph, std::stack<Node> &stack, int &cycle_end, std::vector<unsigned char> &colour) {
    bool ans = false;
    Node v;
    int to;
    bool found_cycle = false;
    std::vector<int> next_vertices;
    while (!stack.empty()) {
        v = stack.top();
        if (colour[v.number] == 0) {
            colour[v.number] = 1;
            next_vertices = graph->GetNextVertices(v.number);
            for (int i = 0; i < next_vertices.size(); ++i) {
                to = next_vertices[i];
                if (to == v.parent)
                    continue;
                else if (colour[to] == 0)
                    stack.push(Node(to, v.number));
                else if (colour[to] == 1) {
                    cycle_end = to;
                    ans = true;
                    found_cycle = true;
                    break;
                }
            }

            if (found_cycle) {
                break;
            }
        } else {
            stack.pop();
            colour[v.number] = 2;
        }
    }
    return ans;
}

//returns false, if there is no cycle
bool FindCycle(Graph *graph, int **matrix_graph, int start, std::vector<int> &cycle) {
    bool ans = false;
    int vertices = graph->GetVertexCount();

    std::vector<unsigned char> colour(vertices);
    for (int i = 0; i < vertices; ++i)
        colour[i] = 0;

    std::stack<Node> stack;

    //first vertex has no parent
    stack.push(Node(start, -1));
    Node v;
    int to = 0;

    //how to remember cycle
    int cycle_end = -1;
    ans = DFSAnyCycle(graph, stack, cycle_end, colour);

    if (ans) {
        RecoverCycle(cycle, stack, cycle_end, colour);
    }
    return ans;
}

struct Area {
    std::vector<int> boarder;
    std::unordered_set<int> vertices;
};

struct Segment {
    std::vector<int> vertices;
    std::unordered_set<int> hash_vertices;
    std::vector<int> contact_vertices;
    std::vector<int> involved_areas;
};

void PlaceCycleOnSurface(int **matrix_graph, std::vector<bool>& placed, std::vector<Area> &areas,
                         std::vector<int> &cycle, std::vector<int> &placed_vertices) {
    int j = 0;
    Area area_1;
    Area area_2;
    while (j < cycle.size() - 1) {
        placed[cycle[j]] = true;
        placed_vertices.push_back(cycle[j]);

        area_1.boarder.push_back(cycle[j]);

        matrix_graph[cycle[j]][cycle[j + 1]] = 2;
        matrix_graph[cycle[j + 1]][cycle[j]] = 2;

        ++j;
    }
    placed[cycle[j]] = true;
    placed_vertices.push_back(cycle[j]);

    area_1.boarder.push_back(cycle[j]);
    area_2.boarder = area_1.boarder;
    std::reverse(area_2.boarder.begin(), area_2.boarder.end());
    areas.push_back(area_1);
    areas.push_back(area_2);

    for (int j = 0; j < areas.size(); ++j) {
        for (int i = 0; i < areas[j].boarder.size(); ++i) {
            areas[j].vertices.insert(areas[j].boarder[i]);
        }
    }
    matrix_graph[cycle[j]][cycle[0]] = 2;
    matrix_graph[cycle[0]][cycle[j]] = 2;
}

void FindAreasForSegment(std::vector<Area> &areas, Segment &seg) {
    int v;
    bool belongs_to_area = true;
    for (int i = 0; i < areas.size(); ++i) {
        belongs_to_area = true;
        for (int j = 0; j < seg.contact_vertices.size(); ++j) {
            v = seg.contact_vertices[j];
            if (areas[i].vertices.count(v) == 0) {
                belongs_to_area = false;
                break;
            }
        }
        if (belongs_to_area) {
            seg.involved_areas.push_back(i);
        }
    }
}

void IsSegmentBest(bool &any_segment, Segment &best_seg, Segment &seg) {
    if (!any_segment) {
        any_segment = true;
        best_seg = seg;
    } else if (seg.involved_areas.size() < best_seg.involved_areas.size())
        best_seg = seg;
}

void FindComplexSegment(Graph *graph, std::vector<int> &comp_vertices, std::queue<int> &q, Segment &seg,
                        int &current_component, int j, std::vector<bool>& segmentated, std::vector<bool>& placed) {
    segmentated[current_component] = true;
    seg.vertices.clear();
    seg.involved_areas.clear();
    seg.contact_vertices.clear();
    seg.vertices.push_back(current_component);
    seg.contact_vertices.push_back(current_component);

    q.push(comp_vertices[j]);
    int v, to = 0;
    std::vector<int> next_vertices;
    while (!q.empty()) {
        v = q.front();
        q.pop();
        if (placed[v] && v != current_component && !segmentated[v]) {
            segmentated[v] = true;
            seg.vertices.push_back(v);
            seg.contact_vertices.push_back(v);
        } else if (!placed[v] && !segmentated[v]) {
            segmentated[v] = true;
            seg.vertices.push_back(v);
            next_vertices = graph->GetNextVertices(v);

            for (int i = 0; i < next_vertices.size(); ++i) {
                to = next_vertices[i];
                if (!segmentated[to] && to != current_component) {
                    q.push(to);
                }
            }
        }
    }
}

void ComplexSegment(Graph *graph, std::vector<int> &comp_vertices, std::vector<int> &placed_vertices,
                    std::vector<Area> &areas, std::vector<bool>& segmentated, std::vector<bool>& placed, Segment &best_seg, bool &any_segment) {
    int current_component;
    for (int i = 0; i < placed_vertices.size(); ++i) {
        current_component = placed_vertices[i];

        Segment seg;
        std::queue<int> q;
        comp_vertices = graph->GetNextVertices(current_component);
        for (int j = 0; j < comp_vertices.size(); ++j) {
            if (!segmentated[comp_vertices[j]] && !placed[comp_vertices[j]]) {
                FindComplexSegment(graph, comp_vertices, q, seg, current_component, j, segmentated, placed);

                FindAreasForSegment(areas, seg);

                for (int i = 0; i < placed_vertices.size(); ++i) {
                    segmentated[placed_vertices[i]] = false;
                }

                IsSegmentBest(any_segment, best_seg, seg);

                if (best_seg.contact_vertices.size() == 0) {
                    return;
                }
            }
        }
    }
}

void
PrimitiveSegment(Graph *graph, int **matrix_graph, std::vector<int> &comp_vertices, std::vector<int> &placed_vertices,
                 std::vector<Area> &areas, std::vector<bool>& segmentated, std::vector<bool>& placed, Segment &best_seg, bool &any_segment) {
    Segment seg;
    int current_component, to;
    std::vector<int> next_vertices;
    for (int i = 0; i < placed_vertices.size(); ++i) {
        current_component = placed_vertices[i];
        segmentated[current_component] = true;
        next_vertices = graph->GetNextVertices(current_component);
        for (int i = 0; i < next_vertices.size(); ++i) {
            to = next_vertices[i];

            if (!segmentated[to] && placed[to] && matrix_graph[to][current_component] != 2) {
                seg.vertices.clear();
                seg.contact_vertices.clear();
                seg.involved_areas.clear();
                seg.contact_vertices.push_back(to);
                seg.contact_vertices.push_back(current_component);
                seg.vertices = seg.contact_vertices;

                FindAreasForSegment(areas, seg);
                IsSegmentBest(any_segment, best_seg, seg);

                if (best_seg.contact_vertices.empty()) {
                    return;
                }
            }
        }
    }
}

Segment FindBestSegment(Graph *graph, int **matrix_graph, std::vector<bool>& placed, std::vector<Area> &areas,
                        std::vector<int> &placed_vertices, std::vector<bool>& segmentated, bool &any_segment) {
    Segment best_seg;
    int verts_count = graph->GetVertexCount();
    // check if the vert was taken as a segment Node
    for (int i = 0; i < verts_count; ++i) {
        segmentated[i] = false;
    }

    std::vector<int> comp_vertices;

    // cycle from contact vertex// complex one
    ComplexSegment(graph, comp_vertices, placed_vertices,
                   areas, segmentated, placed, best_seg, any_segment);

    if (best_seg.contact_vertices.size() == 0) {
        return best_seg;
    }

    PrimitiveSegment(graph, matrix_graph, comp_vertices, placed_vertices,
                     areas, segmentated, placed, best_seg, any_segment);
    return best_seg;
}

void RecoverChain(std::vector<int> &chain, std::stack<Node> &stack, std::vector<unsigned char>& colour, int finish) {
    Node v;
    chain.push_back(finish);
    while (!stack.empty()) {
        v = stack.top();
        stack.pop();
        if (colour[v.number] == 1) {
            colour[v.number] = 2;
            chain.push_back(v.number);
        }
    }
}

void DFSForChain(Graph *graph, std::vector<int> &chain, std::stack<Node> &stack, std::vector<unsigned char>& colour,
                 std::vector<bool>& placed, Segment &seg, int &finish) {
    std::vector<int> next_vertices;
    Node v;
    int to;
    bool break_from_dfs = false;
    while (!stack.empty()) {
        v = stack.top();
        if (colour[v.number] == 0) {
            colour[v.number] = 1;
            next_vertices = graph->GetNextVertices(v.number);
            for (int i = 0; i < next_vertices.size(); ++i) {
                to = next_vertices[i];
                if (seg.hash_vertices.count(to) > 0) {
                    if (v.parent == to) {
                        continue;
                    }
                    if (!placed[to] && colour[to] == 0) {
                        stack.push(Node(to, v.number));
                    }
                    else if (placed[to] && v.parent != -1) {
                        finish = to;
                        break_from_dfs = true;
                        break;
                    }
                }
            }
        } else if (colour[v.number] == 1) {
            stack.pop();
            colour[v.number] = 2;
        }
        if (break_from_dfs) {
            break;
        }
    }
}

void GetChainFromSegment(Graph *graph, Segment &seg, std::vector<bool>& placed, std::vector<int> &chain) {
    int vert = graph->GetVertexCount();
    std::vector<unsigned char> colour(vert, 0);

    std::stack<Node> stack;
    int to, start = seg.contact_vertices[0], finish;
    Node v;
    stack.push(Node(start, -1));

    if (seg.vertices.size() == 2) {
        chain.push_back(seg.vertices[0]);
        chain.push_back(seg.vertices[1]);
        return;
    }

    // find chain using dfs
    DFSForChain(graph, chain, stack, colour, placed, seg, finish);
    // recover chain
    RecoverChain(chain, stack, colour, finish);
}

void MarkEdgesAsPlaced(int **matrix_graph, std::vector<Area> &areas, std::vector<int> &new_boarder,
                       std::vector<int> &chain, std::vector<bool>& placed, std::vector<int> &placed_vertices, int area_number) {
    int v;
    int contact = chain[0];
    for (int i = 0; i < areas[area_number].boarder.size(); ++i) {
        v = areas[area_number].boarder[i];
        if (v != contact) {
            new_boarder.push_back(v);
        }
        else {
            new_boarder.push_back(contact);
            int j;
            for (j = 1; j < chain.size() - 1; ++j) {
                matrix_graph[chain[j]][chain[j - 1]] = 2;
                matrix_graph[chain[j - 1]][chain[j]] = 2;

                new_boarder.push_back(chain[j]);
                placed[chain[j]] = true;
                placed_vertices.push_back(chain[j]);
                areas[area_number].vertices.insert(chain[j]);
            }
            matrix_graph[chain[j]][chain[j - 1]] = 2;
            matrix_graph[chain[j - 1]][chain[j]] = 2;

            new_boarder.push_back(contact);
        }
    }
}

void ChainWithOneContact(int **matrix_graph, std::vector<Area> &areas,
                         int area_number, std::vector<int> &chain, std::vector<bool>& placed, std::vector<int> &placed_vertices) {
    std::vector<int> new_boarder;
    MarkEdgesAsPlaced(matrix_graph, areas, new_boarder, chain, placed, placed_vertices, area_number);
    areas[area_number].boarder = new_boarder;

    Area new_area;
    for (int i = chain.size() - 1; i > 0; --i) {
        new_area.vertices.insert(chain[i]);
        new_area.boarder.push_back(chain[i]);
    }
    areas.push_back(new_area);
}

void AddVerticesAndEdgesToPlaced(int **matrix_graph, std::vector<int> &chain,
                                 std::vector<bool>& placed, std::vector<int> &placed_vertices) {
    int k;
    for (k = 1; k < chain.size() - 1; ++k) {
        placed_vertices.push_back(chain[k]);
        placed[chain[k]] = true;

        matrix_graph[chain[k]][chain[k - 1]] = 2;
        matrix_graph[chain[k - 1]][chain[k]] = 2;
    }
    matrix_graph[chain[k]][chain[k - 1]] = 2;
    matrix_graph[chain[k - 1]][chain[k]] = 2;
}

void MakeAreaOne(Area &area_one, std::vector<int> &chain, std::vector<Area> &areas, int area_number,
                 int &begin_pos, int &end_pos) {
    int v, begin_contact = chain[0], end_contact = chain[chain.size() - 1];
    bool found_begin = false, found_end = false;
    for (int i = 0; i < areas[area_number].boarder.size(); ++i) {
        v = areas[area_number].boarder[i];
        if (!found_begin && !found_end && v != begin_contact && v != end_contact) {
            area_one.boarder.push_back(v);
            area_one.vertices.insert(v);
        } else if (v == begin_contact && !found_end) {
            found_begin = true;
            begin_pos = i;
            for (int j = 0; j < chain.size(); ++j) {
                area_one.boarder.push_back(chain[j]);
                area_one.vertices.insert(chain[j]);
            }
        } else if (v == begin_contact && found_end) {
            begin_pos = i;
            found_begin = true;
        } else if (v == end_contact && !found_begin) {
            found_end = true;
            end_pos = i;
            for (int j = chain.size() - 1; j >= 0; --j) {
                area_one.boarder.push_back(chain[j]);
                area_one.vertices.insert(chain[j]);
            }
        } else if (v == end_contact && found_begin) {
            end_pos = i;
            found_end = true;
        } else if (found_begin && found_end) {
            area_one.boarder.push_back(v);
            area_one.vertices.insert(v);
        }
    }
}

void MakeAreaTwo(Area &area_two, std::vector<int> &chain, std::vector<Area> &areas,
                 int area_number, int begin_pos, int end_pos) {
    if (begin_pos < end_pos) {
        for (int i = begin_pos; i <= end_pos; ++i) {
            area_two.boarder.push_back(areas[area_number].boarder[i]);
            area_two.vertices.insert(areas[area_number].boarder[i]);
        }

        for (int i = chain.size() - 2; i > 0; --i) {
            area_two.boarder.push_back(chain[i]);
            area_two.vertices.insert(chain[i]);
        }
    } else {
        for (int i = end_pos; i <= begin_pos; ++i) {
            area_two.boarder.push_back(areas[area_number].boarder[i]);
            area_two.vertices.insert(areas[area_number].boarder[i]);
        }

        for (int i = 1; i < chain.size() - 1; ++i) {
            area_two.boarder.push_back(chain[i]);
            area_two.vertices.insert(chain[i]);
        }
    }
}

void PlaceChainOnSurface(int **matrix_graph, std::vector<Area> &areas,
                         int area_number, std::vector<int> &chain, std::vector<bool>& placed, std::vector<int> &placed_vertices) {
    if (chain[0] == chain[chain.size() - 1]) {
        ChainWithOneContact(matrix_graph, areas, area_number, chain, placed, placed_vertices);
        return;
    }

    AddVerticesAndEdgesToPlaced(matrix_graph, chain, placed, placed_vertices);

    Area area_one;
    int begin_pos = -1;
    int end_pos = -2;
    MakeAreaOne(area_one, chain, areas, area_number, begin_pos, end_pos);

    Area area_two;
    MakeAreaTwo(area_two, chain, areas, area_number, begin_pos, end_pos);

    areas[area_number] = area_one;
    areas.push_back(area_two);
}

bool IsConnectionComponentPlanar(Graph *graph, int **matrix_graph, int any_vertex, std::vector<bool>& checked, std::vector<bool>& placed) {
    int verts_count = graph->GetVertexCount();
    std::vector<int> component;
    std::vector<int> placed_vertices;
    std::vector<int> chain;
    BFSMarkingConnectionComponentToWhichVertexBelongs(matrix_graph, any_vertex, checked, verts_count,
                                                      component);
    std::vector<bool> segmentated(verts_count);
    std::vector<int> cycle;
    if (!FindCycle(graph, matrix_graph, any_vertex, cycle))
        return true;
    else {
        std::vector<Area> areas;
        PlaceCycleOnSurface(matrix_graph, placed, areas, cycle, placed_vertices);
        while (placed_vertices.size() <= component.size())// <=
        {
            bool any_segment = false;
            Segment best_segment = FindBestSegment(graph, matrix_graph, placed, areas, placed_vertices,
                                                   segmentated, any_segment);
            for (int i = 0; i < best_segment.vertices.size(); ++i)
                best_segment.hash_vertices.insert(best_segment.vertices[i]);
            if (any_segment) {
                if (best_segment.involved_areas.empty())
                    return false;
            } else
                return true;
            chain.clear();
            GetChainFromSegment(graph, best_segment, placed, chain);
            int area_number = best_segment.involved_areas[0];
            PlaceChainOnSurface(matrix_graph, areas, area_number, chain, placed, placed_vertices);
        }
    }
    return true;
}

bool IsGraphPlanar(int vert, int edge, std::vector<int>& from, std::vector<int>& to) {
    bool ans = true;

    int **matrix_graph = new int *[vert];
    for (int i = 0; i < vert; ++i) {
        matrix_graph[i] = new int[vert];
    }
    for (int i = 0; i < vert; ++i) {
        for (int j = 0; j < vert; ++j) {
            matrix_graph[i][j] = 0;
        }
    }
    auto *graph = new Graph(vert);
    BuildGraph(graph, edge, from, to, matrix_graph);

    RemoveBridgesFromGraph(graph, matrix_graph);

    std::vector<bool> checked(vert);
    std::vector<bool> placed(vert);
    for (int i = 0; i < vert; ++i) {
        checked[i] = false;
        placed[i] = false;
    }

    for (int i = 0; i < vert; ++i) {
        if (!checked[i]) {
            if (IsConnectionComponentPlanar(graph, matrix_graph, i, checked, placed)) {
                continue;
            } else {
                ans = false;
                break;
            }
        }
    }

    for (int i = 0; i < vert; ++i) {
        delete[] matrix_graph[i];
    }

    delete[] matrix_graph;
    delete graph;

    return ans;
}

int main() {
    int vert, edge;
    std::cin >> vert >> edge;

    std::vector<int> from(edge);
    std::vector<int> to(edge);

    for (int i = 0; i < edge; ++i) {
        std::cin >> from[i] >> to[i];
    }

    if (IsGraphPlanar(vert, edge, from, to))
        std::cout << "YES" << std::endl;
    else
        std::cout << "NO" << std::endl;

    return 0;
}