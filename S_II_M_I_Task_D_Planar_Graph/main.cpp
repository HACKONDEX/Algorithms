#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>
#include <queue>
#include <unordered_set>

class Graph
{
public:
    explicit Graph( int vert_count = 0 ) : vertex_count( vert_count ), edge_count( 0 ) {

        for( int i = 0; i < vert_count; ++i )
            graph.emplace_back();

    }
    ~Graph() = default;

    Graph( const Graph& ) = delete;
    Graph( Graph&& ) = delete;
    Graph& operator= ( const Graph& ) = delete;
    Graph& operator= ( Graph&& ) = delete;

    void add_vertex();
    void add_edge( int from, int to );
    void remove_edge( int from, int to );
    bool has_edge( int from, int to ) const;
    const std::vector<int>& get_next_vertices( int from ) const;
    int get_vertex_count();
    int get_edge_count();



private:

    int vertex_count;
    int edge_count;
    std::vector< std::vector<int> > graph;

};

void Graph::add_vertex()
{
    ++vertex_count;
    graph.emplace_back();
}

int Graph::get_edge_count()
{
    return edge_count;
}

void Graph::add_edge( int from, int to )
{
    graph[from].push_back( to );
    ++edge_count;
}

void Graph::remove_edge( int from, int to )
{
    std::vector< int > perm;
    for( int i = 0; i < graph[from].size(); ++i )
    {
        if( graph[from][i] != to )
            perm.push_back( graph[from][i] );
    }
    --edge_count;
    graph[from] = perm;
}

bool Graph::has_edge( int from, int to ) const
{
    for( int i = 0; i < graph[from].size(); ++i )
        if( graph[from][i] == to )
            return true;
    return false;
}

const std::vector<int>& Graph::get_next_vertices( int from ) const
{
    return graph[from];
}

int Graph::get_vertex_count()
{
    return vertex_count;
}


int min( int a, int b )
{
    return a < b ? a : b;
}

//// Solution

void create_graph( Graph* graph, int edge, int* from, int* to, int** matrix_graph )
{
    for( int i = 0; i < edge; ++i )
    {
        if( from[i] != to[i] )
        {
            graph->add_edge(from[i], to[i]);
            graph->add_edge(to[i], from[i]);
            matrix_graph[from[i]][to[i]] = 1;
            matrix_graph[to[i]][from[i]] = 1;
        }
    }
}

////little structure dfs
struct node
{
    int number;
    int parent;

    node( int v = - 1, int p = - 1 ) : number( v ), parent( p ){}
};

void do_if_not_coloured( std::stack< node >&stack, std::vector< int >& next_vertices, node& current_node,
        unsigned char* colour, int* in_time,  int* low )
{
    int to;
    for( int i = 0; i < next_vertices.size(); ++i )
    {
        to = next_vertices[i];
        if( to == current_node.parent )
            continue;
        else if( colour[to] == 0 )
            stack.push( node( to, current_node.number ) );
        else
            low[current_node.number] = min( low[current_node.number], in_time[to] );
    }
}

void do_if_coloured( std::vector< int >& next_vertices, node& current_node, int* in_time, int* low, int** matrix_graph )
{
    int to;
    for( int i = 0; i < next_vertices.size(); ++i )
    {
        to = next_vertices[i];
        if( current_node.parent != to  && low[to] < low[current_node.number] )
            low[current_node.number] = low[to];
    }
    for( int i = 0; i < next_vertices.size(); ++i )
    {
        to = next_vertices[i];
        if ( current_node.parent != to && low[to] > in_time[current_node.number])
            if( matrix_graph[to][current_node.number] == 1 )
            {
                matrix_graph[to][current_node.number] = - 1;
                matrix_graph[current_node.number][to] = - 1;
            }
    }
}

void find_bridges( Graph* graph, int** matrix_graph )
{
    int vertices = graph->get_vertex_count();
    auto* colour = new unsigned char[vertices];
    int* low = new int[vertices];
    auto* in_time = new int[vertices];
    for( int i = 0; i < vertices; ++i )
        colour[i] = 0;

    std::stack< node > stack;
    std::vector< int > next_vertices;
    int  timer = 0;
    node current_node( 0, 0 );

    for( int i = 0; i < vertices; ++i )
    {
        if( colour[i] == 0 )
        {
            stack.push( node( i, -1 ) );

            while( !stack.empty() )
            {
                current_node = stack.top();
                next_vertices = graph->get_next_vertices( current_node.number );
                if( colour[current_node.number] == 0 )
                {
                    colour[current_node.number] = 1;
                    low[current_node.number] = in_time[current_node.number] = timer++;
                    do_if_not_coloured( stack, next_vertices, current_node, colour, in_time, low );
                }
                else
                {
                    stack.pop();
                    colour[current_node.number] = 2;
                    do_if_coloured( next_vertices, current_node, in_time, low, matrix_graph );
                }
            }
        }
    }

    delete[] in_time;
    delete[] low;
    delete[] colour;
}

void remove_bridges_from_graph( Graph* graph, int** matrix_graph )
{
    int vertices = graph->get_vertex_count();
    find_bridges( graph, matrix_graph );

    for( int i = 0; i < vertices; ++i )
    {
        for( int j = i; j < vertices; ++j )
        {
            if( matrix_graph[i][j] == -1 )
            {
                matrix_graph[i][j] = 0;
                matrix_graph[j][i] = 0;
                graph->remove_edge( i, j );
                graph->remove_edge( j, i );
            }
        }
    }
}

void bfs_for_marking_connection_component_to_which_vertex_belongs( int** matrix_graph, int start, bool* checked,
        int vertices, std::vector< int >& component )
{
    std::queue< int > q;
    checked[start] = true;
    q.push( start );
    int vertex = 0;
    while( !q.empty() )
    {
        vertex = q.front();
        component.push_back( vertex );
        q.pop();
        for( int i = 0; i < vertices; ++i )
        {
            if( matrix_graph[vertex][i] == 1 && !checked[i] )
            {
                checked[i] = true;
                q.push( i );
            }
        }
    }
}

void recover_cycle( std::vector< int >& cycle, std::stack< node >& stack, int cycle_end, unsigned char* colour )
{
    cycle.push_back( cycle_end );
    node v;

    while( !stack.empty() )
    {
        v = stack.top();

        stack.pop();
        if( v.number == cycle_end )
        {
            break;
        }
        else
        {
            if( colour[v.number] == 1 )
            {
                cycle.push_back( v.number );
                colour[v.number] = 2;
            }
        }
    }
}

bool dfs_for_cycle( Graph* graph ,std::stack< node >& stack, int& cycle_end, unsigned char* colour )
{
    bool ans = false;
    node v;
    int to;
    bool found_cycle = false;
    std::vector< int > next_vertices;
    while( !stack.empty() )
    {
        v = stack.top();
        if( colour[v.number] == 0 )
        {
            colour[v.number] = 1;
            next_vertices = graph->get_next_vertices( v.number );
            for( int i = 0; i < next_vertices.size(); ++i )
            {
                to = next_vertices[i];
                if( to == v.parent )
                    continue;
                else if( colour[to] == 0 )
                    stack.push( node( to, v.number ) );
                else if( colour[to] == 1 )
                {
                    cycle_end = to;
                    ans = true;
                    found_cycle = true;
                    break;
                }
            }

            if( found_cycle )
                break;
        }
        else
        {
            stack.pop();
            colour[v.number] = 2;
        }
    }
    return ans;
}

//returns <false> if there is no cycle
bool find_any_cycle( Graph* graph, int** matrix_graph, int start, std::vector< int >& cycle )
{
    bool ans = false;
    int vertices = graph->get_vertex_count();

    auto* colour = new unsigned char[vertices];
    for( int i = 0; i < vertices; ++i )
        colour[i] = 0;

    std::stack< node > stack;

    //first vertex has no parent
    stack.push( node( start, -1 ) );
    node v; int to = 0;

    //how to remember cycle
    int cycle_end = - 1;

    ans = dfs_for_cycle( graph, stack, cycle_end, colour );

    if( ans )
        recover_cycle( cycle, stack, cycle_end, colour );

    delete[] colour;
    return ans;
}


struct Area
{
    std::vector< int > boarder;
    std::unordered_set< int > vertices;
};

struct Segment
{
    std::vector< int > vertices;
    std::unordered_set< int > hash_vertices;
    std::vector< int > contact_vertices;
    std::vector< int > involved_areas;
};

void place_cycle_on_surface( int** matrix_graph, bool* placed, std::vector< Area > & areas,
        std::vector< int >& cycle, std::vector< int >& placed_vertices )
{
    int j = 0;
    Area area_1;
    Area area_2;
    while( j < cycle.size() - 1 )
    {
        placed[cycle[j]] = true;
        placed_vertices.push_back( cycle[j] );

        area_1.boarder.push_back( cycle[j] );

        matrix_graph[cycle[j]][cycle[j + 1]] = 2;
        matrix_graph[cycle[j + 1]][cycle[j]] = 2;

        ++j;
    }
    placed[cycle[j]] = true;
    placed_vertices.push_back( cycle[j] );

    area_1.boarder.push_back( cycle[j] );
    area_2.boarder = area_1.boarder;
    std::reverse( area_2.boarder.begin(), area_2.boarder.end() );
    areas.push_back( area_1 );
    areas.push_back( area_2 );

    for( int j = 0; j < areas.size(); ++j )
        for( int i = 0; i < areas[j].boarder.size(); ++i )
            areas[j].vertices.insert( areas[j].boarder[i] );

    matrix_graph[cycle[j]][cycle[0]] = 2;
    matrix_graph[cycle[0]][cycle[j]] = 2;
}

void find_areas_for_segment( std::vector< Area >& areas, Segment& seg )
{
    int v;
    bool belongs_to_area = true;
    for (int i = 0; i < areas.size(); ++i)
    {
        belongs_to_area = true;
        for (int j = 0; j < seg.contact_vertices.size(); ++j) {
            v = seg.contact_vertices[j];
            if (areas[i].vertices.count(v) == 0) {
                belongs_to_area = false;
                break;
            }
        }
        if (belongs_to_area)
            seg.involved_areas.push_back(i);
    }
}

void is_segment_the_best( bool& any_segment, Segment& best_seg, Segment& seg )
{
    if( !any_segment )
    {
        any_segment = true;
        best_seg = seg;
    }
    else if( seg.involved_areas.size() < best_seg.involved_areas.size() )
        best_seg = seg;
}

void find_complex_segment( Graph* graph, std::vector< int >& comp_vertices, std::queue< int >& q, Segment& seg,
        int& current_component, int j, bool* segmentated, bool* placed  )
{
    segmentated[current_component] = true;
    seg.vertices.clear();
    seg.involved_areas.clear();
    seg.contact_vertices.clear();
    seg.vertices.push_back( current_component );
    seg.contact_vertices.push_back( current_component );

    q.push(comp_vertices[j]);
    int v, to = 0;
    std::vector< int > next_vertices;
    while (!q.empty()) {
        v = q.front();
        q.pop();
        if (placed[v] && v != current_component && !segmentated[v] ) {
            segmentated[v] = true;
            seg.vertices.push_back(v);
            seg.contact_vertices.push_back(v);
        } else if ( !placed[v] && !segmentated[v] ) {
            segmentated[v] = true;
            seg.vertices.push_back(v);
            next_vertices = graph->get_next_vertices(v);

            for (int i = 0; i < next_vertices.size(); ++i)
            {
                to = next_vertices[i];
                if (!segmentated[to] && to != current_component)
                    q.push(to);
            }
        }
    }
}

void complex_segment( Graph* graph, std::vector< int >& comp_vertices, std::vector< int >& placed_vertices,
        std::vector< Area >& areas, bool* segmentated, bool* placed, Segment& best_seg, bool& any_segment )
{
    int current_component;
    for( int i = 0; i < placed_vertices.size(); ++i )
    {
        current_component = placed_vertices[i];

        Segment seg;
        std::queue< int > q;
        comp_vertices = graph->get_next_vertices( current_component );
        for( int j = 0; j < comp_vertices.size(); ++j )
        {
            if( !segmentated[comp_vertices[j]] && !placed[comp_vertices[j]] )
            {
                find_complex_segment( graph, comp_vertices, q, seg, current_component, j, segmentated, placed );

                find_areas_for_segment( areas, seg );

                for( int i = 0; i < placed_vertices.size(); ++i )
                    segmentated[placed_vertices[i]] = false;

                is_segment_the_best( any_segment, best_seg, seg );

                if( best_seg.contact_vertices.size() == 0 )
                    return;
            }
        }
    }
}

void primitive_segment( Graph* graph, int** matrix_graph, std::vector< int >& comp_vertices, std::vector< int >& placed_vertices,
                        std::vector< Area >& areas, bool* segmentated, bool* placed, Segment& best_seg, bool& any_segment )
{
    Segment seg;
    int current_component, to;
    std::vector< int > next_vertices;
    for( int i = 0; i < placed_vertices.size(); ++i )
    {
        current_component = placed_vertices[i];
        segmentated[current_component] = true;
        next_vertices = graph->get_next_vertices( current_component );
        for( int i = 0; i < next_vertices.size(); ++i )
        {
            to = next_vertices[i];

            if( !segmentated[to] && placed[to] && matrix_graph[to][current_component] != 2 )
            {
                seg.vertices.clear();
                seg.contact_vertices.clear();
                seg.involved_areas.clear();
                seg.contact_vertices.push_back( to );
                seg.contact_vertices.push_back( current_component );
                seg.vertices = seg.contact_vertices;

                find_areas_for_segment( areas, seg );
                is_segment_the_best( any_segment, best_seg, seg );

                if( best_seg.contact_vertices.size() == 0 )
                    return;
            }
        }
    }
}

Segment find_the_best_segment( Graph* graph, int** matrix_graph, bool* placed, std::vector< Area >& areas,
        std::vector< int >& placed_vertices, bool* segmentated, bool& any_segment )
{
    Segment best_seg;
    int verts_count = graph->get_vertex_count();
    // check if the vert was taken as a segment node
    for( int i = 0; i < verts_count; ++i )
        segmentated[i] = false;

    std::vector< int > comp_vertices;

    // cycle from contact vertex// complex one
    complex_segment( graph, comp_vertices, placed_vertices,
            areas, segmentated, placed, best_seg, any_segment );

    if( best_seg.contact_vertices.size() == 0 )
        return best_seg;

    primitive_segment( graph, matrix_graph, comp_vertices, placed_vertices,
            areas, segmentated, placed, best_seg, any_segment );

    return best_seg;
}

void recover_chain( std::vector< int >& chain, std::stack< node >& stack, unsigned char* colour , int finish )
{
    node v;
    chain.push_back( finish );
    while( !stack.empty() )
    {
        v = stack.top();
        stack.pop();
        if( colour[v.number] == 1 )
        {
            colour[v.number] = 2;
            chain.push_back( v.number );
        }
    }
}

void dfs_for_chain( Graph* graph, std::vector< int >& chain, std::stack< node >& stack, unsigned char* colour,
        bool* placed, Segment& seg, int& finish )
{
    std::vector< int > next_vertices;
    node v;
    int to;
    bool break_from_dfs = false;
    while( !stack.empty() )
    {
        v = stack.top();
        if( colour[v.number] == 0 )
        {
            colour[v.number] = 1;
            next_vertices = graph->get_next_vertices( v.number );
            for( int i = 0; i < next_vertices.size(); ++i )
            {
                to = next_vertices[i];
                if( seg.hash_vertices.count( to ) > 0 )
                {
                    if( v.parent == to )
                        continue;
                    if (!placed[to] && colour[to] == 0)
                        stack.push(node(to, v.number));
                    else if (placed[to] && v.parent != -1 )
                    {
                        finish = to;
                        break_from_dfs = true;
                        break;
                    }
                }
            }
        }
        else if( colour[v.number] == 1 )
        {
            stack.pop();
            colour[v.number] = 2;
        }
        if( break_from_dfs )
            break;
    }
}

void get_chain_from_segment( Graph* graph, Segment& seg, bool* placed, std::vector< int >& chain )
{
    int vert = graph->get_vertex_count();
    auto* colour = new unsigned char[vert];
    for( int i = 0; i < vert; ++i )
        colour[i] = 0;

    std::stack< node > stack;
    int to, start = seg.contact_vertices[0], finish;
    node v;
    stack.push( node( start, -1 ) );

    if( seg.vertices.size() == 2 )
    {
        chain.push_back( seg.vertices[0] );
        chain.push_back( seg.vertices[1] );
        return;
    }

    // find chain using dfs
    dfs_for_chain( graph, chain, stack, colour, placed, seg, finish );

    // recover chain
    recover_chain( chain, stack, colour, finish );

    delete[] colour;
}

void mark_edges_as_placed( int** matrix_graph, std::vector< Area >& areas, std::vector<int>& new_boarder,
        std::vector< int >& chain, bool* placed, std::vector< int >& placed_vertices, int area_number )
{
    int v;
    int contact = chain[0];
    for (int i = 0; i < areas[area_number].boarder.size(); ++i)
    {
        v = areas[area_number].boarder[i];
        if( v != contact )
            new_boarder.push_back( v );
        else
        {
            new_boarder.push_back( contact );
            int j;
            for( j = 1; j < chain.size() - 1; ++j )
            {
                matrix_graph[chain[j]][chain[j - 1]] = 2;
                matrix_graph[chain[j - 1]][chain[j]] = 2;

                new_boarder.push_back( chain[j] );
                placed[chain[j]] = true;
                placed_vertices.push_back( chain[j] );
                areas[area_number].vertices.insert(chain[j]);
            }
            matrix_graph[chain[j]][chain[j - 1]] = 2;
            matrix_graph[chain[j - 1]][chain[j]] = 2;

            new_boarder.push_back( contact );
        }
    }
}

void chain_with_one_contact( int** matrix_graph, std::vector< Area >& areas,
        int area_number, std::vector< int >& chain, bool* placed, std::vector< int >& placed_vertices  )
{
    std::vector<int> new_boarder;

    mark_edges_as_placed( matrix_graph, areas, new_boarder, chain, placed, placed_vertices, area_number );

    areas[area_number].boarder = new_boarder;

    Area new_area;
    for( int i = chain.size() - 1; i > 0; --i )
    {
        new_area.vertices.insert( chain[i] );
        new_area.boarder.push_back( chain[i] );
    }

    areas.push_back( new_area );
}

void add_vertices_and_edges_to_placed_ones( int** matrix_graph, std::vector< int >& chain,
        bool* placed, std::vector< int >& placed_vertices )
{
    int k;
    for( k = 1; k < chain.size() - 1; ++k )
    {
        placed_vertices.push_back(chain[k]);
        placed[chain[k]] = true;

        matrix_graph[chain[k]][chain[k - 1]] = 2;
        matrix_graph[chain[k - 1]][chain[k]] = 2;
    }
    matrix_graph[chain[k]][chain[k - 1]] = 2;
    matrix_graph[chain[k - 1]][chain[k]] = 2;
}

void make_area_one( Area& area_one, std::vector< int >& chain, std::vector< Area >& areas, int area_number,
        int& begin_pos, int& end_pos )
{
    int v, begin_contact = chain[0], end_contact = chain[chain.size() - 1];
    bool found_begin = false, found_end = false;
    for( int i = 0; i < areas[area_number].boarder.size(); ++i )
    {
        v = areas[area_number].boarder[i];
        if( !found_begin && !found_end && v != begin_contact && v != end_contact )
        {
            area_one.boarder.push_back( v );
            area_one.vertices.insert( v );
        }
        else if( v == begin_contact && !found_end )
        {
            found_begin = true;
            begin_pos = i;
            for( int j = 0; j < chain.size(); ++j )
            {
                area_one.boarder.push_back( chain[j] );
                area_one.vertices.insert( chain[j] );
            }
        }
        else if( v == begin_contact && found_end )
        {
            begin_pos = i;
            found_begin = true;
        }
        else if( v == end_contact && !found_begin )
        {
            found_end = true;
            end_pos = i;
            for( int j = chain.size() - 1; j >= 0; --j )
            {
                area_one.boarder.push_back( chain[j] );
                area_one.vertices.insert( chain[j] );
            }
        }
        else if( v == end_contact && found_begin )
        {
            end_pos = i;
            found_end = true;
        }
        else if( found_begin && found_end )
        {
            area_one.boarder.push_back( v );
            area_one.vertices.insert( v );
        }
    }
}

void make_area_two( Area& area_two, std::vector< int >& chain, std::vector< Area >& areas,
        int area_number, int begin_pos, int end_pos)
{
    if( begin_pos < end_pos )
    {
        for( int i = begin_pos; i <= end_pos; ++i )
        {
            area_two.boarder.push_back( areas[area_number].boarder[i] );
            area_two.vertices.insert( areas[area_number].boarder[i] );
        }

        for( int i = chain.size() - 2; i > 0; --i )
        {
            area_two.boarder.push_back( chain[i] );
            area_two.vertices.insert( chain[i] );
        }
    }
    else
    {
        for( int i = end_pos; i <=begin_pos; ++i )
        {
            area_two.boarder.push_back( areas[area_number].boarder[i] );
            area_two.vertices.insert( areas[area_number].boarder[i] );
        }

        for( int i = 1; i < chain.size() - 1; ++i )
        {
            area_two.boarder.push_back( chain[i] );
            area_two.vertices.insert( chain[i] );
        }
    }
}

void place_chain_on_surface(  int** matrix_graph, std::vector< Area >& areas,
        int area_number, std::vector< int >& chain, bool* placed, std::vector< int >& placed_vertices )
{
    if( chain[0] == chain[chain.size() - 1] )
    {
        chain_with_one_contact( matrix_graph, areas, area_number, chain, placed, placed_vertices );
        return;
    }

    add_vertices_and_edges_to_placed_ones( matrix_graph, chain, placed, placed_vertices );

    Area area_one;
    int begin_pos = -1;
    int end_pos = -2;
    make_area_one( area_one, chain, areas, area_number, begin_pos, end_pos );

    Area area_two;
    make_area_two( area_two, chain, areas, area_number, begin_pos, end_pos );

    areas[area_number] = area_one;
    areas.push_back( area_two );
}

bool is_connection_component_planar( Graph* graph, int** matrix_graph, int any_vertex, bool* checked, bool* placed )
{
    int verts_count = graph->get_vertex_count();
    std::vector< int > component;
    std::vector< int > placed_vertices;
    std::vector< int > chain;
    bfs_for_marking_connection_component_to_which_vertex_belongs( matrix_graph, any_vertex, checked, verts_count, component );
    bool* segmentated = new bool[verts_count];
    std::vector< int > cycle;
    if( !find_any_cycle( graph, matrix_graph, any_vertex, cycle ) )
        return true;
    else
    {
        std::vector< Area > areas;
        place_cycle_on_surface( matrix_graph, placed, areas, cycle, placed_vertices );
        while( placed_vertices.size() <= component.size() )// <=
        {
            bool any_segment = false;
            Segment best_segment = find_the_best_segment(graph, matrix_graph, placed, areas, placed_vertices,
                                                         segmentated, any_segment);
            for (int i = 0; i < best_segment.vertices.size(); ++i)
                best_segment.hash_vertices.insert(best_segment.vertices[i]);
            if (any_segment)
            {
                if (best_segment.involved_areas.empty())
                    return false;
            }
            else
                return true;
            chain.clear();
            get_chain_from_segment(graph, best_segment, placed, chain);
            int area_number = best_segment.involved_areas[0];
            place_chain_on_surface( matrix_graph, areas, area_number, chain, placed, placed_vertices);
        }
    }

    delete[] segmentated;
    return true;
}

bool is_graph_planar( int vert, int edge, int* from, int* to )
{
    bool ans = true;

    int** matrix_graph = new int*[vert];
    for( int i = 0; i < vert; ++i )
        matrix_graph[i] = new int[vert];
    for( int i = 0; i < vert; ++i )
        for( int j = 0; j < vert; ++j )
            matrix_graph[i][j] = 0;
    auto* graph = new Graph( vert );
   create_graph( graph, edge, from, to, matrix_graph );

    remove_bridges_from_graph( graph, matrix_graph );

    bool* checked = new bool[vert];
    bool* placed = new bool[vert];
    for( int i = 0; i < vert; ++i )
    {
        checked[i] = false;
        placed[i] = false;
    }

    for( int i = 0; i < vert; ++i )
        if( !checked[i] )
        {
            if( is_connection_component_planar( graph, matrix_graph, i, checked, placed ) )
                continue;
            else
            {
                ans = false;
                break;
            }
        }

    for( int i = 0; i < vert; ++i )
        delete[] matrix_graph[i];

    delete[] matrix_graph;
    delete[] checked;
    delete[] placed;
    delete graph;

    return ans;
}

int main()
{
    int vert, edge;
    std::cin >> vert >> edge;

    int* from = new int[edge];
    int* to = new int[edge];

    for( int i = 0; i < edge; ++i )
        std::cin >> from[i] >> to[i];

    if( is_graph_planar( vert, edge, from, to ) )
        std::cout << "YES"<< std::endl;
    else
        std::cout << "NO" << std::endl;

    delete[] to;
    delete[] from;

    return 0;
}