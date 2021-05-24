/* Надо найти транзитивное замыкание графа. Мой алгоритм решения отличается от того, что предполагалось.
 * Лектор предполагал написать алгоритма Флойда Варшалла с оптимизацией при помощи битмасок.
 * Я написал алгоритм Пурдома, который сжимает граф к компонентам сильной связности, находит транзитивное замыкание для неё
 * и востанавливает ответ для исходного графа, его асимптотика автором не доказана но предполагается О( (V+E)^2 ). */
#include <iostream>
#include <vector>
#include <stack>

class Graph
{
public:
    explicit Graph( int vert_count = 0 );
    ~Graph();

    Graph( const Graph& ) = delete;
    Graph( Graph&& ) = delete;
    Graph& operator= ( const Graph& ) = delete;
    Graph& operator= ( Graph&& ) = delete;

    void add_vertex();
    void add_edge( int from, int to );
    bool has_edge( int from, int to ) const;
    const std::vector<int>& get_next_vertices( int from ) const;
    int get_vertex_count();

    void print_matrix() const;
    void dfs_transit_closer();
    void set_matrix_value( int i, int j, int value );
    int get_matrix_value( int i, int j );

private:
    std::vector< std::vector<int> > graph;
    std::vector< std::vector<int> > matrix;
    int vertex_count;

    void dfs( int s, int v, bool first_time );

};

Graph::Graph( int vert_count ) : vertex_count( vert_count ) {
    std::vector<int> vec(vertex_count, 0);
    for( int i = 0; i < vert_count; ++i )
    {
        graph.emplace_back();
        matrix.push_back(vec);
    }
}

Graph::~Graph()
{
    vertex_count = -1;
}

void Graph::add_edge( int from, int to )
{
    graph[from].push_back( to );
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


void Graph::print_matrix() const
{
    for(int i = 0; i < vertex_count; ++i )
    {
        for( int j = 0; j < vertex_count; ++j )
            std::cout << matrix[i][j];
        std::cout << std::endl;
    }
}

void Graph::dfs( int s, int v, bool first_time )
{
    if( !first_time )
        matrix[s][v] = 1;
    for( int i(0); i < graph[v].size(); ++i )
    {
        if( !matrix[s][graph[v][i]] )
            dfs( s, graph[v][i], false );
    }
}

void Graph::dfs_transit_closer()
{
    for( int i(0); i < vertex_count; ++i )
        dfs( i, i, true );
}

void Graph::set_matrix_value( int i, int j, int value )
{
    matrix[i][j] = value;
}

int Graph::get_matrix_value( int i, int j )
{
    return matrix[i][j];
}

//// Solution

void sort_vertices_by_out( Graph* graph, int* out )
{
    unsigned char* colour = new unsigned char[graph->get_vertex_count()];
    for( int i = 0; i < graph->get_vertex_count(); ++i)
        colour[i] = 0;

    int vertex = 0, time = 0;
    std::stack< int > stack;
    std::vector< int > next_vertices;
    for( int i = 0; i < graph->get_vertex_count(); ++i )
    {
        if( colour[i] == 0 )
        {
            stack.push( i );
            while( !stack.empty() )
            {
                vertex = stack.top();
                if( colour[vertex] == 0 )
                {
                    colour[vertex] = 1;
                   next_vertices = graph->get_next_vertices( vertex );
                    for( int i = 0; i < next_vertices.size(); ++i )
                        if( colour[next_vertices[i]] == 0 )
                            stack.push( next_vertices[i] );
                }
                else if( colour[vertex] == 1 )
                {
                    colour[vertex] = 2;
                    out[time++] = vertex;
                    stack.pop();
                }
                else
                    stack.pop();
            }
        }
    }

    delete[] colour;
}

void set_initial_values( int* connected_components, unsigned char* colour, int n )
{
    for( int i = 0; i < n; ++i )
    {
        connected_components[i] = 0;
        colour[i] = 0;
    }
    connected_components[n] = 0;
}

void break_into_connected_components( Graph* graph, int* connected_components, Graph* inverted_graph )
{
    int n = graph->get_vertex_count();
    int* sorted_by_out_time_vertices = new int[n];
    sort_vertices_by_out( graph, sorted_by_out_time_vertices );
    unsigned char* colour = new unsigned char[n];
    set_initial_values( connected_components, colour, n );

    int vertex = 0, number_of_component = -1;
    std::stack< int > stack;
    std::vector< int > next_vertices;
    for( int i = n - 1; i >= 0; --i )
    {
        if( colour[sorted_by_out_time_vertices[i]] == 0 )
        {
            ++number_of_component;
            stack.push( sorted_by_out_time_vertices[i] );
            while( !stack.empty() )
            {
                vertex = stack.top();
                if( colour[vertex] == 0 )
                {
                    colour[vertex] = 1;
                    next_vertices = inverted_graph->get_next_vertices( vertex );
                    for( int i = 0; i < next_vertices.size(); ++i )
                        if( colour[next_vertices[i]] == 0 )
                            stack.push( next_vertices[i] );
                }
                else if( colour[vertex] == 1 )
                {
                    colour[vertex] = 2;
                    connected_components[vertex] = number_of_component;
                    stack.pop();
                }
                else
                    stack.pop();
            }
        }
    }
    connected_components[graph->get_vertex_count()] = number_of_component + 1;
    delete[] sorted_by_out_time_vertices;
    delete[] colour;
}

void make_graph_connected_components_graph( Graph* const graph, Graph* connected_graph, int* connection_components ) {
    int n = graph->get_vertex_count(), to;
    std::vector< int > next_vertices;
    for ( int i = 0; i < n; ++i ) {
        next_vertices = graph->get_next_vertices( i );
        for ( int j = 0; j < next_vertices.size(); ++j ) {
            to = next_vertices[j];
            if ( connection_components[i] != connection_components[to] )
                if (!connected_graph->has_edge( connection_components[i], connection_components[to] ) )
                    connected_graph->add_edge( connection_components[i], connection_components[to] );
            else
                if (i != to)
                    if (!connected_graph->has_edge( connection_components[i], connection_components[to] ) )
                        connected_graph->add_edge( connection_components[i], connection_components[to] );
        }
    }
}

void recover_transition_of_initial_graph( Graph* graph, Graph* connected_graph, int* components )
{
    int n = graph->get_vertex_count();
    int h_comp, j_comp;
    for( int h = 0; h < n; ++h )
    {
        for( int j = 0; j < n; ++j )
        {
            h_comp = components[h];
            j_comp = components[j];
            if( h_comp == j_comp )
            {
                if( connected_graph->get_matrix_value( h_comp, j_comp ) == 1 )
                    graph->set_matrix_value( h, j, 1 );
                else if( h != j )
                    graph->set_matrix_value( h, j, 1 );
            }
            else
                if( connected_graph->get_matrix_value( h_comp, j_comp ) == 1 )
                    graph->set_matrix_value( h, j, 1 );
        }
    }
}

void make_graph_transitively_closed(Graph* graph, Graph* inverted_graph )
{
    int n = graph->get_vertex_count();
    int* connection_components = new int[n + 1];

    break_into_connected_components( graph, connection_components, inverted_graph );
    Graph* connected_components_graph = new Graph( connection_components[n] );

    make_graph_connected_components_graph( graph, connected_components_graph, connection_components );
    connected_components_graph->dfs_transit_closer();

    recover_transition_of_initial_graph( graph, connected_components_graph, connection_components );
    graph->print_matrix();

    delete connected_components_graph;
    delete[] connection_components;
}

int main()
{
    int n;
    std::cin >> n;
    std::string str;
    Graph* graph = new Graph( n );
    Graph* inverted_graph = new Graph( n );
    for( int i = 0; i < n; ++i )
    {
        std::cin >> str;
        for( int j = 0; j < n; ++j )
        {
            if( str[j] == '1' )
            {
                graph->add_edge( i, j );
                graph->set_matrix_value( i, j, 1 );
                inverted_graph->add_edge( j, i );
            }
        }
    }

    make_graph_transitively_closed( graph, inverted_graph );

    delete inverted_graph;
    delete graph;

    return 0;
}