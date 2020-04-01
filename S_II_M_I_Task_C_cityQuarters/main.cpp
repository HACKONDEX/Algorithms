#include <iostream>
#include <vector>
#include <cstring>
#include <stack>
#include <algorithm>
#include <queue>

class Graph
{
public:
    explicit Graph( int vert_count = 0 ) : vertex_count( vert_count ) {

        for( int i = 0; i < vert_count; ++i )
            graph.emplace_back();

    }
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



private:

    int vertex_count;
    std::vector< std::vector<int> > graph;

};

Graph::~Graph()
{
    vertex_count = -1;
}

void Graph::add_vertex()
{
    ++vertex_count;
    graph.emplace_back();
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

//// Solution

void create_graph( Graph* graph, int edge_count, int* from, int* to )
{
    for( int i = 0; i < edge_count; ++i )
        graph->add_edge( from[i] - 1, to[i] - 1 );
}

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

//возвращает массив показывающий принадлежность вершин к компонентам связности
void break_into_connected_components( Graph* graph, int* connected_components, Graph* inverted_graph )
{
    int n = graph->get_vertex_count();
    int* sorted_by_out_time_vertices = new int[n];
    sort_vertices_by_out( graph, sorted_by_out_time_vertices );
    auto* colour = new unsigned char[n];
    set_initial_values( connected_components, colour, n );

    int vertex = 0, number_of_component = 0;
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
    connected_components[graph->get_vertex_count()] = number_of_component;
    delete[] sorted_by_out_time_vertices;
    delete[] colour;
}

void find_starts_ends( Graph* graph, int* connected_components, std::pair< int, int >* starts_ends )
{
    // first - количество пудей в вершину, second -количество путей из вершины
    for( int i = 0; i <= connected_components[graph->get_vertex_count()]; ++i )
    {
        starts_ends->first = 0;
        starts_ends->second = 0;
    }

    std::vector< int > next_vertices;

    for( int i = 0; i < graph->get_vertex_count(); ++i )
    {
        next_vertices = graph->get_next_vertices( i );
        for( int j = 0; j < next_vertices.size(); ++j )
        {
            if( connected_components[i] != connected_components[next_vertices[j]] )
            {
                ++starts_ends[connected_components[i]].second;// есть путь из вершины i
                ++starts_ends[connected_components[next_vertices[j]]].first;//есть путь в вершину next_vertices[j]
            }
        }
    }
}

void find_in_out_nodes_count( int& incoming, int& outcoming, std::pair< int, int >* starts_ends, int count  )
{
    for( int i = 1; i <= count; ++i )
    {
        if( starts_ends[i].first == 0 )

            ++incoming;
        if( starts_ends[i].second == 0 )
            ++outcoming;
    }
}

int find_number_of_streets_to_add( int offices_count, int streets_count, int* from, int* to )
{
    int ans = 0;

    Graph* graph = new Graph( offices_count );
    create_graph( graph, streets_count, from, to );
    Graph* inverted_graph = new Graph( offices_count );
    create_graph( inverted_graph, streets_count, to, from );

    int* connected_components = new int[offices_count + 1];
    break_into_connected_components( graph, connected_components, inverted_graph );
    auto* starts_ends = new std::pair< int, int >[connected_components[graph->get_vertex_count()] + 1];
    find_starts_ends( graph, connected_components, starts_ends );

    if( connected_components[offices_count] == 1 )
        ans = 0;
    else
    {
        int incoming = 0, outcoming = 0;
        find_in_out_nodes_count( incoming, outcoming, starts_ends, connected_components[offices_count] );
        ans = (incoming >= outcoming ? incoming : outcoming);
    }

    delete[] starts_ends;
    delete[] connected_components;
    delete graph;
    delete inverted_graph;

    return ans;
}

int main()
{
    int offices_count, streets_count;
    std::cin >> offices_count >> streets_count;

    int* from = new int[streets_count];
    int* to = new int[streets_count];

    for( int i = 0; i < streets_count; ++i )
        std::cin >> from[i] >> to[i];

    std::cout << find_number_of_streets_to_add( offices_count, streets_count, from, to ) << std::endl;

    delete[] from;
    delete[] to;

    return 0;
}