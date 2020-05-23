/* Если вам интерестно описание, вы знаете автора)
 * Кратко, это задача на топологическую сортировку.*/
#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>

class Graph
{
public:
    explicit Graph( int vert_count = 0 ) : vertex_count( vert_count ) {

        for( int i = 0; i < vert_count; ++i )
            graph.emplace_back( );

    }
    ~Graph() = default;

    Graph( const Graph& ) = delete;
    Graph( Graph&& ) = delete;
    Graph& operator= ( const Graph& ) = delete;
    Graph& operator= ( Graph&& ) = delete;

    void add_vertex();
    void add_edge( int from, int to );
    bool has_edge( int from, int to ) const;
    const std::vector<int>& get_next_vertices( int from ) const;
    int get_vertex_count();
    void sort_vertices();

private:

    int vertex_count;
    std::vector< std::vector<int> > graph;

};

void Graph::add_vertex()
{
    ++vertex_count;
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

void Graph::sort_vertices()//думал надо сортировать, из за описания
{
    for( int i = 0; i < vertex_count; ++i )
        std::sort( graph[i].begin(), graph[i].end() );
}

template < class T >
void swap( T& a, T& b )
{
    T c = a;
    a = b;
    b = c;
}

Graph* make_graph( int vertex_count, int edge_count, int* from, int* to )
{
    Graph* graph = new Graph( vertex_count );

    for( int i = 0; i < edge_count; ++i )
        graph->add_edge(from[i], to[i]);
//    graph->sort_vertices();

    return graph;
}

bool is_operation_executable( Graph* graph, int* calls_order )
{

    unsigned char* colour = new unsigned char[graph->get_vertex_count()];
    for( int i = 0; i < graph->get_vertex_count(); ++i )
        colour[i] = 0;

    std::stack< int > stack;
    std::vector< int > next_vertices;
    int zeros_count = 0, vertex = 0, index = 0;

    for( int i = 0; i < graph->get_vertex_count(); ++i )//dfs
    {
        if( colour[i] == 0 )
        {
            ++zeros_count;
            stack.push( i );

            while( !stack.empty() ) {
                vertex = stack.top();

                if ( colour[vertex] == 0)
                {
                    colour[vertex] = 1;
                    next_vertices = graph->get_next_vertices( vertex );

                    for( int i = next_vertices.size() - 1; i >= 0; --i )
                    {
                        if( colour[next_vertices[i]] == 0 )
                            stack.push( next_vertices[i] );
                        else if( colour[next_vertices[i]] == 1 )
                            return  false;
                    }
                }
                else if( colour[vertex] == 1 )
                {
                    colour[vertex] = 2;
                    calls_order[index++] = vertex;
                    stack.pop();
                }
                else
                {
                    stack.pop();
                }
            }
        }
    }

    if( zeros_count == 0 )
        return false;

    delete[] colour;



    return true;
}

void find_calls_order( int vertex_count, int edge_count, int* from, int* to )
{
    Graph* graph = make_graph( vertex_count, edge_count, from, to );

    int* calls_order = new int[vertex_count];//порядок номеров, для ответа YES
    for( int i = 0; i < vertex_count; ++i )
        calls_order[i] = 0;

    std::vector< int > right_order;

   if( is_operation_executable( graph, calls_order ) )
   {
       std::cout << "YES" << std::endl;

       for( int i = 0; i < vertex_count; ++i )
       {
           std::cout << calls_order[vertex_count - 1 - i] << " ";
           right_order.push_back( calls_order[vertex_count -1 - i] );
       }
           std::cout << std::endl;
   }
   else
       std::cout << "NO" << std::endl;

    delete[] calls_order;
    delete graph;
}

int main() {
    int policeman_count, channels_count;
    std::cin >> policeman_count >> channels_count;

    int *from = new int[channels_count];
    int *to = new int[channels_count];

    for (int i = 0; i < channels_count; ++i)
        std::cin >> from[i] >> to[i];

    find_calls_order( policeman_count, channels_count, from, to );

    delete[] from;
    delete[] to;

    return 0;
}