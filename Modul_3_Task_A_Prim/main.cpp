#include <iostream>
#include <vector>
#include <set>

template < class T >
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
    void add_edge( int from, const T& to );
    bool has_edge( int from, const T& to ) const;
    const std::vector< T >& get_next_vertices( int from ) const;
    int get_vertex_count();


private:

    int vertex_count;
    std::vector< std::vector< T > > graph;

};

template < class T >
void Graph< T >::add_vertex()
{
    ++vertex_count;
}

template < class T >
void Graph< T >::add_edge( int from, const T& to )
{
    graph[from].push_back( to );
}

template < class T >
bool Graph< T >::has_edge( int from, const T& to ) const
{
    for( int i = 0; i < graph[from].size(); ++i )
        if( graph[from][i] == to )
            return true;
    return false;
}

template < class T >
const std::vector< T >& Graph< T >::get_next_vertices( int from ) const
{
    return graph[from];
}

template < class T >
int Graph< T >::get_vertex_count()
{
    return vertex_count;
}

//// Solution
void make_graph( Graph< std::pair< int, int > >* graph, int m, int* from, int* to, int* cost )
{
    for( int i = 0; i < m; ++i )
    {
        graph->add_edge( from[i], std::make_pair( cost[i], to[i] ) );
        graph->add_edge( to[i], std::make_pair( cost[i], from[i] ) );
    }
}

long long Prim( Graph< std::pair< int, int > >* graph )
{
    long long sum = 0;
    int vertex_count = graph->get_vertex_count() - 1,  current_vertex = 1;
    bool* colour = new bool[vertex_count + 1];
    for( int i = 0; i <= vertex_count; ++i )
        colour[i] = false;

    std::set< std::pair< int, int > > min_heap;
    std::vector< std::pair< int, int > > next_vertices;
    std::pair< int, int > current_edge;

    while( vertex_count > 0 )
    {
        colour[current_vertex] = true;
        next_vertices = graph->get_next_vertices( current_vertex );

        for( int i = 0; i < next_vertices.size(); ++i )
            min_heap.insert( next_vertices[i] );

        while( !min_heap.empty() )
        {
            current_edge = *min_heap.begin();
            min_heap.erase( min_heap.begin() );

            if( ! colour[current_edge.second] )
            {
                sum += current_edge.first;
                current_vertex = current_edge.second;
                break;
            }
        }
        --vertex_count;
    }

    delete[] colour;
    return sum;
}

long long get_MST_weight( int n, int m, int* from, int* to, int* cost )
{
    Graph< std::pair< int, int > >* graph = new Graph< std::pair< int, int > >( n + 1 );
    make_graph( graph, m, from, to, cost );
    long long answer = Prim( graph );
    delete graph;
    return answer;
}

int main() {
    int n,m;
    std ::cin >> n >> m;
    int* from = new int[m];
    int* to = new int[m];
    int* cost = new int[m];
    for( int i = 0; i < m; ++i )
        std::cin >> from[i] >> to[i] >> cost[i];

    std::cout << get_MST_weight( n, m, from, to, cost ) << std::endl;

    delete[] cost;
    delete[] to;
    delete[] from;
    return 0;
}