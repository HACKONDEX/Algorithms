/* Из i й вершины можно попасть в i + 1 по цене a и в i^2 + 1 по цене b,
 * надо найти минимальную стоимость из x в y. Вершин не больше  M <= 10^6.
 * */

#include <iostream>
#include <vector>
#include <set>

//// List Graph

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
    void add_edge( int from, T to );
    bool has_edge( int from, T to ) const;
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
void Graph< T >::add_edge( int from, T to )
{
    graph[from].push_back( to );
}

template < class T >
bool Graph< T >::has_edge( int from, T to ) const
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

//// Problem solution

Graph< std::pair< int, int > >* make_graph( int a, int b, int m )
{
    auto* graph = new Graph< std::pair< int, int > >( m + 1 );

    for( long long i = 0; i < m ; ++i )
    {
        graph->add_edge( i, std::make_pair(  ( i + 1 ) % m , a ) );
        graph->add_edge( i, std::make_pair( ( i * i + 1 ) % m, b ) );
    }


    return graph;
}

void Dijkstra_with_set( Graph< std::pair< int, int > >* graph, int start, int* dist )
{
    std::vector< std::pair< int, int > > next_vertices;
    std::set< std::pair< int, int > > bin_tree;
    int current, to, weight;
    dist[start] = 0;
    bin_tree.insert( std::make_pair( dist[start], start ) );

    while( !bin_tree.empty() )
    {
        current = bin_tree.begin()->second;
        bin_tree.erase( bin_tree.begin() );

        next_vertices = graph->get_next_vertices( current );
        for( int i = 0; i < next_vertices.size(); ++i )
        {
            to = next_vertices[i].first;
            weight = next_vertices[i].second;
            if( dist[current] + weight < dist[to] )
            {
                bin_tree.erase( std::make_pair( dist[to], to ) );
                dist[to] = dist[current] + weight;
                bin_tree.insert( std::make_pair( dist[to], to ) );
            }
        }
    }
}

const int MAX = 1000000000;

int find_cost_of_teleportation( int a, int b, int m, int x, int y )
{
    if( x == y )
        return 0;

    int ans = 0;

    Graph< std::pair< int, int > >* graph = make_graph( a, b, m );

    int* dist = new int[m + 2];
    for( int i = 0; i <= m; ++i )
        dist[i] = MAX;

    Dijkstra_with_set( graph, x, dist );
    ans = dist[y];

    delete[] dist;
    delete graph;

    return ans;
}

int main() {

    int a,b,m,x,y;
    std::cin >> a >> b >> m >> x >> y;

    std::cout << find_cost_of_teleportation( a, b, m, x, y ) << std::endl;

    return 0;
}