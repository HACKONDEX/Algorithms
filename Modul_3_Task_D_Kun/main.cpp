#include <iostream>
#include <vector>
#include <string>

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
    bool any_edge( int vertex ) const;
    void print() const;


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

template < class T >
void Graph< T >::print() const
{
    for( int i = 0; i < vertex_count; ++i )
    {
        std::cout << "( " << i << " ) -->";
        for( int j = 0; j < graph[i].size(); ++j )
            std::cout << " " << graph[i][j];
        std::cout << std::endl;
    }
}

template < class T >
bool Graph< T >::any_edge(int vertex) const
{
    return !graph[vertex].empty();
}

//// Solution

int get_int_matrix_bridge( int** int_bridge_matrix, int n, int m, std::string* bridge )
{
    int number = 1, share = 1;
    for (int i = 0; i < n; ++i)
    {
        share = ( i % 2 == 0 ) ? 1 : -1;
        for (int j = 0; j < m; ++j)
        {
            int_bridge_matrix[i][j] = (bridge[i][j] == '*') ? share * (number++) : 0;
            share *= -1;
        }
    }
    return number;
}

void make_graph( Graph< int >* graph, int n, int m, int** bridge )
{
    for( int i = 0; i < n; ++i )
    {
        for( int j = 0; j < m; ++j )
        {
            if( bridge[i][j] >= 0 )
                continue;
            if( i - 1 >= 0 && bridge[i - 1][j] != 0 )
                graph->add_edge( -bridge[i][j], bridge[i - 1][j] );
            if( i + 1 < n && bridge[i + 1][j] != 0 )
                graph->add_edge( -bridge[i][j], bridge[i + 1][j] );
            if( j - 1 >= 0 && bridge[i][j - 1] != 0 )
                graph->add_edge( -bridge[i][j], bridge[i][j - 1] );
            if( j + 1 < m && bridge[i][j + 1] != 0 )
                graph->add_edge( -bridge[i][j], bridge[i][j + 1] );
        }
    }
}

bool Kuhn( int vertex, Graph< int >* graph, std::vector< int >& pairs, bool* colour )
{
    if( colour[vertex] )
        return false;
    colour[vertex] = true;
    std::vector< int > next_vertices = graph->get_next_vertices( vertex );
    int to;
    for( int i = 0; i < next_vertices.size(); ++i )
    {
        to = next_vertices[i];
        if( pairs[to] == -1 || Kuhn( pairs[to], graph, pairs, colour ) )
        {
            pairs[to] = vertex;
            return true;
        }
    }
    return false;
}

void get_maximal_pairs( Graph< int >* graph, std::vector< int >& pairs )
{
    int size = pairs.size();
    bool* colour = new bool[size];
    for( int j = 1; j < size; ++j )
    {
        if( !graph->any_edge( j ) )
            continue;
        for( int i = 0; i < size; ++i )
            colour[i] = false;
        Kuhn( j, graph, pairs, colour );
    }

    delete[] colour;
}

long long count_minimal_time( std::vector< int >& pairs, int a, int b )
{
    long long answer;
    int count = 0, size = pairs.size();
    for( int i = 1; i< size; ++i )
    {
        if( pairs[i] != -1 )
            ++count;
    }
    answer = count * a + ( ( size - 1 ) - 2 * count ) * b;
    return answer;
}

long long get_minimal_time( int n, int m, int a, int b, std::string* bridge )
{
    int** int_bridge_matrix = new int*[n];
    for( int i = 0; i < n; ++i )
        int_bridge_matrix[i] = new int[m];

    long long answer = 0;
    int max_vertex_count = get_int_matrix_bridge( int_bridge_matrix, n, m, bridge );
    Graph< int >* graph = new Graph< int>( max_vertex_count );
    if( a >= b * 2 )
        answer = ( max_vertex_count - 1 ) * b;
    else
    {
        make_graph(graph, n, m, int_bridge_matrix);
        std::vector<int> maximal_pairs(max_vertex_count, -1);
        get_maximal_pairs(graph, maximal_pairs);
        answer = count_minimal_time( maximal_pairs, a, b );
    }

    delete graph;
    for( int i = 0; i < n; ++i )
        delete[] int_bridge_matrix[i];
    delete[] int_bridge_matrix;

    return answer;
}

int main()
{
    int n, m, a, b;
    std::cin >> n >> m >> a >> b;
    std::string* bridge = new std::string[n];
    for( int i = 0; i < n; ++i )
        std::cin >> bridge[i];

    std::cout << get_minimal_time( n, m, a, b, bridge ) << std::endl;

    delete[] bridge;
    return 0;
}