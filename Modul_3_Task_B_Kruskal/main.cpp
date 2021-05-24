#include <iostream>
#include <vector>
#include <algorithm>

//// Disjoint set union

class DisjointSetUnion
{
public:
    DisjointSetUnion( int n ) : vertex_count( n )
    {
        for( int i = 0; i < n; ++i )
            roots.push_back( std::make_pair( i, 1 ) );
    }
    ~DisjointSetUnion() = default;

    bool are_in_same_component( int v, int u );
    void join( int u, int v );
    int get_vertex_count() const;

private:
    std::vector< std::pair< int, int > > roots;
    int vertex_count;

    int root( int v );
};

int DisjointSetUnion::get_vertex_count() const
{
    return vertex_count;
}

bool DisjointSetUnion::are_in_same_component( int v, int u )
{
    return root( v ) == root( u );
}

int DisjointSetUnion::root( int v )
{
    if( v == roots[v].first )
        return v;
    int current_root = root( roots[v].first );
    roots[v].first  = current_root;
    return current_root;
}

void DisjointSetUnion::join( int u, int v )
{
    int head = root( u ), tail = root( v );
    if( roots[head].second >= roots[tail].second )
        std::swap( head, tail );
    roots[head].second += roots[tail].second;
    roots[tail].first = head;
}

//// Solution

struct Edge
{
    int left;
    int right;
    int weight;

    explicit Edge( int l = 0, int r = 0, int w = 0 ) : left( l ), right( r ), weight( w ) {}

    bool operator < ( const Edge& other ) const;
};

bool Edge::operator < ( const Edge& other ) const
{
    if( weight < other.weight )
        return true;
    else if( weight == other.weight )
    {
        if( left < other.left )
            return true;
        else if( left == other.left )
            return right < other.right;
        return false;
    }
    return false;
}

bool comparator( const Edge& first, const Edge& second )
{
    return first < second;
}

void get_sorted_edges( std::vector< Edge >& edges, int m, int* from, int* to, int* cost )
{
    for( int i = 0; i < m; ++i )
        edges.push_back( Edge( from[i], to[i], cost[i] ) );
    std::sort( edges.begin(), edges.end(), comparator );
}

long long Kruskal(std::vector< Edge >& edges, DisjointSetUnion& dsu )
{
    long long minimum_weight_sum = 0;
    int n = dsu.get_vertex_count(), u, v, marked = 0;
    bool* colour = new bool[n];
    for( int i = 0; i < n; ++i )
        colour[i] = false;
    for( int i = 0; i < edges.size(); ++i )
    {
        u = edges[i].left;
        v = edges[i].right;
        if( dsu.are_in_same_component( v, u ) )
            continue;
        minimum_weight_sum += edges[i].weight;
        dsu.join( u, v );
        marked += ( !colour[u] ? 1 : 0 ) + ( !colour[v] ? 1 : 0 );
        colour[u] = colour[v] = true;
        if( marked >= n - 1 )
            break;
    }
    delete[] colour;
    return minimum_weight_sum;
}

long long get_MST_weight( int n, int m, int* from, int* to, int* cost )
{
    std::vector< Edge > edges;
    get_sorted_edges( edges, m, from, to, cost );
    DisjointSetUnion dsu( n + 1 );
    return Kruskal(edges, dsu);
}

int main()
{
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
