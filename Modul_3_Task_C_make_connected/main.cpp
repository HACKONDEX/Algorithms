#include <iostream>
#include <vector>
#include <algorithm>

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

//// Edge structure

struct Edge
{
    int left;
    int right;
    unsigned long long weight;

    explicit Edge( int l = 0, int r = 0, unsigned long long w = 0 ) : left( l ), right( r ), weight( w ) {}

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

//// Kruskal

unsigned long long Kruskal(std::vector< Edge >& edges, DisjointSetUnion& dsu )
{
    unsigned long long minimum_weight_sum = 0;
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

//// Solution

void make_sorted_edges( std::vector< Edge >& edges ,int n, int m, unsigned long long* a, unsigned long long* cost,  int* from, int* to )
{
    unsigned long long minimum = 1000000000009;
    int min_number = 1;
    for( int i = 0; i < n; ++i )
    {
        if( a[i] < minimum )
        {
            minimum = a[i];
            min_number = i + 1;
        }
    }

    for( int i = 0; i < n; ++i )
    {
        if( i + 1 != min_number)
            edges.push_back(Edge(min_number, i + 1, minimum + a[i]));
    }
    for( int i = 0; i < m; ++i )
        edges.push_back( Edge( from[i], to[i], cost[i] ) );
    std::sort( edges.begin(), edges.end(), comparator );
}

unsigned long long get_MST_weight( int n, int m, unsigned long long* a, unsigned long long* cost,  int* from, int* to )
{
    std::vector< Edge > edges;
    make_sorted_edges( edges, n, m, a, cost, from, to );
    DisjointSetUnion dsu( n + 1 );
    return Kruskal( edges, dsu );
}

int main()
{
    int n, m;
    std::cin >> n >> m;
    auto* a = new unsigned long long[n];
    auto* special_cost = new unsigned long long[m];
    int* from = new int[m];
    int* to = new int[m];
    for( int i = 0; i < n; ++i )
        std::cin >> a[i];
    for( int i = 0; i < m; ++i )
        std::cin >> from[i] >> to[i] >> special_cost[i];

    std::cout << get_MST_weight( n, m, a, special_cost, from, to ) << std::endl;

    delete[] to;
    delete[] from;
    delete[] special_cost;
    delete[] a;
    return 0;
}