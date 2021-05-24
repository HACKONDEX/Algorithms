#include <iostream>
#include <iomanip>
#include <vector>

//// Structure for edges

struct Edge
{
    int from;
    int to;
    double cost;

    Edge( int f, int t, double c ) : from( f ), to( t ), cost( c ){}
};

const double MIN = -2;

void make_edge_graph( std::vector< Edge >& graph ,int m, const int* from, const int* to, const double* cost )
{
    double c;
    for( int i = 0; i < m; ++i )
    {
        c = 100 - cost[i];
        c /= 100;
        graph.emplace_back( Edge( from[i], to[i], c ) );
        graph.emplace_back( Edge( to[i], from[i], c ) );
    }
}

void Ford_Bellman( std::vector< Edge >& graph ,int vert_count, double* dist, int start, int finish )
{

    int to, from;
    dist[start] = 1;


    for( int i = 0; i < vert_count; ++i )
    {
        for(auto & j : graph)
        {
            from = j.from;
            to = j.to;

            if( from == finish || to == start )
                continue;

            if( dist[from] != MIN && dist[to] < dist[from] * j.cost )
                dist[to] = dist[from] * j.cost;
        }
    }
}

double find_minimum_probability( int vert_count, int m, int start, int finish, const int* from, const int* to, const double* cost )
{
    std::vector< Edge > graph;
    make_edge_graph( graph, m, from, to, cost );

    auto* prob = new double[vert_count + 1];
    for( int i = 1; i <= vert_count; ++i )
        prob[i] = MIN;

    Ford_Bellman( graph, vert_count, prob, start, finish );

    double ans = 1 - prob[finish];

    delete[] prob;

    return ans;
}

int main() {

    std::setprecision( 6 );
    int n, m, s, f;
    std::cin >> n >> m >> s >> f;

    int* from = new int[m];
    int* to = new int[m];
    auto* cost = new double[m];

    for( int i = 0; i < m; ++i )
        std::cin >> from[i] >> to[i] >> cost[i];

    std::cout << find_minimum_probability( n, m, s, f, from, to, cost ) << std::endl;

    delete[] cost;
    delete[] to;
    delete[] from;
    return 0;
}