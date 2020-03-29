/* Дан взевшанный граф из n <= 300 вершин, m <= 10^5 ребер, нужно найти путь минимальной длины из вершины s в вершину f
 * который проходит максимум через k ребер. W[i] <= 10^6 */

#include <iostream>
#include <vector>

//// Solution

const long long MAX = 1000000000000000;

//// Edge graph structure

struct Edge
{
    int from;
    int to;
    int cost;
};

Edge* make_edge_graph( int edge_count, const int* from, const int* to, const int* cost )
{
    Edge* graph = new Edge[edge_count];

    for( int i = 0; i < edge_count; ++i )
    {
        graph[i].from = from[i];
        graph[i].to = to[i];
        graph[i].cost = cost[i];
    }

    return graph;
}

long long* make_dist_with_infinity( int ver_count )
{
    auto* dist = new long long[ver_count + 1];

    for( int i = 0; i <= ver_count; ++i )
        dist[i] = MAX;

    return dist;
}

void Ford_Bellman_modified( long long* dist, Edge* graph , int vert_count, int edge_count, int start, int max_verts )
{
    bool* colour = new bool[vert_count + 1];

    int to, from, cost;
    dist[start] = 0;
    for( int i = 0; i < max_verts; ++i )
    {
        for( int j = 0; j <= vert_count; ++j )
            colour[j] = false;

        for( int j = 0; j < edge_count; ++j )
        {
            from = graph[j].from;
            if( dist[from] < MAX && !colour[from] )
            {
                to = graph[j].to;
                cost = graph[j].cost;
                if( dist[to] > dist[from] + cost )
                {
                    dist[to] = dist[from] + cost;
                    colour[to] = true;
                }
            }
        }
    }

    delete[] colour;
}

long long find_minimum_cost( int vert_count, int edge_count, int max_vert, int start, int finish, const int* from, const int* to, const int* cost )
{
    if( start == finish )
        return 0;

    Edge* graph = make_edge_graph( edge_count, from, to, cost );
    long long* dist = make_dist_with_infinity( vert_count );

    Ford_Bellman_modified( dist, graph, vert_count, edge_count, start, max_vert );

    long long ans = dist[finish] < MAX ? dist[finish] : -1;

    delete[] dist;
    delete[] graph;

    return ans;
}

int main()
{
    int n, m, k, s, f;
    std::cin >> n >> m >> k >> s >> f;

    int* from = new int[m];
    int* to = new int[m];
    int* cost = new int[m];
    for( int i = 0; i < m; ++i )
        std::cin >> from[i] >> to[i] >> cost[i];

    std::cout << find_minimum_cost( n, m, k, s, f, from, to, cost ) << std::endl;

    delete[] cost;
    delete[] to;
    delete[] from;

    return 0;
}