#include <iostream>
#include <vector>

void build_graph_tree( int n, std::vector < std::vector<int> >& graph_tree )
{
    int parent;
    for( size_t i = 1; i < n; ++i )
    {
        std::cin >> parent;
        graph_tree[parent].emplace_back(i);
    }
}

void reprocessing_dfs( std::vector< std::vector< int > >& graph_tree, std::vector< std::vector< int > >& jump,
                       std::vector< int >& in_time, std::vector< int >& out_time,  int vertex, int parent, int& timer )
{
    in_time[vertex] = ++timer;
    jump[vertex][0] = parent;

    for( int i = 1; i < jump[vertex].size(); ++i )
        jump[vertex][i] = jump[jump[vertex][i - 1]][i - 1];

    for( int i = 0; i < graph_tree[vertex].size(); ++i )
        if( graph_tree[vertex][i] != parent )
            reprocessing_dfs( graph_tree, jump, in_time, out_time,
                              graph_tree[vertex][i] , vertex, timer );

    out_time[vertex] = ++timer;

}

int find_log( int n )
{
    int log = 0, value = 1;
    while( value <= n )
    {
        value <<= 1;
        ++log;
    }
    return log;

}

void set_jump_size( int log, std::vector < std::vector<int> >& jump )
{
    for( int i = 0; i < jump.size(); ++i )
        jump[i].resize( log + 1, 0 );
}

int get_least_common_ancestor( int u, int v, std::vector< int >& in_time, std::vector< int >& out_time, std::vector< std::vector < int > >& jump )
{
    if( in_time[u] <= in_time[v] && out_time[u] >= out_time[v] )
        return u;
    if( in_time[v] <= in_time[u] && out_time[v] >= out_time[u] )
        return v;
    for( int i = jump[u].size() - 1; i >= 0; --i )
    {
        if( ! ( in_time[jump[u][i]] <= in_time[v] && out_time[jump[u][i]] >= out_time[v] ) )
            u = jump[u][i];
    }
    return jump[u][0];
}

void process_requests( int n, int m, std::vector< int >& in_time, std::vector< int >& out_time, std::vector< std::vector < int > >& jump )
{
    long long a1, a2, x, y, z;
    std::cin >> a1 >> a2 >> x >> y >> z;
    long long sum, v = sum = get_least_common_ancestor( a1, a2, in_time, out_time, jump );
    while ( m >= 2 )
    {
        a1 =( x * a1 + y * a2 + z ) % n;
        a2 = ( x * a2 + y * a1 + z ) % n;
        v = get_least_common_ancestor( ( a1 + v ) % n, a2, in_time, out_time, jump );
        sum += v;
        --m;
    }
    std::cout << sum << std::endl;
}

void process_LCA_requests( int n, int m )
{
    std::vector < std::vector<int> > graph_tree( n ), jump(n);
    build_graph_tree( n, graph_tree );
    int log_n = find_log( n ), timer = 0;
    set_jump_size( log_n, jump );
    std::vector< int > in_time( n ), out_time( n );
    reprocessing_dfs( graph_tree, jump, in_time, out_time, 0, 0, timer );
    process_requests( n, m, in_time, out_time, jump );
}

int main()
{
    int n, m;
    std::cin >> n >> m;
    process_LCA_requests( n, m );
    return 0;
}