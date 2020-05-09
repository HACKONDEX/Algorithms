#include <iostream>
#include <vector>

using std::vector;

struct Edge
{
    int to;
    int capacity;
    int flow;
    int friend_edge_number;

    explicit Edge( int t = 0, int c = 0, int f = 0, int i = 0 ) : to( t ), capacity( c ), flow( f ), friend_edge_number( i ) {}
};

int s, t;

void make_graph( int n, int m, int* from, int* to, vector< vector< Edge> >& graph, vector< vector< Edge > >& indirect_graph )
{
    vector< Edge > a;
    for( int i = 0; i <= n; ++i )
    {
        graph.push_back( a );
        indirect_graph.push_back( a );
    }
    for( int i = 0; i < m; ++i )
    {
        if( from[i] == to[i] )
            continue;
        int size = graph[from[i]].size();
        int indirect_size = indirect_graph[to[i]].size();
        graph[from[i]].push_back( Edge(to[i], 1, 0, indirect_size) );
        indirect_graph[to[i]].push_back( Edge( from[i], 0, 0,  size) );
    }
}

int dfs( int vertex, int min_flow, vector< bool >& marked, vector< vector< Edge> >& graph, vector< vector< Edge > >& indirect_graph )
{
    if( marked[vertex] )
        return 0;
    marked[vertex] = true;
    if( vertex == t )
    {
        return min_flow;
    }
    for( int i = 0; i < graph[vertex].size(); ++i )
    {
        Edge u = graph[vertex][i];
        if( u.flow < u.capacity )
        {
            int new_delta = dfs(u.to, std::min( min_flow, u.capacity - u.flow ), marked, graph, indirect_graph );
            if( new_delta > 0 )
            {
                graph[vertex][i].flow += new_delta;
                indirect_graph[u.to][u.friend_edge_number].flow -= new_delta;
                return new_delta;
            }
        }
    }
    for( int i = 0; i < indirect_graph[vertex].size(); ++i )
    {
        Edge u = indirect_graph[vertex][i];
        if( u.flow < u.capacity )
        {
            int new_delta = dfs( u.to, std::min( min_flow, u.capacity - u.flow ), marked, graph, indirect_graph );
            if( new_delta > 0 )
            {
                indirect_graph[vertex][i].flow += new_delta;
                graph[u.to][u.friend_edge_number].flow -= new_delta;
                return new_delta;
            }
        }
    }
    return 0;
}

int dfs_for_answer( int vertex, vector<int>& path, vector< bool >& marked, vector< vector < Edge > >& graph )
{
    if( marked[vertex] )
        return 0;
    marked[vertex] = true;
    if( vertex == t )
        return t;
    for( int i = 0; i < graph[vertex].size(); ++i )
    {
        Edge u = graph[vertex][i];
        if( u.capacity == u.flow )
        {
            int new_vertex = dfs_for_answer( u.to, path, marked, graph );
            if( new_vertex > 0 )
            {
                path.push_back(new_vertex);
                graph[vertex][i].capacity += 1;
                return vertex;
            }
        }
    }
    return 0;
}

int Ford_Falkerson( vector< bool >& marked, vector< vector < Edge > >& graph, vector< vector< Edge > >& indirect_graph )
{
    int answer = 0;
    while( true )
    {
        marked.assign( marked.size(), false );
        int flow = dfs( s, 1000000, marked, graph, indirect_graph );
        if( flow > 0 )
            answer += flow;
        else
            break;
    }
    return answer;
}

void print_path( vector< bool >& marked, vector< vector < Edge > >& graph )
{
    vector< int > path;
    marked.assign( marked.size(), false );
    dfs_for_answer( s, path, marked, graph );
    path.push_back(s);
    for( int i = path.size() - 1; i >= 0; --i )
        std::cout << path[i] << " ";
    std::cout << std::endl;
}

void solve( int n, int m, int s, int t, int* to, int* from )
{
    vector< vector < Edge > > graph;
    vector< vector < Edge > > indirect_graph;
    make_graph( n, m, from, to, graph, indirect_graph );
    vector< bool > marked( n + 1, false );
    if( Ford_Falkerson( marked, graph, indirect_graph ) < 2 )
    {
        std::cout << "NO" << std::endl;
        return;
    }
    std::cout << "YES" << std::endl;
    print_path( marked, graph );
    print_path( marked, graph );
}



int main() {
    int n, m;
    std::cin >> n >> m >> s >> t;
    int* from = new int[m];
    int* to = new int[m];

    for(int i = 0; i < m; ++i)
        std::cin >> from[i] >> to[i];

    solve(n, m, s, t, to, from);

    delete[] to;
    delete[] from;
    return 0;
}