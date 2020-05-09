#include <iostream>
#include <vector>
#include <string>

using std::vector;

void make_graph( int n, std::string* matrix, vector< vector< int > >& graph )
{
    for( int i = 0; i < n; ++i )
    {
        for( int j = 0; j < n; ++j )
        {
            if( matrix[i][j] == '1' )
                graph[i].push_back( 1 );
            else
                graph[i].push_back( 0 );
        }
    }
}

void Shtor_Wagner( int n, vector< vector< int > >& graph, vector< int >& minimal_cut )
{
    int minimal_cut_cost = 625000;
    vector< vector< int > > compressed_vertex(n);
    for( int i = 0; i < n; ++i )
        compressed_vertex[i].assign( 1, i );
    vector< int > weight_to_the_set( n, 0 );
    vector< bool > exist( n, true );
    vector< bool > is_in_set( n, false );
    for( int phase = 0; phase < n - 1; ++phase )
    {
        is_in_set.assign( n, false );
        weight_to_the_set.assign( n, 0 );
        for( int k = 0, penul; k < n - phase; ++k )
        {
            int current_vertex = -1;
            for( int i = 0; i < n; ++i )
                if( exist[i] && !is_in_set[i] && ( current_vertex == -1 || weight_to_the_set[current_vertex] < weight_to_the_set[i] ) )
                    current_vertex = i;
            if( k == n - phase - 1 )
            {
                if( weight_to_the_set[current_vertex] < minimal_cut_cost )
                {
                    minimal_cut_cost = weight_to_the_set[current_vertex];
                    minimal_cut = compressed_vertex[current_vertex];
                }
                compressed_vertex[penul].insert(compressed_vertex[penul].end(), compressed_vertex[current_vertex].begin(), compressed_vertex[current_vertex].end());
                for ( int i=0; i<n; ++i )
                    graph[penul][i] = graph[i][penul] += graph[current_vertex][i];
                exist[current_vertex] = false;
            }
            else
            {
                is_in_set[current_vertex] = true;
                for( int i = 0; i < n; ++i )
                    weight_to_the_set[i] += graph[current_vertex][i];
                penul = current_vertex;
            }
        }
    }
}

void recover_second_set( int n, vector< int >& minimal_cut, vector< int >& second_set )
{
    vector< bool > colour( n, false );
    for( int i = 0; i < minimal_cut.size(); ++i )
        colour[minimal_cut[i]] = true;
    for( int i = 0; i < colour.size(); ++i )
        if( !colour[i] )
            second_set.push_back( i );
}

void print_the_set( vector< int >& set )
{
    for( int i = 0; i < set.size(); ++i )
        std::cout << set[i] + 1 << " ";
    std::cout << std::endl;
}

void minimal_cut( int n, std::string* matrix )
{
    vector< vector < int > > graph(n);
    make_graph( n, matrix, graph );
    vector< int > minimal_cut, second_set;
    Shtor_Wagner( n, graph, minimal_cut );
    recover_second_set( n, minimal_cut, second_set );
    print_the_set( minimal_cut );
    print_the_set( second_set );
}

int main() {
    int n;
    std::cin >> n;
    std::string* matrix = new std::string[n];
    for( int i = 0; i< n; ++i )
        std::cin >> matrix[i];

    minimal_cut( n, matrix );

    delete[] matrix;
    return 0;
}