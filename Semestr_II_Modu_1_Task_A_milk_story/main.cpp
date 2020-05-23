#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <queue>

class Graph
{
public:
    explicit Graph( int vert_count = 0 ) : vertex_count( vert_count ) {

        for( int i = 0; i < vert_count; ++i )
        {
            graph.emplace_back( );
        }

    }
    ~Graph();

    Graph( const Graph& ) = delete;
    Graph( Graph&& ) = delete;
    Graph& operator= ( const Graph& ) = delete;
    Graph& operator= ( Graph&& ) = delete;

    void add_vertex();
    void add_edge( int from, int to );
    bool has_edge( int from, int to ) const;
    const std::vector<int>& get_next_vertices( int from ) const;
    int get_vertex_count();


private:

    int vertex_count;
    std::vector< std::vector<int> > graph;

};

Graph::~Graph()
{
    vertex_count = -1;
}

void Graph::add_vertex()
{
    ++vertex_count;
}

void Graph::add_edge( int from, int to )
{
    graph[from].push_back( to );
}

bool Graph::has_edge( int from, int to ) const
{
    for( int i = 0; i < graph[from].size(); ++i )
        if( graph[from][i] == to )
            return true;
    return false;
}

const std::vector<int>& Graph::get_next_vertices( int from ) const
{
    return graph[from];
}

int Graph::get_vertex_count()
{
    return vertex_count;
}

Graph* create_graph( int vertex_count, int edge_count, int* from, int* to )
{
    Graph* graph = new Graph( vertex_count );

    for( int i = 0; i < edge_count; ++i ) {
        graph->add_edge(from[i], to[i]);
        graph->add_edge(to[i], from[i]);
    }
    return  graph;
}

struct pair
{
    int vertex;
    int dist;

    explicit pair( int v = -1, int d = -1 ) : vertex( v ), dist( d ){}
};

void bfs_for_dists( int start, int* distances, Graph* graph )
{
    std::queue< pair > queue;

    bool* colour = new bool[graph->get_vertex_count()];
    std::memset( colour,false, graph->get_vertex_count() );

    pair top;
    colour[start] = true;
    queue.push( pair( start, 0 ) );
    std::vector< int > neighbours;

    while( !queue.empty() )
    {
        top = queue.front();
        queue.pop();
        distances[top.vertex] = top.dist;

        neighbours = graph->get_next_vertices( top.vertex );
        for( int i =  0; i < neighbours.size(); ++i )
            if( !colour[neighbours[i]] )
            {
                colour[neighbours[i]] = true;
                queue.push(pair(neighbours[i], top.dist + 1));
            }

    }

    delete[] colour;

}

int find_min_edges_count( int vertex_count, int edge_count, int* from, int* to, int leon, int matilda, int milk )
{
    int ans = 0;

    Graph* graph = create_graph( vertex_count, edge_count, from, to );

    int* leon_dist = new int[vertex_count];
    int* matilda_dist = new int[vertex_count];
    int* milk_dist = new int[vertex_count];

    bfs_for_dists( leon, leon_dist, graph );
    bfs_for_dists( matilda, matilda_dist, graph );
    bfs_for_dists( milk, milk_dist, graph );



    for( int i = 1; i < vertex_count; ++i )
    {
        if( i==1 || leon_dist[i] + matilda_dist[i] + milk_dist[i] <= ans )
            ans = leon_dist[i] + matilda_dist[i] + milk_dist[i];
    }

    delete[] leon_dist;
    delete[] matilda_dist;
    delete[] milk_dist;

    delete graph;


    return ans;
}

int main() {

    int vertex_count, edge_count, leon_pos, matilda_pos, milk_pos;
    std::cin >> vertex_count >> edge_count >> leon_pos >> matilda_pos >> milk_pos;
    ++vertex_count;

    int* from = new int[edge_count];
    int* to = new int[edge_count];

    for( int i = 0; i < edge_count; ++i )
        std::cin >> from[i] >> to[i];

    std::cout << find_min_edges_count( vertex_count, edge_count, from, to, leon_pos, matilda_pos, milk_pos ) << std::endl;

    delete[] from;
    delete[] to;




    return 0;
}