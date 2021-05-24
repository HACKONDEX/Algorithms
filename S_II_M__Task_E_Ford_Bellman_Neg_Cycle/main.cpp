/* Оригинальную задачу сложно нормально описать, надо вывести есть ли отрицательный цикл в смысле задачи  */

#include <iostream>
#include <vector>
#include <iomanip>

const double MINIM = -1000000;
struct Edge
{
    int from_vertex;
    int to_vertex;
    double tarif;
    double commis;

    explicit Edge( int from = 0, int ver = 0, double t = 0, double c = 0) : from_vertex(from), to_vertex(ver), tarif(t), commis(c) {}
};

template < typename T >
void copy( T* this_one, T* into_this, int size )
{
    for( int i = 0; i < size; ++i )
        into_this[i] = this_one[i];
}

void make_graph_make_edges( Edge* edges, int m, int* from, int* to, double* a_b, double* b_a, double* c_a_b, double* c_b_a )
{
    int i = 0;
    for( int j = 0; j < m; ++j )
    {
        edges[i++] = Edge( from[j], to[j], a_b[j], c_a_b[j] );
        edges[i++] = Edge( to[j], from[j], b_a[j], c_b_a[j] );
    }
}

void Ford_Bellman_negative_cicle( Edge* edges, int n, int edges_count, double* d_last ) {
    for ( int j(0); j < n; ++j ) {
        for ( int i(0); i <edges_count; ++i ) {
            int from = edges[i].from_vertex;
            if (d_last[from] == MINIM)
                continue;
            int to = edges[i].to_vertex;
            double tarif = edges[i].tarif;
            double commis = edges[i].commis;
            if (d_last[to] == MINIM || tarif * (d_last[from] - commis) > d_last[to]) {
                d_last[to] = tarif * (d_last[from] - commis);
            }
        }
    }
}

void get_answer( double* d_next, double* d_last, int n )
{
    for( int i(1); i <= n; ++i )
    {
        if( d_next[i] > d_last[i] )
        {
            std::cout << "YES" << std::endl;
            return;
        }
    }
    std::cout << "NO" << std::endl;
}


void can_we_increase_money( int n, int m, int s, double value, int* from, int* to, double* from_tarif_to, double* from_commission_to, double* to_tarif_from, double* to_commission_from  )
{
    Edge* edges = new Edge[2 * m];

    make_graph_make_edges( edges, m, from, to, from_tarif_to, to_tarif_from, from_commission_to, to_commission_from );

    double* d_last = new double[n + 1];
    double* d_next = new double[n + 1];
    for( int i(1); i <= n; ++i )
        d_last[i] = MINIM;

    d_last[s] = value;
    Ford_Bellman_negative_cicle( edges, n, 2 * m, d_last );
    copy( d_next, d_last, n + 1 );
    Ford_Bellman_negative_cicle( edges, 1, 2 * m, d_next );
    get_answer( d_next, d_last, n );
    delete[] d_next;
    delete[] d_last;
    delete[] edges;
}

int main()
{
    int n, m, s;
    double value;
    std::setprecision(8);
    std::cin >> n >> m >> s >> value;
    int* from = new int[m];
    int* to = new int[m];
    double* from_tarif_to = new double[m];
    double* from_commission_to = new double[m];
    double* to_tarif_from = new double[m];
    double* to_commission_from = new double[m];
    for( int i(0); i < m; ++i )
    {
        std::cin >> from[i] >> to[i];
        std::cin >> from_tarif_to[i] >> from_commission_to[i];
        std::cin >> to_tarif_from[i] >> to_commission_from[i];

    }

    can_we_increase_money( n, m, s, value, from, to, from_tarif_to, from_commission_to, to_tarif_from, to_commission_from );

    delete[] to_commission_from;
    delete[] from_commission_to;
    delete[] from_tarif_to;
    delete[] to_tarif_from;
    delete[] from;
    delete[] to;

    return 0;
}