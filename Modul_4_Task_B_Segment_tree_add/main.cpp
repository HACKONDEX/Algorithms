#include <iostream>
#include <vector>

const int minus_inf = -100000000;

class Segment_tree
{
private:
    struct Node
    {
        int max;
        int add;
        explicit Node( int m = minus_inf, int a = 0 ) : max( m ), add( a ) {}
    };

    int nodes_count;
    int vertices_count;
    std::vector< Node > nodes;

    static int get_min_degree_of_two( int n );
    void proselytize( int index, int left, int right );
    int max_( int left, int right, int index, int left_border, int right_border );
    void update_add( int left, int right, int index, int value, int left_border, int right_border );

public:

    Segment_tree( const int* data, int n );

    void add_value_on_segment( int left, int right, int value );
    int get_max_on_segment( int left, int right );

    ~Segment_tree() = default;
    Segment_tree( const Segment_tree& ) = delete;
    Segment_tree( Segment_tree&& ) = delete;
    Segment_tree& operator = ( const Segment_tree& ) = delete;
    Segment_tree& operator = ( Segment_tree&& ) = delete;
};

Segment_tree::Segment_tree( const int* data, int n )
{
    vertices_count = get_min_degree_of_two( n );
    nodes_count = 2 * vertices_count - 1;
    nodes.assign( nodes_count, Node() );
    for( int i = n - 1; i >= 0; --i )
        nodes[nodes_count - vertices_count + i].max = data[i];
    for( int i = nodes_count - vertices_count - 1; i >= 0; --i )
        nodes[i].max = std::max( nodes[ 2 * i + 1].max, nodes[2 * i + 2].max );
}

int Segment_tree::get_min_degree_of_two( int n )
{
    int value = 1;
    while( value < n )
        value *= 2;
    return value;
}

void Segment_tree::proselytize( int index, int left_border, int right_border )
{
    if( right_border - left_border == 0 )
        return;
    if( nodes[index].add == 0 )
        return;
    nodes[2 * index + 1].add += nodes[index].add;
    nodes[2 * index + 2].add += nodes[index].add;
    nodes[2 * index + 1].max += nodes[index].add;
    nodes[2 * index + 2].max += nodes[index].add;
    nodes[index].add = 0;
}

int Segment_tree::max_( int left, int right, int index, int left_border, int right_border )
{
    proselytize( index, left_border, right_border );
    if( left_border > right || left > right_border )
        return minus_inf;
    if( left_border >= left && right_border <= right )
        return nodes[index].max;
    int separation_index = ( left_border + right_border ) / 2;
    int left_sub_segment = max_( left, right, 2 * index + 1, left_border, separation_index );
    int right_sub_segment = max_( left, right, 2 * index + 2, separation_index + 1, right_border);
    return std::max( left_sub_segment, right_sub_segment );
}

void Segment_tree::update_add( int left, int right, int index, int value, int left_border, int right_border )
{
    proselytize( index, left_border, right_border );
    if( left_border > right || left > right_border)
        return;
    if( left_border >= left && right_border <= right )
    {
        nodes[index].max += value;
        nodes[index].add += value;
        return;
    }
    int separation_index = ( left_border + right_border ) / 2;
    update_add( left, right, 2 * index + 1, value, left_border, separation_index );
    update_add( left, right, 2 * index + 2, value, separation_index + 1, right_border );
    nodes[index].max = std::max( nodes[2 * index + 1].max, nodes[2 * index + 2].max );
}

int Segment_tree::get_max_on_segment(int left, int right)
{
    return max_( left, right, 0, 0, vertices_count - 1 );
}

void Segment_tree::add_value_on_segment(int left, int right, int value)
{
    update_add( left, right, 0, value, 0, vertices_count - 1 );
}

void process_requests(int n, int max_tickets_count, int requests_count, int* places, int* left, int* right, int* tickets )
{
    Segment_tree segment_tree( places, n - 1 );
    for( int i = 0; i < requests_count; ++i )
        if( segment_tree.get_max_on_segment( left[i], right[i] - 1 ) + tickets[i] <= max_tickets_count )
            segment_tree.add_value_on_segment( left[i], right[i] - 1, tickets[i] );
        else
            std::cout << i << " ";
}

int main()
{
    int n, requests_count, max_tickets_count;
    std::cin >> n;
    int* places = new int[n - 1];
    for( int i = 0; i < n - 1; ++i )
        std::cin >> places[i];
    std::cin >> max_tickets_count >> requests_count;
    int* left = new int[requests_count];
    int *right = new int[requests_count];
    int* tickets = new int[requests_count];
    for( int i = 0; i < requests_count; ++i )
        std::cin >> left[i] >> right[i] >> tickets[i];

    process_requests(n, max_tickets_count, requests_count, places, left, right, tickets );

    delete[] tickets;
    delete[] right;
    delete[] left;
    delete[] places;
    return 0;
}