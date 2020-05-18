#include <iostream>
#include <vector>

const int inf = 100000000;

class Segment_tree
{
private:
    struct Node
    {
        int min;
        int assign;
        explicit Node( int m = inf, int a = -1 ) : min( m ), assign( a ) {}
    };

    int nodes_count;
    int vertices_count;
    int n;
    std::vector< Node > nodes;

    static int get_min_degree_of_two( int n);
    void proselytize(int index, int left, int right );
    int min_( int left, int right, int index, int left_border, int right_border );
    void update( int left, int right, int index, int value, int left_border, int right_border );

public:

    Segment_tree( const int* data, int m );

    void set_segment_value( int left, int right, int value );
    int get_min_on_segment( int left, int right );

    ~Segment_tree() = default;
    Segment_tree( const Segment_tree& ) = delete;
    Segment_tree( Segment_tree&& ) = delete;
    Segment_tree& operator = ( const Segment_tree& ) = delete;
    Segment_tree& operator = ( Segment_tree&& ) = delete;
};

Segment_tree::Segment_tree( const int* data, int m )
{
    n = m;
    vertices_count = get_min_degree_of_two( n );
    nodes_count = 2 * vertices_count - 1;
    nodes.assign( nodes_count, Node() );
    for( int i = n - 1; i >= 0; --i )
        nodes[nodes_count - vertices_count + i].min = data[i];
    for( int i = nodes_count - vertices_count - 1; i >= 0; --i )
        nodes[i].min = std::min( nodes[ 2 * i + 1].min, nodes[2 * i + 2].min );
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
    if( nodes[index].assign == -1 )
        return;
    nodes[2 * index + 1].assign = nodes[2 * index + 2].assign = nodes[2 * index + 1].min = nodes[2 * index + 2].min = nodes[index].assign;
    nodes[index].assign = -1;
}

int Segment_tree::min_( int left, int right, int index, int left_border, int right_border )
{
    proselytize( index, left_border, right_border );
    if( left_border > right || left > right_border )
        return inf;
    if( left_border >= left && right_border <= right )
        return nodes[index].min;
    int separation_index = ( left_border + right_border ) / 2;
    int left_sub_segment = min_( left, right, 2 * index + 1, left_border, separation_index );
    int right_sub_segment = min_( left, right, 2 * index + 2, separation_index + 1, right_border);
    return std::min( left_sub_segment, right_sub_segment );
}

void Segment_tree::update( int left, int right, int index, int value, int left_border, int right_border )
{
    proselytize( index, left_border, right_border );
    if( left_border > right || left > right_border)
        return;
    if( left_border >= left && right_border <= right )
    {
        nodes[index].min = nodes[index].assign = value;
        return;
    }
    int separation_index = ( left_border + right_border ) / 2;
    update( left, right, 2 * index + 1, value, left_border, separation_index );
    update( left, right, 2 * index + 2, value, separation_index + 1, right_border );
    nodes[index].min = std::min( nodes[2 * index + 1].min, nodes[2 * index + 2].min );
}

int Segment_tree::get_min_on_segment( int left, int right )
{
    return min_( left, right, 0, 0, vertices_count - 1 );
}

void Segment_tree::set_segment_value(int left, int right, int value)
{
    update( left, right, 0, value, 0, vertices_count - 1);
}

void process_requests( int n, int* numbers )
{
    Segment_tree segment_tree( numbers, n );
    int k, r, g, b, change_request, change_left, change_right, min_left, min_right;
    std::cin >> k;
    std::vector< int > answers;
    for( int i = 0; i < k; ++i )
    {
        std::cin >> change_left >> change_right >> r >> g >> b >> min_left >> min_right;
        change_request = r + g + b;
        segment_tree.set_segment_value( change_left, change_right, change_request );
        answers.push_back( segment_tree.get_min_on_segment( min_left, min_right ) );
    }
    for( int i = 0; i < answers.size(); ++i )
        std::cout << answers[i] << " ";
}

int main()
{
    int n;
    std::cin >> n;
    int* numbers = new int[n];
    int r, g, b;
    for( int i = 0; i < n; ++i )
    {
        std::cin >> r >> g >> b;
        numbers[i] = r + g + b;
    }
    process_requests( n, numbers );
    delete[] numbers;
    return 0;
}