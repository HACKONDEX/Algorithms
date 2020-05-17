#include <iostream>
#include <vector>
#include <cmath>

void count_logarithms_degrees_of_two( std::vector< int >& logarithms, std::vector< int >& degrees_of_two )
{
    int log = 0, value = 2;
    for( int i = 2; i < logarithms.size(); ++i )
    {
        if( value == i )
            log += 1;
        else if( value < i )
        {
            degrees_of_two.push_back( value );
            value *= 2;
        }
        logarithms[i] = log;
    }
    if( degrees_of_two.back() != value )
        degrees_of_two.push_back( value );
}

void create_sparse_table( int** sparse_table, int n, int max_log, const long long* const numbers )
{
    for( int i = 0; i < n; ++i )
        sparse_table[0][i] = i;
    for( int log = 1, range = 1; log <= max_log; ++log )
    {
        for( int i = 0; i < n; ++i )
        {
            if( i + range < n )
                sparse_table[log][i] = ( numbers[sparse_table[log - 1][i]] <= numbers[sparse_table[log - 1][i + range]] ? sparse_table[log - 1][i] : sparse_table[log - 1][i + range] );
            else
                sparse_table[log][i] = sparse_table[log - 1][i];
        }
        range *= 2;
    }
}

int range_minimum_query(int** const sparse_table, const std::vector< int >& logarithms, const std::vector< int >& degrees_of_two, const long long* const numbers, int left, int right )
{
    int first = sparse_table[logarithms[right - left + 1]][left], second = sparse_table[logarithms[right - left + 1]][right - degrees_of_two[logarithms[right - left + 1]] + 1];
    return  ( numbers[first] <= numbers[second] ) ? first : second;
}

void process_requests( int n, int m, long long* numbers, int* left, int* right )
{
    std::vector< int > logarithms(n + 5, 0), degrees_of_two;
    degrees_of_two.push_back(1);
    count_logarithms_degrees_of_two(logarithms, degrees_of_two);
    auto** sparse_table = new int*[logarithms[n - 1] + 1];
    for( int i = 0; i <= logarithms[n - 1]; ++i )
        sparse_table[i] = new int[n];
    create_sparse_table( sparse_table, n, logarithms[n - 1], numbers );
    int first_min;
    long long first, second;
    for( int i = 0; i < m; ++i )
    {
        if (right[i] == left[i])
        {
            std::cout << numbers[right[i]] << std::endl;
            continue;
        }
        first = second = 1000000009;
        first_min = range_minimum_query(sparse_table, logarithms, degrees_of_two, numbers, left[i], right[i]);
        if (first_min + 1 <= right[i])
            first = numbers[range_minimum_query(sparse_table, logarithms, degrees_of_two, numbers, first_min + 1, right[i])];
        if (first_min - 1 >= left[i])
            second = numbers[range_minimum_query(sparse_table, logarithms, degrees_of_two, numbers, left[i],first_min - 1)];
        std::cout << std::min(first, second) << std::endl;
    }
    for( int i = 0; i <= logarithms[n - 1]; ++i )
        delete[] sparse_table[i];
    delete[] sparse_table;
}

int main()
{
    int n, m;
    std::cin >> n >> m;
    auto* numbers = new long long[n];
    int* left = new int[m];
    int* right = new int[m];
    for( int i = 0; i < n; ++i )
        std::cin >> numbers[i];
    for( int j = 0; j < m; j++ )
    {
        std::cin >> left[j] >> right[j];
        --left[j];
        --right[j];
    }

    process_requests( n, m, numbers, left, right );

    delete[] right;
    delete[] left;
    delete[] numbers;
    return 0;
}