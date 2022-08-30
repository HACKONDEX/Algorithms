#include <iostream>
#include <vector>

void CountLogarithmsAndDegreesOfTwo(std::vector<int> &logarithms, std::vector<int> &degrees_of_two) {
    int log = 0, value = 2;
    for (int i = 2; i < logarithms.size(); ++i) {
        if (value == i) {
            log += 1;
        } else if (value < i) {
            degrees_of_two.push_back(value);
            value *= 2;
        }
        logarithms[i] = log;
    }
    if (degrees_of_two.back() != value) {
        degrees_of_two.push_back(value);
    }
}

void CreateSparseTable(std::vector<std::vector<int>> &sparse_table, int n, int max_log,
                       const std::vector<int64_t> &numbers) {
    for (int i = 0; i < n; ++i) {
        sparse_table[0][i] = i;
    }
    for (int log = 1, range = 1; log <= max_log; ++log) {
        for (int i = 0; i < n; ++i) {
            if (i + range < n) {
                sparse_table[log][i] = (numbers[sparse_table[log - 1][i]] <= numbers[sparse_table[log - 1][i + range]]
                                        ? sparse_table[log - 1][i] : sparse_table[log - 1][i + range]);
            } else {
                sparse_table[log][i] = sparse_table[log - 1][i];
            }
        }
        range *= 2;
    }
}

int RangeMinimumQuery(const std::vector<std::vector<int>> &sparse_table, const std::vector<int> &logarithms,
                      const std::vector<int> &degrees_of_two, const std::vector<int64_t> &numbers, int left,
                      int right) {
    int first = sparse_table[logarithms[right - left + 1]][left];
    int second = sparse_table[logarithms[right - left + 1]][right - degrees_of_two[logarithms[right - left + 1]] + 1];
    return (numbers[first] <= numbers[second]) ? first : second;
}

void ProcessRequests(int n, int m, std::vector<int64_t> &numbers, const std::vector<int> &left,
                     const std::vector<int> &right) {
    std::vector<int> logarithms(n + 5, 0);
    std::vector<int> degrees_of_two;
    degrees_of_two.push_back(1);
    CountLogarithmsAndDegreesOfTwo(logarithms, degrees_of_two);
    std::vector<std::vector<int>> sparse_table;
    for (int i = 0; i <= logarithms[n - 1]; ++i) {
        sparse_table.emplace_back(std::vector<int>(n));
    }
    CreateSparseTable(sparse_table, n, logarithms[n - 1], numbers);
    int first_min = 0;
    int64_t first;
    int64_t second;
    for (int i = 0; i < m; ++i) {
        if (right[i] == left[i]) {
            std::cout << numbers[right[i]] << std::endl;
            continue;
        }
        first = second = 1000000009;
        first_min = RangeMinimumQuery(sparse_table, logarithms, degrees_of_two, numbers, left[i], right[i]);
        if (first_min + 1 <= right[i]) {
            first = numbers[RangeMinimumQuery(sparse_table, logarithms, degrees_of_two, numbers, first_min + 1,
                                              right[i])];
        }
        if (first_min - 1 >= left[i]) {
            second = numbers[RangeMinimumQuery(sparse_table, logarithms, degrees_of_two, numbers, left[i],
                                               first_min - 1)];
        }
        std::cout << std::min(first, second) << std::endl;
    }
}

int main() {
    int n, m;
    std::cin >> n >> m;
    std::vector<int64_t> numbers(n);
    std::vector<int> left(m);
    std::vector<int> right(m);
    for (int i = 0; i < n; ++i) {
        std::cin >> numbers[i];
    }
    for (int j = 0; j < m; j++) {
        std::cin >> left[j] >> right[j];
        --left[j];
        --right[j];
    }

    ProcessRequests(n, m, numbers, left, right);
    return 0;
}