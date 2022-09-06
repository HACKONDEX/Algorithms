#include <iostream>
#include <string>
#include <vector>

void CalculatePrefixValues(const std::string &string, std::vector<int> &prefix_values) {
    int current_substr_length = 0;
    prefix_values[0] = -1;
    for (int i = 2; i <= string.length(); ++i) {
        current_substr_length = prefix_values[i - 1];
        while (current_substr_length >= 0) {
            if (string[i - 1] == string[current_substr_length]) {
                prefix_values[i] = current_substr_length + 1;
                break;
            }
            current_substr_length = prefix_values[current_substr_length];
        }
    }
}

void FindEntryPositions(const std::string &input_string, const std::string &search_template,
                        std::vector<int> &entry_positions) {
    int template_length = search_template.length();
    std::vector<int> prefix_values(template_length + 1, 0);
    CalculatePrefixValues(search_template, prefix_values);
    int previous_prefix_length = 0;
    int current_prefix_length = 0;
    bool suffix_found = false;
    for (int i = 0; i < input_string.length(); ++i) {
        previous_prefix_length = suffix_found ? current_prefix_length : 0;
        suffix_found = false;
        while (previous_prefix_length >= 0) {
            if (input_string[i] == search_template[previous_prefix_length]) {
                suffix_found = true;
                current_prefix_length = previous_prefix_length + 1;
                if (current_prefix_length == template_length)
                    entry_positions.push_back(i - template_length + 1);
                break;
            }
            previous_prefix_length = prefix_values[previous_prefix_length];
        }
    }
}

int main() {
    std::string search_template;
    std::string input_string;

    std::cin >> search_template;
    std::cin >> input_string;

    std::vector<int> entry_positions;
    FindEntryPositions(input_string, search_template, entry_positions);

    for (int i: entry_positions) {
        std::cout << i << " ";
    }
    std::cout << '\n';

    return 0;
}