#include <iostream>
#include <string>
#include <vector>

void CalculatePalindromesLengthForEachCentre(std::vector<size_t> &length, const std::string &text) {
    size_t symmetric_pos, border_centre = 0; // centre of the palindrome with rightest border
    for (size_t i = 1; i < 2 * text.length(); ++i) {
        border_centre = length[i - 1] + i - 1 > border_centre + length[border_centre] ? i - 1 : border_centre;
        if (i <= border_centre + length[border_centre]) {
            symmetric_pos = 2 * border_centre - i;
            if (symmetric_pos - length[symmetric_pos] > border_centre - length[border_centre]) {
                length[i] = length[symmetric_pos];
                continue;
            } else {
                length[i] = symmetric_pos - (border_centre - length[border_centre]);
            }
        }
        for (int l = i - length[i] - 1, j = i + length[i] + 1;
             l >= 0 && j <= 2 * text.length() - 1; --l, ++j, ++length[i]) {
            if (!(j % 2 == 0 || (j % 2 == 1 && text[(l - 1) / 2] == text[(j - 1) / 2]))) {
                break;
            }
        }
    }
}

unsigned long long FindSubstringPalindromesCount(const std::string &text) {
    std::vector<size_t> palindrome_lengths(text.length() * 2, 0);
    CalculatePalindromesLengthForEachCentre(palindrome_lengths, text);
    unsigned long long palindrome_count = 0;
    for (size_t i = 0; i < palindrome_lengths.size(); ++i) {
        palindrome_count += i % 2 == 0 ? (palindrome_lengths[i] + 1) / 2 : (palindrome_lengths[i] / 2);
    }
    return palindrome_count;
}

int main() {
    std::string input_text;
    std::cin >> input_text;

    std::cout << FindSubstringPalindromesCount(input_text) << std::endl;
    return 0;
}