#include <array>
#include <iostream>
#include <memory>
#include <queue>
#include <string>
#include <vector>

using std::shared_ptr;
using std::weak_ptr;

const int ALPHABET_SIZE = 26;
const int MINIMAL_SYMBOL_NUMBER = 97;

class Trie {
  private:
    struct Node {
        unsigned char
                vertex_status_;  /// 0 - root, 1 - not terminal, 2 - terminal
        weak_ptr<Node> parent_;
        std::vector<int> word_index_;
        weak_ptr<Node> suffix_link_ = shared_ptr<Node>(nullptr);
        weak_ptr<Node> compressed_suffix_link_ = shared_ptr<Node>(nullptr);
        std::array<shared_ptr<Node>, ALPHABET_SIZE> next_vertices_ = {nullptr};
        std::array<weak_ptr<Node>, ALPHABET_SIZE> transition_suffix_link_ = {
                shared_ptr<Node>(nullptr)};

        Node(unsigned char status,
             const weak_ptr<Node> &parent = shared_ptr<Node>(nullptr))
                : vertex_status_(status), parent_(parent) {}
    };

  public:
    Trie() {
        root_ = std::make_shared<Node>(0);
        ConstructTrie();
    }
    Trie(const std::vector<std::string> &words) {
        words_set_.reserve(words.size());
        for (const auto &str: words) {
            words_set_.push_back(str);
        }

        root_ = std::make_shared<Node>(0);
        ConstructTrie();
    }

    void AddString(const std::string &str);
    void FindAllEntriesInText(std::vector<int> &words_count_in_each_text_position,
                              const std::vector<int> &word_position_in_template,
                              const std::string &text);

  private:
    void AddStringPrivate(const std::string_view &string,
                          int index_in_words_vector);
    void ConstructTrie();
    void ConstructDirectAndCompressedSuffixLinks();
    void FindSuffixLink(
            weak_ptr<Node> &current_node, weak_ptr<Node> &current_link,
            const std::queue<std::pair<size_t, weak_ptr<Node>>> &queue);
    void GoSuffixLink(weak_ptr<Node> &answer, const weak_ptr<Node> &suff_link,
                      int symbol);

  private:
    std::shared_ptr<Node> root_;
    std::vector<std::string_view> words_set_;

};

void Trie::ConstructTrie() {
    root_->parent_ = root_;
    for (size_t i = 0; i < words_set_.size(); ++i) {
        AddStringPrivate(words_set_[i], i);
    }
    ConstructDirectAndCompressedSuffixLinks();
}

void Trie::FindSuffixLink(
        weak_ptr<Node> &current_node, weak_ptr<Node> &current_link,
        const std::queue<std::pair<size_t, weak_ptr<Node>>> &queue) {
    while (true) {
        if (current_link.lock()->vertex_status_ == 0) {
            current_node.lock()->suffix_link_ =
                    root_->next_vertices_[queue.front().first] != nullptr &&
                    current_node.lock()->parent_.lock() != root_
                    ? root_->next_vertices_[queue.front().first]
                    : root_;
            break;
        }
        if (current_link.lock()->next_vertices_[queue.front().first] == nullptr) {
            current_link = current_link.lock()->suffix_link_;
        } else {
            current_node.lock()->suffix_link_ =
                    current_link.lock()->next_vertices_[queue.front().first];
            break;
        }
    }
}

void Trie::ConstructDirectAndCompressedSuffixLinks() {
    std::queue<std::pair<size_t, weak_ptr<Node>>> queue;
    weak_ptr<Node> current_node;
    weak_ptr<Node> current_link;
    root_->suffix_link_ = root_;
    for (size_t i = 0; i < root_->next_vertices_.size(); ++i) {
        if (root_->next_vertices_[i] != nullptr) {
            queue.push(std::make_pair(i, root_->next_vertices_[i]));
        }
    }
    while (!queue.empty()) {
        current_node = queue.front().second;
        for (size_t i = 0; i < current_node.lock()->next_vertices_.size(); ++i) {
            if (current_node.lock()->next_vertices_[i] != nullptr) {
                queue.push(std::make_pair(i, current_node.lock()->next_vertices_[i]));
            }
        }
        current_link = current_node.lock()->parent_.lock()->suffix_link_;
        FindSuffixLink(current_node, current_link, queue);
        current_node.lock()->compressed_suffix_link_ =
                current_node.lock()->suffix_link_.lock()->vertex_status_ == 2
                ? current_node.lock()->suffix_link_
                : current_node.lock()->suffix_link_.lock()->compressed_suffix_link_;
        queue.pop();
    }
}

void Trie::AddStringPrivate(const std::string_view &string,
                            int index_in_words_vector) {
    shared_ptr<Node> current_node = root_;
    for (auto symbol: string) {
        if (current_node->next_vertices_[symbol - MINIMAL_SYMBOL_NUMBER].get() ==
            nullptr) {
            current_node->next_vertices_[symbol - MINIMAL_SYMBOL_NUMBER] =
                    std::make_shared<Node>(1, current_node);
        }
        current_node = current_node->next_vertices_[symbol - MINIMAL_SYMBOL_NUMBER];
    }
    current_node->vertex_status_ = 2;
    current_node->word_index_.push_back(index_in_words_vector);
}

void Trie::AddString(const std::string &str) {
    words_set_.push_back(str);
    AddStringPrivate(str, words_set_.size() - 1);
}

void Trie::GoSuffixLink(weak_ptr<Node> &answer, const weak_ptr<Node> &suff_link,
                        int symbol) {
    if (suff_link.lock()->next_vertices_[symbol] != nullptr) {
        answer = suff_link.lock()->next_vertices_[symbol];
    } else if (suff_link.lock()->transition_suffix_link_[symbol].lock() !=
               nullptr) {
        answer = suff_link.lock()->transition_suffix_link_[symbol];
    } else if (suff_link.lock()->vertex_status_ == 0) {
        answer = shared_ptr<Node>(nullptr);
    } else {
        GoSuffixLink(suff_link.lock()->transition_suffix_link_[symbol],
                     suff_link.lock()->suffix_link_, symbol);
        answer = suff_link.lock()->transition_suffix_link_[symbol];
    }
}

void Trie::FindAllEntriesInText(
        std::vector<int> &words_count_in_each_text_position,
        const std::vector<int> &word_position_in_template,
        const std::string &text) {
    weak_ptr<Node> current_node = root_;
    weak_ptr<Node> suff_link;
    weak_ptr<Node> compressed_link;
    int symbol;
    int possible_template_start;
    for (int i = 0; i < text.length(); ++i) {
        symbol = text[i] - MINIMAL_SYMBOL_NUMBER;
        /// find next vertex
        if (current_node.lock()->next_vertices_[symbol] != nullptr) {
            current_node = current_node.lock()->next_vertices_[symbol];
        } else if (current_node.lock()->transition_suffix_link_[symbol].lock() !=
                   nullptr) {
            current_node = current_node.lock()->transition_suffix_link_[symbol];
        } else {
            suff_link = current_node.lock()->suffix_link_;
            GoSuffixLink(current_node.lock()->transition_suffix_link_[symbol],
                         suff_link, symbol);
            current_node =
                    current_node.lock()->transition_suffix_link_[symbol].lock() != nullptr
                    ? current_node.lock()->transition_suffix_link_[symbol]
                    : root_;
        }

        /// terminal
        compressed_link = current_node;
        while (true) {
            if (compressed_link.lock() == nullptr) {
                break;
            }
            for (auto word_index: compressed_link.lock()->word_index_) {
                possible_template_start = i - words_set_[word_index].length() + 1 -
                                          word_position_in_template[word_index];
                if (possible_template_start >= 0) {
                    ++words_count_in_each_text_position[possible_template_start];
                }
            }
            compressed_link = compressed_link.lock()->compressed_suffix_link_;
        }
    }
}

/// Task Solution

void GetStringBlocksPositions(std::vector<int> &blocks_start_positions,
                              std::vector<std::string> &words_list,
                              const std::string &template_str) {
    int index = 0;
    std::string block;
    while (index < template_str.length()) {
        block.clear();
        while (index < template_str.length() && template_str[index] != '?') {
            block += template_str[index];
            ++index;
        }
        if (!block.empty()) {
            words_list.emplace_back(block);
            blocks_start_positions.push_back(index - block.length());
        }
        while (index < template_str.length() && template_str[index] == '?') {
            ++index;
        }
    }
}

void RecoverEntryPositions(
        std::vector<int> &entries_positions,
        std::vector<int> &blocks_count_in_each_text_position,
        std::pair<int, int> &&template_text_length, int blocks_count) {
    if (blocks_count == 0) {
        for (int i = 0;
             i < template_text_length.second - template_text_length.first + 1;
             ++i) {
            entries_positions.push_back(i);
        }
        return;
    }
    for (int i = 0; i < blocks_count_in_each_text_position.size(); ++i) {
        if (blocks_count_in_each_text_position[i] == blocks_count &&
            i + template_text_length.first - 1 < template_text_length.second) {
            entries_positions.push_back(i);
        }
    }
}

void FindAllEntriesOfTemplateInText(std::vector<int> &entries_positions,
                                    const std::string &template_str,
                                    const std::string &input_text) {
    std::vector<int> blocks_start_positions;
    std::vector<std::string> word_list;
    GetStringBlocksPositions(blocks_start_positions, word_list, template_str);
    std::vector<int> blocks_count_in_each_text_position(input_text.length(), 0);
    Trie bor(word_list);
    bor.FindAllEntriesInText(blocks_count_in_each_text_position,
                             blocks_start_positions, input_text);
    RecoverEntryPositions(
            entries_positions, blocks_count_in_each_text_position,
            std::make_pair(template_str.length(), input_text.length()),
            word_list.size());
}

int main() {
    std::string template_string;
    std::string input_text;
    std::cin >> template_string;
    std::cin >> input_text;

    std::vector<int> entries;
    FindAllEntriesOfTemplateInText(entries, template_string, input_text);
    for (int i: entries) {
        std::cout << i << " ";
    }
    std::cout << '\n';
    return 0;
}
