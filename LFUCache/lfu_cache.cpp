#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template <class U, class V>
using HashMap = std::unordered_map<U, V>;

struct ListNode;

using Node = ListNode*;

struct ListNode {
    Node next;
    Node prev;
    int key;
    int value;

    ListNode(int key, int value)
        : next(nullptr), prev(nullptr), key(key), value(value) {}
};

class DoubleLinkedList {
  public:
    DoubleLinkedList() = default;

    Node PopFront();
    Node PushBack(Node node);
    Node Erase(Node node);

    int Size() const;

  private:
    Node head{nullptr};
    Node tail{nullptr};
    int size = 0;
};

class LFUCache {
  public:
    explicit LFUCache(int capacity)
        : capacity(capacity), size(0), min_frequency(1) {
        freq_map_list.insert({1, DoubleLinkedList()});
    }

    int get(int key);
    void put(int key, int value);

  private:
    void UpdateFreq(int key);

  private:
    int capacity;
    int size;
    int min_frequency;
    HashMap<int, DoubleLinkedList> freq_map_list;
    HashMap<int, Node> key_map_node;
    HashMap<int, int> key_map_freq;
};

////////// List //////////

int DoubleLinkedList::Size() const { return size; }

Node DoubleLinkedList::PopFront() {
    if (size == 0) {
        std::cerr << "PopFront from empty list\n";
        std::abort();
    }
    if (size == 1) {
        --size;
        auto old_head = head;
        head = tail = nullptr;
        return old_head;
    }
    Node new_head = head->next;
    head->next = nullptr;
    new_head->prev = nullptr;
    Node old_head = head;
    head = new_head;
    --size;
    return old_head;
}

Node DoubleLinkedList::PushBack(Node node) {
    if (size == 0) {
        head = tail = node;
        ++size;
        return node;
    }
    tail->next = node;
    node->prev = tail;
    tail = node;
    ++size;
    return node;
}

Node DoubleLinkedList::Erase(Node node) {
    if (node->next == nullptr) {
        if (node != tail) {
            std::cerr << "Erase from wrong List, tail doesn't match\n";
            std::abort();
        }
        if (size == 1) {
            head = tail = nullptr;
            size = 0;
            return node;
        }
        tail = tail->prev;
        tail->next = nullptr;
        node->prev = nullptr;
        --size;
        return node;
    } else if (node->prev == nullptr) {
        if (node != head) {
            std::cerr << "Erase from wrong List, head doesn't match\n";
            std::abort();
        }
        if (size == 1) {
            head = tail = nullptr;
            size = 0;
            return node;
        }
        head = head->next;
        head->prev = nullptr;
        node->next = nullptr;
        --size;
        return node;
    }
    --size;
    auto next = node->next;
    auto prev = node->prev;
    node->next = nullptr;
    node->prev = nullptr;
    prev->next = next;
    next->prev = prev;
    return node;
}

////////// LFUCache //////////

int LFUCache::get(int key) {
    auto iter = key_map_node.find(key);
    if (iter == key_map_node.end()) {
        return -1;
    }

    int value = iter->second->value;
    UpdateFreq(key);
    return value;
}

void LFUCache::put(int key, int value) {
    if (capacity == 0) {
        return;
    }

    {
        auto node = key_map_node.find(key);
        if (node != key_map_node.end()) {
            node->second->value = value;
            UpdateFreq(key);
            return;
        }
    }

    if (size == capacity) {
        /// Remove LeastFrequentNode
        Node node_to_delete = freq_map_list[min_frequency].PopFront();
        key_map_node.erase(node_to_delete->key);
        key_map_freq.erase(node_to_delete->key);
        delete node_to_delete;
        --size;
    }

    DoubleLinkedList& list = freq_map_list[1];
    auto node = new ListNode(key, value);
    list.PushBack(node);
    key_map_node.insert({key, node});
    key_map_freq.insert({key, 1});
    ++size;
    min_frequency = 1;
}

void LFUCache::UpdateFreq(int key) {
    Node prev_node = key_map_node[key];
    int prev_freq = key_map_freq[key];

    DoubleLinkedList& prev_freq_list = freq_map_list[prev_freq];
    prev_freq_list.Erase(prev_node);

    int new_freq = prev_freq + 1;
    if (freq_map_list.find(new_freq) == freq_map_list.end()) {
        freq_map_list.insert({new_freq, DoubleLinkedList()});
    }
    DoubleLinkedList& new_freq_list = freq_map_list[new_freq];
    auto node = new_freq_list.PushBack(prev_node);

    key_map_node[key] = node;
    key_map_freq[key] = new_freq;
    if (prev_freq_list.Size() == 0 && prev_freq == min_frequency) {
        ++min_frequency;
    }
}

int main() {
    {
        std::cout << "Test 1\n";
        LFUCache c(2);
        c.put(1, 1);
        c.put(2, 2);
        std::cout << c.get(1) << endl;
        c.put(3, 3);
        std::cout << c.get(2) << endl;
        std::cout << c.get(3) << endl;
        c.put(4, 4);
        std::cout << c.get(1) << endl;
        std::cout << c.get(3) << endl;
        std::cout << c.get(4) << endl;
    }

    {
        std::cout << "Test 2\n";
        LFUCache c(3);
        c.put(2, 2);
        c.put(1, 1);
        std::cout << c.get(2) << endl;
        std::cout << c.get(1) << endl;
        std::cout << c.get(2) << endl;
        c.put(3, 3);
        c.put(4, 4);
        std::cout << c.get(3) << endl;
        std::cout << c.get(2) << endl;
        std::cout << c.get(1) << endl;
        std::cout << c.get(4) << endl;
    }
    return 0;
}
