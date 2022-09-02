#include <iostream>
#include <string>

template<class T>
class CartesianTree {
  private:
    struct Node {
        T data;
        int children_count;
        long long priority;
        Node *left;
        Node *right;

        Node(const T &string, long long prior, int count = 1, Node *left_ = nullptr, Node *right_ = nullptr)
                : data(string), priority(prior), children_count(count), left(left_), right(right_) {}
    };

  public:
    CartesianTree() : root(nullptr) {}
    ~CartesianTree();
    CartesianTree(const CartesianTree &) = default;
    CartesianTree(CartesianTree &&) = default;

    CartesianTree &operator=(const CartesianTree &other) = delete;
    CartesianTree &operator=(CartesianTree &&other) = delete;

    void InsertElement(int position, const T &element);
    void Remove(int left_position, int right_position);
    T &Get(int position);

  private:
    static void UpdateChildren(Node *leaf);
    void DeleteNode(Node *leaf);
    void Split(Node *&tree, int count, Node *&left_tree, Node *&right_tree);
    void Merge(Node *&return_tree, Node *left_tree, Node *right_tree);
    void Insert(int position, Node *node);
    void GetMostRight(Node *node, Node *&return_node);
    void Get(int position, Node *node, Node *&return_node);

  private:
    Node *root;

};

template<class T>
void CartesianTree<T>::DeleteNode(CartesianTree::Node *leaf) {
    if (leaf->left != nullptr) {
        DeleteNode(leaf->left);
    }
    if (leaf->right != nullptr) {
        DeleteNode(leaf->right);
    }
    delete leaf;
}

template<class T>
CartesianTree<T>::~CartesianTree() {
    DeleteNode(root);
}

template<class T>
void CartesianTree<T>::UpdateChildren(CartesianTree::Node *leaf) {
    if (leaf == nullptr) {
        return;
    }
    leaf->children_count = 1;
    if (leaf->left != nullptr) {
        leaf->children_count += leaf->left->children_count;
    }
    if (leaf->right != nullptr) {
        leaf->children_count += leaf->right->children_count;
    }
}

template<class T>
void CartesianTree<T>::Split(CartesianTree::Node *&tree, int count, CartesianTree::Node *&left_tree,
                             CartesianTree::Node *&right_tree) {
    if (count == 0) {
        right_tree = tree;
        left_tree = nullptr;
    } else if (count == tree->children_count) {
        left_tree = tree;
        right_tree = nullptr;
    } else {
        int left_count;
        if (tree->left == nullptr) {
            left_count = 0;
        } else {
            left_count = tree->left->children_count;
        }
        if (left_count >= count) {
            right_tree = tree;
            Split(tree->left, count, left_tree, right_tree->left);
        } else {
            left_tree = tree;
            Split(tree->right, count - left_count - 1, left_tree->right, right_tree);
        }
    }
    UpdateChildren(left_tree);
    UpdateChildren(right_tree);
}

template<class T>
void CartesianTree<T>::Merge(CartesianTree::Node *&return_tree, CartesianTree::Node *left_tree,
                             CartesianTree::Node *right_tree) {
    if (left_tree == nullptr) {
        return_tree = right_tree;
    } else if (right_tree == nullptr) {
        return_tree = left_tree;
    } else if (left_tree->priority <= right_tree->priority) {
        Merge(right_tree->left, left_tree, right_tree->left);
        return_tree = right_tree;
    } else {
        Merge(left_tree->right, left_tree->right, right_tree);
        return_tree = left_tree;
    }
    UpdateChildren(return_tree);
}

template<class T>
void CartesianTree<T>::Insert(int position, CartesianTree::Node *node) {
    Node *left_tree = nullptr;
    Node *right_tree = nullptr;
    Split(root, position, left_tree, right_tree);
    Node *left_under_tree = nullptr;
    Merge(left_under_tree, left_tree, node);
    Merge(root, left_under_tree, right_tree);
}

template<class T>
void CartesianTree<T>::InsertElement(int position, const T &element) {
    Node *node = new Node(element, rand() * rand() + rand());
    Insert(position, node);
}

template<class T>
void CartesianTree<T>::Get(int position, CartesianTree::Node *node, CartesianTree::Node *&return_node) {
    if (position == node->children_count) {
        GetMostRight(node, return_node);
    } else {
        int left_count;
        if (node->left != nullptr) {
            left_count = node->left->children_count;
        } else {
            left_count = 0;
        }
        if (left_count >= position) {
            Get(position, node->left, return_node);
        } else if (left_count + 1 == position) {
            return_node = node;
        } else {
            Get(position - left_count - 1, node->right, return_node);
        }
    }
}

template<class T>
T &CartesianTree<T>::Get(int position) {
    Node *interesting_node = nullptr;
    Get(position + 1, root, interesting_node);
    return interesting_node->data;
}

template<class T>
void CartesianTree<T>::GetMostRight(CartesianTree::Node *node, CartesianTree::Node *&return_node) {
    if (node->right == nullptr) {
        return_node = node;
    } else {
        GetMostRight(node->right, return_node);
    }
}

template<class T>
void CartesianTree<T>::Remove(int left_position, int right_position) {
    if (left_position > right_position) {
        std::swap(left_position, right_position);
    }

    Node *left_tree_to_merge = nullptr;
    Node *right_tree = nullptr;
    Split(root, left_position, left_tree_to_merge, right_tree);
    Node *to_be_deleted = nullptr;
    Node *right_tree_to_merge = nullptr;
    Split(right_tree, right_position - left_position + 1, to_be_deleted, right_tree_to_merge);
    Merge(root, left_tree_to_merge, right_tree_to_merge);
    DeleteNode(to_be_deleted);
}

void ProcessRequests(int requests_count) {
    CartesianTree<std::string> tree;
    for (int i = 0; i < requests_count; ++i) {
        char command;
        std::cin >> command;
        if (command == '+') {
            int position;
            std::cin >> position;
            std::string string;
            std::cin >> string;
            tree.InsertElement(position, string);
        } else if (command == '-') {
            int left, right;
            std::cin >> left >> right;
            tree.Remove(left, right);
        } else {
            int position;
            std::cin >> position;
            std::cout << tree.Get(position) << std::endl;
        }
    }
}

int main() {
    int requests_count;
    std::cin >> requests_count;
    ProcessRequests(requests_count);
    return 0;
}