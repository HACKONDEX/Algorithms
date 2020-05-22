#include <iostream>
#include <string>


template < class T >
class Cartesian_tree
{
private:
    struct Node
    {
        T data;
        int children_count;
        long long priority;
        Node* left;
        Node* right;

        Node( const T& string, long long prior, int count = 1, Node* left_ = nullptr, Node* right_ = nullptr )
            : data( string ),priority(prior), children_count( count ), left( left_), right( right_ ) {}
    };

    Node* root;

    void delete_node( Node* leaf );
    static void update_children( Node* leaf );
    void split( Node*& tree, int count, Node*& left_tree, Node*& right_tree );
    void merge( Node*& return_tree, Node* left_tree, Node* right_tree );
    void insert_(int position, Node* node );
    void get_the_most_right( Node* node, Node*& return_node );
    void get_( int position, Node* node, Node*& return_node );


public:

    void insert( int position, const T& element );
    void remove( int left_position, int right_position );
    T& get( int position );

    Cartesian_tree() : root( nullptr ) {}
    ~Cartesian_tree();

    Cartesian_tree( const Cartesian_tree& ) = delete;
    Cartesian_tree( Cartesian_tree&& ) = delete;
    Cartesian_tree& operator = ( const Cartesian_tree& other ) = delete;
    Cartesian_tree& operator = ( Cartesian_tree&& other ) = delete;
};

template < class T >
void Cartesian_tree< T >::delete_node( Cartesian_tree::Node* leaf )
{
    if( leaf->left != nullptr )
        delete_node( leaf->left );
    if( leaf->right != nullptr )
        delete_node( leaf->right );
    delete leaf;
}

template < class T >
Cartesian_tree< T >::~Cartesian_tree()
{
    delete_node( root );
}

template < class T >
void Cartesian_tree< T >::update_children( Cartesian_tree::Node* leaf )
{
    if( leaf == nullptr )
        return;
    leaf->children_count = 1;
    if( leaf->left != nullptr )
        leaf->children_count += leaf->left->children_count;
    if( leaf->right != nullptr )
        leaf->children_count += leaf->right->children_count;
}

template < class T >
void Cartesian_tree< T >::split(Cartesian_tree::Node *&tree, int count, Cartesian_tree::Node *&left_tree, Cartesian_tree::Node *&right_tree)
{
    if( count == 0 )
    {
        right_tree = tree;
        left_tree = nullptr;
    }
    else if( count == tree->children_count )
    {
        left_tree = tree;
        right_tree = nullptr;
    }
    else
    {
        int left_count;
        if (tree->left == nullptr)
            left_count = 0;
        else
            left_count = tree->left->children_count;
        if (left_count >= count)
        {
            right_tree = tree;
            split(tree->left, count, left_tree, right_tree->left);
        }
        else
        {
            left_tree = tree;
            split(tree->right, count - left_count - 1, left_tree->right, right_tree);
        }
    }
    update_children( left_tree );
    update_children( right_tree );
}

template < class T >
void Cartesian_tree< T >::merge( Cartesian_tree::Node*& return_tree, Cartesian_tree::Node *left_tree, Cartesian_tree::Node *right_tree )
{
    if( left_tree == nullptr )
        return_tree = right_tree;
    else if( right_tree == nullptr )
        return_tree = left_tree;
    else if( left_tree->priority <= right_tree->priority )
    {
        merge( right_tree->left, left_tree, right_tree->left );
        return_tree = right_tree;
    }
    else
    {
        merge( left_tree->right, left_tree->right, right_tree );
        return_tree = left_tree;
    }
    update_children( return_tree );
}

template < class T >
void Cartesian_tree< T >::insert_( int position, Cartesian_tree::Node* node )
{
    Node* left_tree = nullptr;
    Node* right_tree = nullptr;
    split( root, position, left_tree, right_tree );
    Node* left_under_tree = nullptr;
    merge( left_under_tree, left_tree, node );
    merge( root, left_under_tree, right_tree );
}

template < class T >
void Cartesian_tree< T >::insert( int position, const T& element )
{
    Node* node = new Node( element, rand() * rand() + rand() );
    insert_( position, node );
}

template < class T >
void Cartesian_tree< T >::get_(int position, Cartesian_tree::Node *node, Cartesian_tree::Node *&return_node)
{
    if( position == node->children_count )
        get_the_most_right( node, return_node );
    else
    {
        int left_count;
        if (node->left != nullptr)
            left_count = node->left->children_count;
        else
            left_count = 0;
        if( left_count >= position )
            get_( position, node->left, return_node );
        else if( left_count + 1 == position )
            return_node = node;
        else
            get_( position - left_count - 1, node->right, return_node );
    }
}

template < class T >
T& Cartesian_tree< T >::get( int position )
{
    Node* interesting_node = nullptr;
    get_( position + 1, root, interesting_node );
    return interesting_node->data;
}

template < class T >
void Cartesian_tree< T >::get_the_most_right( Cartesian_tree::Node *node, Cartesian_tree::Node *&return_node )
{
    if( node->right == nullptr )
        return_node = node;
    else
        get_the_most_right( node->right, return_node );
}

template < class T >
void Cartesian_tree< T >::remove( int left_position, int right_position )
{
    if( left_position > right_position )
        std::swap( left_position, right_position );

    Node* left_tree_to_merge = nullptr;
    Node* right_tree = nullptr;
    split( root, left_position, left_tree_to_merge, right_tree );
    Node* to_be_deleted = nullptr;
    Node* right_tree_to_merge = nullptr;
    split(right_tree, right_position - left_position + 1, to_be_deleted, right_tree_to_merge);
    merge( root, left_tree_to_merge, right_tree_to_merge );
    delete_node( to_be_deleted );
}

//// Process requests

void process_requests( int requests_count )
{
    Cartesian_tree< std::string > tree;
    for( int i = 0; i < requests_count; ++i )
    {
        char command;
        std::cin >> command;
        if( command == '+' )
        {
            int position;
            std::cin >> position;
            std::string string;
            std::cin >> string;
            tree.insert( position, string );
        }
        else if( command == '-' )
        {
            int left, right;
            std::cin >> left >> right;
            tree.remove( left, right );
        }
        else
        {
            int position;
            std::cin >> position;
            std::cout << tree.get( position ) << std::endl;
        }
    }
}

int main()
{
    int requests_count;
    std::cin >> requests_count;
    process_requests( requests_count );
    return 0;
}