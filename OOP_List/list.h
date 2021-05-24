#ifndef OOP_LIST_LIST_H
#define OOP_LIST_LIST_H

#include <type_traits>
#include <iostream>

template <class T> class List;

template <class T>
class List
{
private:

    struct Node;
    struct Value_node;

    struct Node
    {
        Node* next;
        Node* previous;
        Value_node* value_node;

        Node():next(nullptr), previous(nullptr), value_node(nullptr){}
    };

    struct Value_node
    {
        T value;

//        Value_node( const T& value_ ) : value( value_ ) { std::cout << "COOOOOPY" << std::endl;}
//        Value_node( T&& other ) : value( std::move( other ) ) { //std::cout << "MOOOOOOOVE" << std::endl;}
        template<class... Args>
                Value_node( Args&&... args ):value( std::forward<Args>(args)...){}
    };

public:
    template <bool is_const>
    class List_iterator;

    template <bool is_const = true>
    class List_iterator
    {
    private:
        friend class List;
        Node* node;
    public:
        //// const iterator
        typedef typename std::conditional<is_const, const T&, T&>::type Reference;
        //// non const iterator
        typedef typename std::conditional<is_const, const T*, T*>::type Pointer;

        using difference_type = std::ptrdiff_t;
        using value_type = T;
        using reference = Reference ;
        using pointer = Pointer ;
        using iterator_category = std::bidirectional_iterator_tag;

        List_iterator(Node* nd = nullptr) : node(nd) {}
        List_iterator(const List_iterator<false>& other ) : node( other.node ) {}
        List_iterator(const List_iterator<true>& other ) : node( other.node ) {}
        T& operator * () const;// dereferencing operator
        List_iterator& operator=(const List_iterator<true>& other);
        List_iterator& operator=(const List_iterator<false>& other);
        bool operator==(const List_iterator<true>& other) const;
        bool operator==(const List_iterator<false>& other) const;
        bool operator!=(const List_iterator<true>& other) const;
        bool operator!=(const List_iterator<false>& other) const;
        List_iterator& operator++();//prefix
        List_iterator& operator--();//prefix
        List_iterator operator++(int);
        List_iterator operator--(int);
    };
private:

    size_t size_;
    Node* head;
    Node* tail;
public:

    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using reference =  T&;
    using pointer =  T*;
    using iterator_category = std::bidirectional_iterator_tag;
    using iterator = List_iterator<false>;
    using const_iterator = List_iterator<true>;

    //Good story to tell
//    using Iterator = List_iterator<false>;
//    using Const_iterator = List_iterator<true>;


    List();
    List(size_t size_, const T& value = T() );

    ~List();
    List& operator = ( const List& other );
    List& operator = ( List&& other );
    List( const List& other );
    List( List&& other );

    size_t size() const;
    T& front();
    T& back();
    T front() const;
    T back() const;
    bool empty() const;
    void clear();
    void reverse();
    void unique();

    void push_back(const T& value );
    void push_back( T&& value );
    void push_front( const T& value );
    void push_front( T&& value );
    void pop_back();
    void pop_front();

    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;
    std::reverse_iterator<iterator> rbegin();
    std::reverse_iterator<iterator> rend();
    std::reverse_iterator<const_iterator> crbegin() const;
    std::reverse_iterator<const_iterator> crend() const;

    iterator insert(const_iterator position, const T& value );
    iterator insert(const_iterator position, T&& value );
    template< class InputIt >
    iterator insert( iterator pos, InputIt first, InputIt last);

    template <class... Args>
    void emplace( const_iterator position, Args&&... args);
    template <class... Args>
    void emplace_back(Args&&... args);

    template <class... Args>
    void emplace_front(Args&&... args);

    void erase( iterator pos );
    void erase( iterator first, iterator last );

private:
    void clean_body_memory();
    void copy( const List& other );
    void insert_before_position(const_iterator position, Node* new_node);
    void delete_node_in_position(const_iterator position);
};

//// Iterator functions

template <class T>
template <bool is_const>
T& List<T>::List_iterator<is_const>::operator*() const
{
    return (node->value_node->value);
}

template <class T>
template <bool is_const>
bool List<T>::List_iterator<is_const>::operator==(const List<T>::List_iterator<true> &other) const
        {
            return this->node == other.node;
        }

template <class T>
template <bool is_const>
bool List<T>::List_iterator<is_const>::operator==(const List<T>::List_iterator<false> &other) const
{
    return this->node == other.node;
}

template <class T>
template <bool is_const>
bool List<T>::List_iterator<is_const>::operator!=(const List<T>::List_iterator<true> &other) const
{
    return this->node != other.node;
}

template <class T>
template <bool is_const>
bool List<T>::List_iterator<is_const>::operator!=(const List<T>::List_iterator<false> &other) const
{
    return this->node != other.node;
}

template <class T>
template <bool is_const>
typename List<T>::template List_iterator<is_const>& List<T>::List_iterator<is_const>::operator++()
{
    if( node->next == nullptr )
        return *this;
    node = node->next;
    return *this;
}

template <class T>
template <bool is_const>
typename List<T>::template List_iterator<is_const>& List<T>::List_iterator<is_const>::operator--()
{
    node = node->previous;
    return *this;
}

template <class T>
template <bool is_const>
typename List<T>::template List_iterator<is_const> List<T>::List_iterator<is_const>::operator++(int)
{
    Node* tmp = node;
    node = node->next;
    return List_iterator<is_const>(tmp);
}

template <class T>
template <bool is_const>
typename List<T>::template List_iterator<is_const> List<T>::List_iterator<is_const>::operator--(int)
{
    Node* tmp = node;
    node = node->previous;
    return List_iterator<is_const>(tmp);
}

template <class T>
template <bool is_const>
typename List<T>::template List_iterator<is_const>& List<T>::List_iterator<is_const>::operator =(const List<T>::List_iterator<false>& other)
{
    node = other.node;
    return *this;
}

template <class T>
template <bool is_const>
typename List<T>::template List_iterator<is_const>& List<T>::List_iterator<is_const>::operator = ( const List<T>::List_iterator<true>& other )
{
    node = other.node;
    return *this;
}

//// List functions

template <class T>
List<T>::List() : size_(0)
{
    head = new Node();
    tail = new Node();
    head->next = tail;
    tail->previous = head;
}

template <class T>
List<T>::List(size_t size_, const T &value) : size_( size_ )
{
    head = new Node();
    tail = new Node();
    head->next = tail;
    tail->previous = head;
    if(size_ == 0)
        return ;
    Node* body = new Node();
    Node* prev = nullptr;
    body->value_node = new Value_node( value );
    head->next = body;
    body->previous = head;
    prev = body;
    for( size_t i = 1; i < size_; ++i )
    {
        body = new Node();
        body->value_node = new Value_node( value );
        body->previous = prev;
        prev->next = body;
        prev = body;
    }
    prev->next = tail;
    tail->previous = prev;
}

template <class T>
void List<T>::copy(const List &other)
{
    size_ = other.size_;
    if( size_ == 0 )
        return;
    Node* body = new Node();
    Node* other_body = other.head->next;
    body->value_node = new Value_node(other_body->value_node->value);
    head->next = body;
    body->previous = head;
    Node* prev = body;
    for( size_t i = 1; i < size_; ++i )
    {
        other_body = other_body->next;
        body = new Node();
        body->value_node = new Value_node(other_body->value_node->value);
        body->previous = prev;
        prev->next = body;
        prev = body;
    }
    tail->previous = prev;
    prev->next = tail;
}

template <class T>
List<T>::List( const List<T>& other )
{
    head = new Node();
    tail = new Node();
    head->next = tail;
    tail->previous = head;
    copy( other );
}

template <class T>
List<T>::List(List &&other) : size_( other.size_ )
{
    head = new Node();
    tail = new Node();
    head->next = tail;
    tail->previous = head;
    if( size_ == 0 )
        return ;
    head->next = other.head->next;
    other.head->next->previous = head;
    tail->previous = other.tail->previous;
    other.tail->previous->next = tail;
    other.head->next = other.tail;
    other.tail->previous = other.head;
    other.size_ = 0;
}

template <class T>
void List<T>::clean_body_memory()
{
    if( size_ == 0 )
        return;
    Node* body = head->next;
    Node* prev = nullptr;
    for( size_t i = 0; i < size_; ++i )
    {
        prev = body;
        body = body->next;
        delete prev->value_node;
        delete prev;
    }
    size_ = 0;
    head->next = tail;
    tail->previous = head;
}

template <class T>
List<T>::~List()
{
    clean_body_memory();
    delete head;
    delete tail;
}

template <class T>
List<T>& List<T>::operator=(const List& other)
{
    clean_body_memory();
    copy( other );
    return *this;
}

template <class T>
List<T>& List<T>::operator=( List &&other)
{
    clean_body_memory();
    size_ = other.size_;
    if( size_ == 0 )
        return *this;
    head->next = other.head->next;
    other.head->next->previous = head;
    tail->previous = other.tail->previous;
    other.tail->previous->next = tail;
    other.head->next = other.tail;
    other.tail->previous = other.head;
    other.size_ = 0;
    return *this;
}

template <class T>
size_t List<T>::size() const
{
    return size_;
}

template <class T>
T& List<T>::front()
{
    return head->next->value_node->value;
}

template <class T>
T List<T>::front() const
{
    return head->next->value_node->value;
}

template <class T>
T& List<T>::back()
{
    return tail->previous->value_node->value;
}

template <class T>
T List<T>::back() const
{
    return tail->previous->value_node->value;
}

template <class T>
bool List<T>::empty() const
{
    return size_ == 0;
}

template <class T>
void List<T>::clear()
{
    clean_body_memory();
}

template <class T>
typename List<T>::iterator List<T>::begin()
{
    return iterator( head->next );
}

template <class T>
typename List<T>::iterator List<T>::end()
{
    return iterator( tail );
}

template <class T>
typename List<T>::const_iterator List<T>::cbegin() const
{
    return const_iterator( head->next );
}

template <class T>
typename List<T>::const_iterator List<T>::cend() const
{
    return const_iterator( tail );
}

template <class T>
std::reverse_iterator<typename List<T>::iterator> List<T>::rbegin()
{
    return std::make_reverse_iterator<iterator>(iterator(tail->previous));
}

template <class T>
std::reverse_iterator<typename List<T>::iterator> List<T>::rend()
{
    return std::make_reverse_iterator<iterator>(iterator(head->next));
}

template <class T>
std::reverse_iterator<typename List<T>::const_iterator> List<T>::crbegin() const
{
    return std::make_reverse_iterator<const_iterator>(const_iterator(tail->previous));
}

template <class T>
std::reverse_iterator<typename List<T>::const_iterator> List<T>::crend() const
{
    return std::make_reverse_iterator<const_iterator>(const_iterator(head->next));
}

template <class T>
void List<T>::insert_before_position(List::const_iterator position, Node* new_node)
{
    ++size_;
    Node* prev = position.node->previous;
    prev->next = new_node;
    position.node->previous = new_node;
    new_node->next = position.node;
    new_node->previous = prev;
}

template <class T>
typename List<T>::iterator List<T>::insert(List::const_iterator position, const T& value)
{
    Node* body = new Node();
    body->value_node = new Value_node( value );
    insert_before_position(position, body);
    return position;
}

template <class T>
typename List<T>::iterator List<T>::insert(List::const_iterator position, T&& value)
{
    Node* body = new Node();
    body->value_node = new Value_node( std::move(value) );
    insert_before_position(position, body);
    return position;
}

template <class T>
template <class InputIt>
typename List<T>::iterator List<T>::insert(List::iterator pos, InputIt first, InputIt last)
{
    while ( first != last )
    {
        insert( pos, *first );
        ++first;
    }
    return pos;
}

template <class T>
void List<T>::delete_node_in_position(const_iterator position )
{
    --size_;
    Node* next_node = position.node->next;
    Node* prev_node = position.node->previous;
    next_node->previous = prev_node;
    prev_node->next = next_node;
    delete position.node->value_node;
    delete position.node;
}

template <class T>
void List<T>::erase(List::iterator pos)
{
    delete_node_in_position(pos);
}

template <class T>
void List<T>::erase(List::iterator first, List::iterator last)
{
    const_iterator tmp;
    while(first != last)
    {
        tmp = first;
        ++first;
        delete_node_in_position(tmp);
    }
}

template <class T>
void List<T>::push_back(const T &value)
{
    Node* body = new Node();
    body->value_node = new Value_node( value );
    insert_before_position(iterator(tail), body);
}

template<class T>
void List<T>::push_back(T &&value)
{
    Node* body = new Node();
    body->value_node = new Value_node( std::move( value ) );
    insert_before_position(iterator(tail), body);
}

template <class T>
void List<T>::push_front(const T &value)
{
    Node* body = new Node();
    body->value_node = new Value_node( value );
    insert_before_position(iterator(head->next), body);
}

template <class T>
void List<T>::push_front(T&& value)
{
    Node* body = new Node();
    body->value_node = new Value_node( std::move(value) );
    insert_before_position(iterator(head->next), body);
}

template <class T>
void List<T>::pop_back()
{
    delete_node_in_position(const_iterator(tail->previous));
}

template <class T>
void List<T>::pop_front()
{
    delete_node_in_position(const_iterator(head->next));
}

template <class T>
void List<T>::reverse()
{
    if( size_ == 0 || size_ == 1 )
        return;
    Node* body = tail->previous;
    body->next = nullptr;
    delete tail;
    head->next->previous = nullptr;
    delete head;
    Node* new_head = new Node;
    Node* new_tail = new Node;
    Node* next_node = body->previous;
    Node* prev = new_head;
    for( size_t i = 0; i < size_; ++i )
    {
        prev->next = body;
        body->previous = prev;
        body->next = next_node;
        prev = body;
        body = next_node;
        if( next_node->previous != nullptr )
            next_node = next_node->previous;
        else
            next_node = new_tail;
    }
    prev->next = new_tail;
    new_tail->previous = prev;
    head = new_head;
    tail = new_tail;
}

template <class T>
void List<T>::unique()
{
    if( size_ == 0 || size_ == 1 )
        return;
    Node* current_node = head->next;
    Node* next_node = current_node->next;
    Node* node_to_be_deleted = nullptr;
    size_t new_size = 1;
    while (next_node != tail)
    {
        if( current_node->value_node->value != next_node->value_node->value )
        {
            current_node = current_node->next;
            next_node = next_node->next;
            ++new_size;
        }
        else
        {
            node_to_be_deleted = next_node;
            next_node = next_node->next;
            next_node->previous = current_node;
            current_node->next = next_node;
            node_to_be_deleted->next = nullptr;
            node_to_be_deleted->previous = nullptr;
            delete node_to_be_deleted->value_node;
            delete node_to_be_deleted;
        }
    }
    size_ = new_size;
}

template <class T>
template <class... Args>
void List<T>::emplace(List::const_iterator position, Args &&... args)
{
    Node* body = new Node();
    body->value_node = new Value_node( T(std::forward<Args>(args)...)  );
    insert_before_position(position, body);
}

template <class T>
template <class... Args>
void List<T>::emplace_back(Args &&... args)
{
    Node* body = new Node();
    body->value_node = new Value_node(T(std::forward<Args>(args)...) );
    insert_before_position(const_iterator(tail), body);
}

template <class T>
template <class... Args>
void List<T>::emplace_front(Args &&... args)
{
    Node* body = new Node();
    body->value_node = new Value_node( T(std::forward<Args>(args)...)  );
    insert_before_position(const_iterator(head->next), body);
}

#endif //OOP_LIST_LIST_H
