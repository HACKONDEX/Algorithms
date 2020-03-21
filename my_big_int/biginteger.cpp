#include <iostream>
#include <string>
#include <cassert>
#include <vector>


class BigInteger
{
public:

    BigInteger() : buffer( "" ), digits_count( 0 ), sign( true ){}
    BigInteger( const int& num );
    BigInteger( const std::string& str );
    BigInteger( const BigInteger& other ) : buffer( other.buffer ), sign( other.sign ), digits_count( other.digits_count ) {}

    ~BigInteger() = default;

    bool operator < ( const BigInteger& other) const;
    bool operator == ( const BigInteger& other ) const;
    bool operator != ( const BigInteger& other ) const;
    bool operator >= ( const BigInteger& other ) const;
    bool operator > ( const BigInteger& other ) const;
    bool operator <= ( const BigInteger& other ) const;

    std::string toString() const;
    void input( const std::string& str );
    BigInteger abs() const;

    explicit operator int() const;
    explicit operator bool() const;

    BigInteger& operator = ( const BigInteger& other );
    BigInteger& operator = ( const int& num );
    BigInteger& operator = ( const std::string& str );

    BigInteger operator - () const;
    BigInteger& operator += ( const BigInteger& other );
    BigInteger& operator -= ( const BigInteger& other );
    BigInteger& operator *= ( const BigInteger& other );
    BigInteger& operator /= ( const BigInteger& other );
    BigInteger& operator %= ( const BigInteger& other );
    BigInteger& operator ++ ();
    BigInteger& operator -- ();
    BigInteger operator ++ ( int );
    BigInteger operator -- ( int );

private:

    std::string buffer;
    int digits_count;
    bool sign;

    void convert_from_int( int num );
    void convert_from_string( const std::string& str );
    void convert_from_reversed_string( const std::string& str );

    void clear();

    void subtraction( const BigInteger& other, int one, int max_digit );

    void multiplication_one_long( const BigInteger& other );
    void multiplication_long_one( const BigInteger& other );
    void multiplication_recursively( const BigInteger& other );

};

std::istream& operator >> ( std::istream& in, BigInteger& num );
std::ostream& operator << (std::ostream& out, const BigInteger& num );
BigInteger operator + ( const BigInteger& a, const BigInteger& b );
BigInteger operator - ( const BigInteger& a, const BigInteger& b );
BigInteger operator * ( const BigInteger& a, const BigInteger& b );
BigInteger operator / ( const BigInteger& a, const BigInteger& b );
BigInteger operator % ( const BigInteger& a, const BigInteger& b );


////Private methods

void BigInteger::convert_from_int ( int num )
{
    buffer.clear();
    digits_count = 0;
    sign = true;
    if( num < 0 )
    {
        sign = false;
        num *= -1;
    }
    do
    {
        buffer += ( num % 10 ) + '0';
        num /= 10;
        ++digits_count;
    }
    while( num > 0 );
}//// V

void BigInteger::convert_from_string( const std::string& str )
{
    buffer.clear();
    digits_count = 0;
    bool any_sign = false;
    sign = true;
    if( str[0] == '-' )
    {
        if( str[1] != '0' )
        {
            sign = false;
        }
        any_sign = true;
    }
    else if( str[0] == '+' )
    {
        any_sign = true;
    }
    if( any_sign )
    {
        digits_count = str.length() - 1;
        for( int i = digits_count; i >=1; --i )
            buffer += str[i];
    }
    else
    {
        digits_count = str.length();
        for( int i = digits_count - 1; i >=0; --i )
            buffer += str[i];
    }

    if( digits_count != 1)
    {
        int i = digits_count - 1;
        while( i >= 1)
        {
            if( buffer[i] == '0')
            {
                --digits_count;
                --i;
                buffer.pop_back();
            }
            else
                break;
        }
    }
}//// V

void BigInteger::convert_from_reversed_string( const std::string& str )
{
    sign = true;
    digits_count = str.length();
    buffer.clear();
    buffer = str;

    if( digits_count != 1)
    {
        int i = digits_count - 1;
        while( i >= 1)
        {
            if( buffer[i] == '0')
            {
                --digits_count;
                --i;
                buffer.pop_back();
            }
            else
                break;
        }
    }

}

void BigInteger::clear()
{
    int length = buffer.length();
    for( int i = length - 1; i >= digits_count; --i )
        buffer.pop_back();
}//// V

void BigInteger::subtraction( const BigInteger& other, int one, int max_digit  )
{
    for( int i = digits_count; i < other.digits_count; ++i )
        buffer += '0';

    int sub = 0, tmp =0, zero_count = 0;
    digits_count = 0;

    for( int i = 0; i < max_digit; ++i )
    {
        if( i < other.digits_count )
            tmp = one * ( ( buffer[i] - '0' ) - ( other.buffer[i] - '0' ) ) + sub;
        else
            tmp = ( buffer[i] - '0' ) + sub;
        if( tmp > 0 )
        {
            buffer[i] = tmp + '0';
            digits_count += zero_count + 1;
            zero_count = 0;
            sub = 0;
        }
        else if ( tmp == 0 )
        {
            buffer[i] = '0';
            ++zero_count;
            sub = 0;
        }
        else
        {
            sub = -1;
            tmp += 10;
            buffer[i] = tmp + '0';
            digits_count += zero_count + 1;
            zero_count = 0;
        }
    }

    clear();
}

void BigInteger::multiplication_one_long( const BigInteger& other )
{
    int und = buffer[0] - '0';
    buffer.pop_back();
    digits_count = other.digits_count;
    int add = 0, tmp = 0;
    for( int i = 0; i < other.digits_count; ++i )
    {
        tmp = und * ( other.buffer[i] - '0' ) + add;
        if( tmp > 9 )
        {
            add = tmp / 10;
            buffer += ( tmp % 10 ) + '0';
        }
        else
        {
            add = 0;
            buffer += tmp + '0';
        }
    }

    if( add > 0 )
    {
        buffer += add + '0';
        ++digits_count;
    }
}

void BigInteger::multiplication_long_one( const BigInteger& other )
{
    int und = other.buffer[0] - '0';
    int add = 0, tmp = 0;

    for( int i = 0; i < digits_count; ++i )
    {
        tmp = ( buffer[i] - '0' ) * und + add;
        if( tmp > 9 )
        {
            add = tmp / 10;
            buffer[i] = ( tmp % 10 ) + '0';
        }
        else
        {
            add = 0;
            buffer[i] = tmp + '0';
        }
    }

    if( add > 0 )
    {
        buffer += add + '0';
        ++digits_count;
    }
}

void BigInteger::multiplication_recursively( const BigInteger& other )
{
    BigInteger a, b, c, d;
    int n_half = std::max( digits_count, other.digits_count ) / 2;
    bool short_ = false;


    std::string str;

    for( int i = 0; i < n_half; ++i )
    {
        if( i < digits_count )
            str += buffer[i];
        else
        {
            short_ = true;
            break;
        }
    }
    b.convert_from_reversed_string( str );
    str.clear();

    if( short_ )
    {
        a.buffer += '0';
        a.digits_count = 1;
    }
    else
    {
        for( int i = n_half; i < digits_count; ++i )
        {
            a.buffer += buffer[i];
            ++a.digits_count;
        }
    }

    short_ = false;
    for( int i = 0; i < n_half; ++i )
    {
        if( i < other.digits_count )
            str += other.buffer[i];
        else
        {
            short_ = true;
            break;
        }
    }

    d.convert_from_reversed_string( str );

    if( short_ )
    {
        c.buffer += '0';
        c.digits_count = 1;
    }
    else
    {
        for( int i = n_half; i < other.digits_count; ++i )
        {
            c.buffer += other.buffer[i];
            ++c.digits_count;
        }
    }

    BigInteger ac = a * c;
    BigInteger bd = b * d;
    BigInteger bc_ad = ( a + b ) * ( c + d ) - ( ac + bd );

    BigInteger mid, front;

    if( bc_ad != 0 ) {
        mid.digits_count = n_half + bc_ad.digits_count;
        for (int i = 0; i < n_half; ++i)
            mid.buffer += '0';
        mid.buffer += bc_ad.buffer;
    }
    else
        mid = 0;

    if( ac != 0 ) {
        front.digits_count = 2 * n_half + ac.digits_count;
        for (int i = 0; i < 2 * n_half; ++i)
            front.buffer += '0';
        front.buffer += ac.buffer;
    }
    else
        ac = 0;

    *this = 0;
    *this += bd + mid + front;

}

//// Constructors

BigInteger::BigInteger( const int &num )
{
    this->sign = ( num >= 0 );
    this->digits_count = 0;
    convert_from_int( num );
}//// V

BigInteger::BigInteger(const std::string& str)
{
    this->digits_count = 0;
    this->sign = false;
    this->convert_from_string( str );
}//// V

//// Comparisons

bool BigInteger::operator < ( const BigInteger& other ) const
{
    if( sign == other.sign )
    {
        if( digits_count < other.digits_count )
            return sign;
        else if( digits_count > other.digits_count )
            return !sign;
        else
        {
            for( int i = digits_count - 1; i >= 0; --i )
            {
                if ( buffer[i] == other.buffer[i] )
                    continue;
                return sign == ( static_cast<int>( buffer[i] ) < static_cast<int>( other.buffer[i] ) );
            }
            return false;
        }
    }
        return !sign;
}//// V

bool BigInteger::operator == ( const BigInteger& other ) const
{
    return !( *this < other ) && !( other < *this );
}//// V

bool BigInteger::operator != ( const BigInteger& other ) const
{
    return *this < other || other < *this;
}//// V

bool BigInteger::operator >= ( const BigInteger& other ) const
{
    return !( *this < other );
}//// V

bool BigInteger::operator > ( const BigInteger& other ) const
{
    return other < *this;
}//// V

bool BigInteger::operator <= ( const BigInteger& other ) const
{
    return !( other < *this );
}//// V

//// Other public methods

void BigInteger::input( const std::string &str )
{
    this->convert_from_string( str );
}//// V

std::string BigInteger::toString() const
{
    std::string str;
    if( !this->sign && this->buffer[digits_count - 1] != '0' )
        str += '-';
    for( int i = digits_count -1; i >= 0; --i )
        str += buffer[i];
    return str;
}//// V

BigInteger BigInteger::abs() const
{
    BigInteger a = *this;
    a.sign = true;
    return a;
}


//// Input, Output

std::istream& operator >> ( std::istream& in, BigInteger& num )
{
    std::string str;
    in >> str;
    num.input( str );
    return in;
}//// V

std::ostream& operator << ( std::ostream& out, const BigInteger& num )
{
    std::string str = num.toString();
    out << str;
    return out;
}//// V

//// Copy operators

BigInteger& BigInteger::operator = ( const BigInteger& other )
{
    this->sign = other.sign;
    this->digits_count = other.digits_count;
    buffer.clear();
    this->buffer = other.buffer;
    return *this;
}//// V

BigInteger& BigInteger::operator = ( const int &num )
{
    this->convert_from_int( num );
    return* this;
}//// V

BigInteger& BigInteger::operator = ( const std::string& str )
{
    this->convert_from_string( str );
    return* this;
}//// V

//// Implicit conversion

BigInteger::operator bool() const
{
    if( this->digits_count == 0 )
        return false;
    return this->buffer[0] != '0';
}//// V

BigInteger::operator int() const
{
    int ans = 0;
    for( int i = this->digits_count - 1, exp = 1; i >= 0; i--, exp *= 10 )
        ans += ( this->buffer[i] + '0' ) * exp;
    if( this->sign )
        return ans;
    return -ans;
}//// V

//// Arithmetical operators

BigInteger BigInteger::operator - () const
{
    BigInteger tmp = *this;
    if( buffer[digits_count - 1] != '0' )
        tmp.sign = !tmp.sign;
    return tmp;
}

BigInteger& BigInteger::operator += ( const BigInteger& other )
{
    if( sign == other.sign )
    {
        int min_digit_count = std::min( digits_count, other.digits_count );
        int add = 0;
        int tmp = 0;
        for( int i = 0; i < min_digit_count; ++i )
        {
            tmp = ( buffer[i] - '0' ) + ( other.buffer[i] - '0') + add;
            if( tmp > 9 )
            {
                add = tmp / 10;
                buffer[i] = ( tmp % 10 ) + '0';
            }
            else
            {
                buffer[i] = tmp + '0';
                add = 0;
            }
        }

        if( digits_count == min_digit_count && other.digits_count != min_digit_count )
        {
            for( int i = min_digit_count; i < other.digits_count; ++i )
            {
                tmp = add + ( other.buffer[i] - '0' );
                if( tmp > 9 )
                {
                    add = tmp / 10;
                    buffer += ( tmp % 10 ) + '0';
                }
                else
                {
                    buffer += tmp + '0';
                    add = 0;
                }
            }
            digits_count = other.digits_count;
        }
        else
        {
            for( int i = min_digit_count; i < digits_count; ++i )
            {
                tmp = add + ( buffer[i] - '0' );
                if( tmp > 9 )
                {
                    add = tmp / 10;
                    buffer[i] = ( tmp % 10 ) + '0';
                }
                else
                {
                    buffer[i] = tmp + '0';
                    add = 0;
                }
            }
        }

        if( add > 0 )
        {
            buffer += add + '0';
            ++digits_count;
        }
    }
    else
    {
        if( other == 0 )
            return *this;
        *this -= -other;
    }
    return *this;
}

BigInteger& BigInteger::operator -= ( const BigInteger& other )
{
    if( sign == other.sign )
    {
        if( other == * this )
        {
            buffer[0] = '0';
            sign = true;
            digits_count = 1;
            clear();
            return *this;
        }
        if( sign )
        {
            if( *this >= other )
                subtraction( other, 1, digits_count );
            else
            {
                sign = false;
                subtraction( other, -1, other.digits_count );
            }
        }
        else
        {
            if( *this >= other )
            {
                sign = true;
                subtraction( other, -1, other.digits_count );
            }
            else
                subtraction(other, 1, digits_count);
        }
    }
    else
    {
        if( other == 0 )
            return *this;
        *this += - other;
    }
    return *this;
}

BigInteger operator + ( const BigInteger& a, const BigInteger& b )
{
    BigInteger c = a;
    c += b;
    return c;
}

BigInteger operator - ( const BigInteger& a, const BigInteger& b )
{
    BigInteger c = a;
    c -= b;
    return c;
}

BigInteger& BigInteger::operator *= ( const BigInteger& other )
{
    if( buffer[digits_count - 1] == '0')
        return *this;
    if( other.buffer[other.digits_count - 1] == '0')
    {
        buffer.clear();
        buffer += '0';
        digits_count = 1;
        sign = true;
        return *this;
    }

     bool new_sign = sign = ( sign == other.sign );

    if( digits_count == 1 )
    {
                multiplication_one_long( other );
                return *this;
    }
    else if( other.digits_count == 1 )
    {
        multiplication_long_one( other );
        return *this;
    }

    this->multiplication_recursively( other );
    sign = new_sign;
    return *this;
}

BigInteger operator * ( const BigInteger& a, const BigInteger& b )
{
    BigInteger c = a;
    c *= b;
    return c;
}

BigInteger& BigInteger::operator /= ( const BigInteger& other )
{
    if( other == 0 )
        assert( false );

    bool this_sign = ( sign == other.sign );
    BigInteger abs_other = other.abs();
    sign = true;

    if( abs_other > *this )
    {
        buffer.clear();
        buffer += '0';
        digits_count = 1;
        return *this;
    }
    else if( abs_other == *this )
    {
        buffer.clear();
        buffer += '1';
        digits_count = 1;
        return *this;
    }

    BigInteger sub = 0, tmp = 0;
    int times = 0;
    std::string ans;
    std::string str;

    for( int i = digits_count - 1; i >= 0; --i )
    {
        str += buffer[i];
        tmp.convert_from_string( str );

        times = 0;
        while( tmp - abs_other >= 0 )
        {
            ++times;
            tmp -= abs_other;
        }

        ans += times + '0';


        str = tmp.toString();
    }

    this->convert_from_string( ans );
    sign = this_sign;
    return *this;
}

BigInteger operator / ( const BigInteger& a, const BigInteger& b )
{
    BigInteger c = a;
    c /= b;
    return c;
}

BigInteger& BigInteger::operator %= ( const BigInteger& other )
{
    if( other == 0 )
        assert( false );

    if( *this == 0 )
        return *this;
    BigInteger one = 1;

    BigInteger abs_this;
    bool neg = false;
    if( *this < 0 )
    {
        abs_this = -*this;
        neg = true;
    }
    else
    {
        abs_this = *this;
    }

    *this = abs_this - ( abs_this / other ) * other;

    if( *this == 0 )
        return *this;
    else if( neg )
        sign = false;
    return *this;
}

BigInteger operator % ( const BigInteger& a, const BigInteger& b )
{
    BigInteger c = a;
    c %= b;
    return c;
}

BigInteger& BigInteger::operator ++ ()
{
    *this += 1;
    return *this;
}

BigInteger BigInteger::operator ++ ( int )
{
    BigInteger c = *this;
    ++*this;
    return c;
}

BigInteger& BigInteger::operator -- ()
{
    *this -= 1;
    return *this;
}

BigInteger BigInteger::operator -- ( int )
{
    BigInteger c = *this;
    --*this;
    return c;
}





//// Rational class

class Rational
{
public:

    Rational() : numerator(), denominator(), sign( true ) {}
    Rational( const int& num );
    Rational( const BigInteger& num );
    Rational( const Rational& other ) : numerator( other.numerator ), denominator( other.denominator ), sign( other.sign ) {}

    ~Rational() = default;

    Rational& operator = ( const int& num );
    Rational& operator = ( const BigInteger& num );
    Rational& operator = ( const Rational& other );

    bool operator < ( const Rational& other ) const;
    bool operator == ( const Rational& other) const;
    bool operator != ( const Rational& other ) const;
    bool operator >= ( const Rational& other ) const;
    bool operator > ( const Rational& other ) const;
    bool operator <= ( const Rational& other ) const;

    Rational operator - () const;
    Rational& operator += ( const Rational& other );
    Rational& operator -= ( const Rational& other );
    Rational& operator *= ( const Rational& other );
    Rational& operator /= ( const Rational& other );

    std::string toString() const;
    std::string asDecimal( size_t precision = 0 ) const;
    explicit operator double() const;

private:

    BigInteger numerator;
    BigInteger denominator;
    bool sign;

    void cut_back();


};

BigInteger GCD( const BigInteger& a, const BigInteger& b );
Rational operator + ( const Rational& a, const Rational& b );
Rational operator - ( const Rational& a, const Rational& b );
Rational operator * ( const Rational& a, const Rational& b );
Rational operator / ( const Rational& a, const Rational& b );

//// Greatest common divisor

BigInteger GCD( const BigInteger& a, const BigInteger& b )
{
    BigInteger c = a.abs();
    BigInteger d = b.abs();
    if( c == 0 )
        return d;
    if( d == 0)
        return c;
    if( d < c )
        return  GCD( d, c % d );
    else
        return  GCD( c, d % c );
}

//// Constructors

Rational::Rational( const int &num )
{
    sign = true;
    numerator = num;
    denominator = 1;
    if ( num < 0 )
    {
        sign = false;
        numerator *= -1;
    }
}

Rational::Rational( const BigInteger& num ) {
    sign = true;
    numerator = num;
    denominator = 1;
    if (num < 0)
    {
        sign = false;
        numerator *= -1;
    }
}

//// Private methods

void Rational::cut_back()
{
    if( denominator == 1 )
        return;
    BigInteger gcd = GCD( denominator, numerator );
    if( gcd != 1 )
    {
        denominator /= gcd;
        numerator /= gcd;
    }
}

//// Copy operators

Rational& Rational::operator = ( const int& num )
{
    sign = true;
    numerator = num;
    denominator = 1;
    if( num < 0 )
    {
        sign = false;
        numerator *= -1;
    }
    return *this;
}

Rational& Rational::operator = ( const BigInteger& num )
{
    sign = true;
    numerator = num;
    denominator = 1;
    if( num < 0 )
    {
        sign = false;
        numerator *= -1;
    }
    return *this;
}

Rational& Rational::operator = ( const Rational& other )
{
    sign = other.sign;
    numerator = other.numerator;
    denominator = other.denominator;
    return *this;
}

//// Comparison operators

bool Rational::operator < ( const Rational& other ) const
{
    if( sign == other.sign )
        return ( ( numerator * other.denominator ) < ( other.numerator * denominator ) == sign );
    else
        return !sign;
}

bool Rational::operator == ( const Rational& other ) const
{
    return !( *this < other ) && !( other < *this );
}

bool Rational::operator != ( const Rational& other ) const
{
    return *this < other || other < *this;
}

bool Rational::operator >= ( const Rational& other ) const
{
    return !( *this < other );
}

bool Rational::operator > ( const Rational& other ) const
{
    return other < *this;
}

bool Rational::operator <= ( const Rational& other ) const
{
    return !( other < *this );
}

//// Arithmetical operators

Rational Rational::operator - () const
{
    Rational a = *this;
    if( a.numerator != 0 )
        a.sign = !sign;
    return a;
}

Rational& Rational::operator += ( const Rational &other )
{
    if( sign == other.sign )
        numerator = numerator * other.denominator + denominator * other.numerator;
    else if( sign )
        numerator = numerator * other.denominator - denominator * other.numerator;
    else
        numerator = denominator * other.numerator - numerator * other.denominator;

    if( numerator < 0 )
    {
        numerator *= -1;
        sign = false;
    }
    denominator *= other.denominator;
    cut_back();
    return *this;
}

Rational& Rational::operator -= ( const Rational& other )
{
    return *this += -other;
}

Rational operator + ( const Rational& a, const Rational& b )
{
    Rational c = a;
    c += b;
    return c;
}

Rational operator - ( const Rational& a, const Rational& b )
{
    Rational c = a;
    c -= b;
    return c;
}

Rational& Rational::operator *= ( const Rational& other )
{
    sign = ( sign == other.sign );
    numerator *= other.numerator;
    denominator *= other.denominator;
    cut_back();
    return *this;
}

Rational operator * ( const Rational& a, const Rational& b )
{
    Rational c = a;
    c *= b;
    return c;
}

Rational& Rational::operator /= ( const Rational& other )
{
    if( other.numerator == 0 )
        assert( false );

    sign = ( sign == other.sign );
    numerator *= other.denominator;
    denominator *= other.numerator;
    cut_back();
    return *this;
}

Rational operator / ( const Rational& a, const Rational& b )
{
    Rational c = a;
    c /= b;
    return c;
}

//// Public methods

std::string Rational::toString() const
{
    std::string ans;
    if ( !sign )
        ans += '-';
    ans += numerator.toString();
    if( denominator != 1 )
    {
        ans += '/';
        ans += denominator.toString();
    }

    return ans;
}

std::string Rational::asDecimal( size_t precision ) const
{
    std::string ans;
    if( !sign )
        ans += '-';
    BigInteger integer_part = numerator / denominator;
    BigInteger fractional_part = numerator % denominator;

    ans += integer_part.toString();

    if( precision > 0 )
    {
        ans += '.';
        while (precision > 0) {
            fractional_part *= 10;
            ans += (fractional_part / denominator).toString();
            fractional_part %= denominator;
            --precision;
        }
    }

    return ans;
}

Rational::operator double() const
{
    double ans = 0;
    std::string num_str = numerator.toString();
    std::string dem_str = denominator.toString();
    int num = 0, dem = 0;
    for( int i = num_str.length() - 1, j = 1; i >= 0 ; --i, j *= 10 )
        num += j * ( num_str[i] - '0' );
    for( int i = dem_str.length() - 1, j = 1; i >= 0; --i, j*=10 )
        dem += j * ( dem_str[i] - '0' );
    double num_db = num, dem_db = dem;
    ans = num_db / dem_db;
    return ans;
}






