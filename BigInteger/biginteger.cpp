#ifndef BIGINTEGER_BIGINTEGER_H
#define BIGINTEGER_BIGINTEGER_H

#include <string>
#include <iostream>
#include <cassert>

class BigInteger {
public:

    BigInteger() : buffer(""), digits_count(0), sign(true) {}

    BigInteger(const int &num);

    BigInteger(const BigInteger &other) : buffer(other.buffer), digits_count(other.digits_count), sign(other.sign) {}

    BigInteger(const std::string str);

    ~BigInteger() = default;

    BigInteger &operator=(const int &num);

    BigInteger &operator=(const BigInteger &other);

    BigInteger &operator=(const std::string &str);

    bool operator<(const BigInteger &other) const;

    bool operator==(const BigInteger &other) const;

    bool operator!=(const BigInteger &other) const;

    bool operator>=(const BigInteger &other) const;

    bool operator>(const BigInteger &other) const;

    bool operator<=(const BigInteger &other) const;

    std::string toString() const;

    void input(const std::string &str);

    BigInteger abs() const;

    explicit operator int() const;

    explicit operator bool() const;

    BigInteger operator-() const;

    BigInteger &operator+=(const BigInteger &other);

    BigInteger &operator-=(const BigInteger &other);

    BigInteger &operator*=(const BigInteger &other);

    BigInteger &operator/=(const BigInteger &other);

    BigInteger &operator%=(const BigInteger &other);

    BigInteger &operator++();

    BigInteger &operator--();

    BigInteger operator++(int);

    BigInteger operator--(int);

private:

    std::string buffer;
    int digits_count;
    bool sign;

    void convert_from_string(const std::string &str);

    void convert_from_int(const int &num);

    void clear();

    void make_longer(int max_digits);

    void addition(const BigInteger &other);

    void subtraction(const BigInteger &other, int one);

    void multiplication_one_long(const BigInteger &other);

    void multiplication_long_one(const BigInteger &other);

    void multiplication_recursively(const BigInteger &other);

    void cut(BigInteger &a, BigInteger &b, int n_half) const;

    void divide(BigInteger &other);

};

void string_reverse(std::string &str);

std::istream &operator>>(std::istream &in, BigInteger &num);

std::ostream &operator<<(std::ostream &out, const BigInteger &num);

BigInteger operator+(const BigInteger &a, const BigInteger &b);

BigInteger operator-(const BigInteger &a, const BigInteger &b);

BigInteger operator*(const BigInteger &a, const BigInteger &b);

BigInteger operator/(const BigInteger &a, const BigInteger &b);

BigInteger operator%(const BigInteger &a, const BigInteger &b);

//// String reverse

void string_reverse(std::string &str) {
    int length = str.length();
    for (int i = 0; i < length / 2; ++i) {
        std::swap(str[i], str[length - 1 - i]);
    }
}

//// Input and output

std::istream &operator>>(std::istream &in, BigInteger &num) {
    std::string str;
    in >> str;
    num.input(str);
    return in;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &num) {
    std::string str = num.toString();
    out << str;
    return out;
}

//// Private methods

void BigInteger::convert_from_int(const int &num) {
    int number = num;
    if (number == 0) {
        digits_count = 1;
        buffer.clear();
        buffer += '0';
        sign = true;
        return;
    }

    if (number < 0) {
        sign = false;
        number *= -1;
    } else {
        sign = true;
    }

    buffer.clear();
    digits_count = 0;
    while (number > 0) {
        buffer += (number % 10) + '0';
        ++digits_count;
        number /= 10;
    }

    clear();

}

void BigInteger::clear() {
    int length = buffer.length();
    for (int i = length - 1; i >= 1; --i) {
        if (buffer[i] == '0') {
            buffer.pop_back();
        } else {
            break;
        }
    }
    digits_count = buffer.length();
}

void BigInteger::convert_from_string(const std::string &str) {
    bool any = false;
    if (str[0] == '-') {
        sign = false;
        any = true;
    } else if (str[0] == '+') {
        sign = true;
        any = true;
    }

    int stop;
    buffer.clear();
    if (any) {
        stop = 1;
        digits_count = str.length() - 1;
    } else {
        sign = true;
        digits_count = str.length();
        stop = 0;
    }

    for (int i = str.length() - 1; i >= stop; --i) {
        buffer += str[i];
    }

    clear();
    if (digits_count == 1 && buffer[0] == '0') {
        sign = true;
    }

}

void BigInteger::make_longer(int max_digit) {
    for (int i = digits_count; i < max_digit; ++i) {
        buffer += '0';
    }
}

void BigInteger::addition(const BigInteger &other) {
    int max_digit = std::max(digits_count, other.digits_count);
    make_longer(max_digit);

    int tmp = 0, add = 0;
    for (int i = 0; i < max_digit; ++i) {
        if (i < other.digits_count) {
            tmp = add + (buffer[i] - '0') + (other.buffer[i] - '0');
        } else {
            tmp = add + (buffer[i] - '0');
        }
        if (tmp > 9) {
            add = tmp / 10;
            buffer[i] = (tmp % 10) + '0';
        } else {
            add = 0;
            buffer[i] = tmp + '0';
        }
    }

    if (add > 0) {
        buffer += (add + '0');
        ++digits_count;
    }
}

void BigInteger::subtraction(const BigInteger &other, int one) {
    int max_digit = std::max(other.digits_count, digits_count);
    make_longer(max_digit);

    int tmp = 0, sub = 0;
    for (int i = 0; i < max_digit; ++i) {
        if (i < other.digits_count) {
            tmp = one * ((buffer[i] - '0') - (other.buffer[i] - '0')) + sub;
        } else {
            tmp = one * (buffer[i] - '0') + sub;
        }
        if (tmp >= 0) {
            sub = 0;
            buffer[i] = tmp + '0';
        } else {
            sub = -1;
            buffer[i] = (tmp + 10) + '0';
        }
    }

    clear();
}

void BigInteger::multiplication_one_long(const BigInteger &other) {
    int mult = buffer[0] - '0';
    buffer.clear();
    int tmp = 0, add = 0;
    digits_count = other.digits_count;

    for (int i = 0; i < other.digits_count; ++i) {
        tmp = (other.buffer[i] - '0') * mult + add;
        if (tmp > 9) {
            add = tmp / 10;
            buffer += (tmp % 10) + '0';
        } else {
            add = 0;
            buffer += tmp + '0';
        }
    }

    if (add > 0) {
        buffer += add + '0';
        ++digits_count;
    }
    clear();
}

void BigInteger::multiplication_long_one(const BigInteger &other) {
    int mult = other.buffer[0] - '0', add = 0, tmp = 0;
    for (int i = 0; i < digits_count; ++i) {
        tmp = (buffer[i] - '0') * mult + add;
        if (tmp > 9) {
            add = tmp / 10;
            buffer[i] = (tmp % 10) + '0';
        } else {
            add = 0;
            buffer[i] = tmp + '0';
        }
    }

    if (add > 0) {
        buffer += add + '0';
        ++digits_count;
    }

    clear();
}

void BigInteger::multiplication_recursively(const BigInteger &other) {
    BigInteger a = 0, b = 0, c = 0, d = 0;
    int n_half = std::max(digits_count, other.digits_count) / 2;
    std::string str;
    this->cut(a, b, n_half);
    other.cut(c, d, n_half);

    BigInteger ac = a * c;
    BigInteger bd = b * d;
    BigInteger bc_ad = (a + b) * (c + d) - (ac + bd);
    BigInteger mid, front;

    if (bc_ad != 0) {
        mid.digits_count = n_half + bc_ad.digits_count;
        for (int i = 0; i < n_half; ++i)
            mid.buffer += '0';
        mid.buffer += bc_ad.buffer;
    } else {
        mid = 0;
    }

    if (ac != 0) {
        front.digits_count = 2 * n_half + ac.digits_count;
        for (int i = 0; i < 2 * n_half; ++i) {
            front.buffer += '0';
        }
        front.buffer += ac.buffer;
    } else {
        front = 0;
    }

    *this = bd + mid + front;
}

void BigInteger::cut(BigInteger &a, BigInteger &b, int n_half) const {
    std::string str;
    bool short_ = false;
    int i;
    for (i = 0; i < n_half; ++i) {
        if (i < digits_count) {
            str += buffer[i];
        } else {
            short_ = true;
            break;
        }
    }

    string_reverse(str);
    b.convert_from_string(str);
    if (short_ || digits_count == n_half) {
        a = 0;
    } else {
        str.clear();
        for (int i = digits_count - 1; i >= n_half; --i) {
            str += buffer[i];
        }
        a.convert_from_string(str);
    }
}

void BigInteger::divide(BigInteger &other) {
    BigInteger sub = 0, tmp = 0;
    int times = 0;
    std::string ans;
    std::string str;

    for (int i = digits_count - 1; i >= 0; --i) {
        str += buffer[i];
        tmp.convert_from_string(str);
        times = 0;
        while (tmp - other >= 0) {
            ++times;
            tmp -= other;
        }
        ans += times + '0';
        str = tmp.toString();
    }

    this->convert_from_string(ans);
}

//// Constructors

BigInteger::BigInteger(const int &num) {
    convert_from_int(num);
}

BigInteger::BigInteger(const std::string str) {
    convert_from_string(str);
}

//// Copy operators

BigInteger &BigInteger::operator=(const int &num) {
    convert_from_int(num);
    return *this;
}

BigInteger &BigInteger::operator=(const BigInteger &other) {
    buffer = other.buffer;
    digits_count = other.digits_count;
    sign = other.sign;
    return *this;
}

BigInteger &BigInteger::operator=(const std::string &str) {
    convert_from_string(str);
    return *this;
}

//// Other public methods

void BigInteger::input(const std::string &str) {
    convert_from_string(str);
}

std::string BigInteger::toString() const {
    std::string str;
    if (!sign) {
        str += '-';
    }
    for (int i = digits_count - 1; i >= 0; --i) {
        str += buffer[i];
    }

    return str;
}

BigInteger BigInteger::abs() const {
    BigInteger a = *this;
    a.sign = true;
    return a;
}

//// Implicit conversion operators

BigInteger::operator bool() const {
    return *this != 0;
}

BigInteger::operator int() const {
    int ans = 0;
    for (int i = 0, exp = 1; i < digits_count; ++i, exp *= 10) {
        ans += exp * (buffer[i] - '0');
    }

    if (!sign) {
        ans *= -1;
    }

    return ans;
}

//// Comparison operators

bool BigInteger::operator<(const BigInteger &other) const {
    if (sign == other.sign) {
        if (digits_count < other.digits_count) {
            return sign;
        } else if (digits_count > other.digits_count) {
            return !sign;
        } else {
            for (int i = digits_count - 1; i >= 0; --i) {
                if (buffer[i] == other.buffer[i]) {
                    continue;
                }
                return sign == (buffer[i] - '0' < other.buffer[i] - '0');
            }
            return false;
        }
    } else {
        return !sign;
    }
}

bool BigInteger::operator==(const BigInteger &other) const {
    return !(*this < other) && !(other < *this);
}

bool BigInteger::operator!=(const BigInteger &other) const {
    return *this < other || other < *this;
}

bool BigInteger::operator>=(const BigInteger &other) const {
    return !(*this < other);
}

bool BigInteger::operator>(const BigInteger &other) const {
    return other < *this;
}

bool BigInteger::operator<=(const BigInteger &other) const {
    return !(*this > other);
}

//// Arithmetical operators

BigInteger BigInteger::operator-() const {
    BigInteger tmp = *this;
    if (tmp != 0) {
        tmp.sign = !sign;
    }
    return tmp;
}

BigInteger &BigInteger::operator+=(const BigInteger &other) {
    if (other == 0) {
        return *this;
    }
    if (sign == other.sign) {
        addition(other);
    } else {
        *this -= -other;
    }

    clear();
    return *this;
}

BigInteger &BigInteger::operator-=(const BigInteger &other) {
    if (other == 0) {
        return *this;
    }

    if (sign == other.sign) {
        int one = 0;
        if (*this == other) {
            *this = 0;
        } else if (*this < other) {
            if (sign) {
                one = -1;
            } else {
                one = 1;
            }
            sign = false;
        } else {
            if (sign) {
                one = 1;
            } else {
                one = -1;
            }
            sign = true;
        }
        subtraction(other, one);
    } else {
        *this += -other;
    }

    clear();
    return *this;
}

BigInteger operator+(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c += b;
    return c;
}

BigInteger operator-(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c -= b;
    return c;
}

BigInteger &BigInteger::operator*=(const BigInteger &other) {
    if (*this == 0 || other == 0) {
        *this = 0;
        return *this;
    }

    bool this_sign = sign = (sign == other.sign);

    if (digits_count == 1) {
        multiplication_one_long(other);
    } else if (other.digits_count == 1) {
        multiplication_long_one(other);
    } else {
        multiplication_recursively(other);
    }

    sign = this_sign;
    return *this;
}

BigInteger operator*(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c *= b;
    return c;
}

BigInteger &BigInteger::operator/=(const BigInteger &other) {
    if (other == 0) {
        assert(false);
    }

    bool this_sign = (sign == other.sign);
    BigInteger abs_other = other.abs();
    sign = true;

    if (abs_other > *this) {
        *this = 0;
        return *this;
    } else if (abs_other == *this) {

        *this = 1;
        sign = this_sign;
        return *this;
    }

    divide(abs_other);
    sign = this_sign;
    return *this;
}

BigInteger operator/(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c /= b;
    return c;
}

BigInteger &BigInteger::operator%=(const BigInteger &other) {
    if (other == 0) {
        assert(false);
    }
    if (*this == 0) {
        return *this;
    }

    BigInteger abs_this;
    bool neg = false;
    if (*this < 0) {
        abs_this = -*this;
        neg = true;
    } else {
        abs_this = *this;
    }

    *this = abs_this - (abs_this / other) * other;

    if (*this == 0) {
        return *this;
    } else if (neg) {
        sign = false;
    }
    return *this;
}

BigInteger operator%(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a;
    c %= b;
    return c;
}

BigInteger &BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger BigInteger::operator++(int) {
    BigInteger c = *this;
    *this += 1;
    return c;
}

BigInteger &BigInteger::operator--() {
    *this -= 1;
    return *this;
}

BigInteger BigInteger::operator--(int) {
    BigInteger c = *this;
    *this -= 1;
    return c;
}



//// Rational

class Rational {
public:

    Rational() : numerator(0), denominator(1), sign(true) {}

    Rational(const int &num);

    Rational(const BigInteger &num);

    Rational(const Rational &other) : numerator(other.numerator), denominator(other.denominator), sign(other.sign) {}

    ~Rational() = default;

    Rational &operator=(const int &num);

    Rational &operator=(const BigInteger &num);

    Rational &operator=(const Rational &other);

    bool operator<(const Rational &other) const;

    bool operator==(const Rational &other) const;

    bool operator!=(const Rational &other) const;

    bool operator>=(const Rational &other) const;

    bool operator>(const Rational &other) const;

    bool operator<=(const Rational &other) const;

    std::string toString() const;

    std::string asDecimal(size_t precision = 0) const;

    explicit operator double() const;

    Rational operator-() const;

    Rational &operator+=(const Rational &other);

    Rational &operator-=(const Rational &other);

    Rational &operator*=(const Rational &other);

    Rational &operator/=(const Rational &other);

private:

    BigInteger numerator;
    BigInteger denominator;
    bool sign;

    void cut_back();

};

BigInteger GCD(const BigInteger &a, const BigInteger &b);

Rational operator+(const Rational &a, const Rational &b);

Rational operator-(const Rational &a, const Rational &b);

Rational operator*(const Rational &a, const Rational &b);

Rational operator/(const Rational &a, const Rational &b);

//// Greatest common divisor

BigInteger GCD(const BigInteger &a, const BigInteger &b) {
    BigInteger c = a.abs();
    BigInteger d = b.abs();
    if (c == 0)
        return d;
    if (d == 0)
        return c;
    if (d < c)
        return GCD(d, c % d);
    else
        return GCD(c, d % c);
}

//// Private methods

void Rational::cut_back() {
    if (numerator == 0)
        sign = true;
    if (denominator == 1)
        return;
    BigInteger gcd = GCD(denominator, numerator);
    if (gcd != 1) {
        denominator /= gcd;
        numerator /= gcd;
    }
}

//// Constructors

Rational::Rational(const int &num) {
    sign = true;
    numerator = num;
    denominator = 1;
    if (num < 0) {
        sign = false;
        numerator *= -1;
    }
}

Rational::Rational(const BigInteger &num) {
    sign = true;
    numerator = num;
    denominator = 1;
    if (num < 0) {
        sign = false;
        numerator *= -1;
    }
}

Rational &Rational::operator=(const int &num) {
    sign = true;
    numerator = num;
    denominator = 1;
    if (num < 0) {
        sign = false;
        numerator *= -1;
    }
    return *this;
}

Rational &Rational::operator=(const BigInteger &num) {
    sign = true;
    numerator = num;
    denominator = 1;
    if (num < 0) {
        sign = false;
        numerator *= -1;
    }
    return *this;
}

Rational &Rational::operator=(const Rational &other) {
    sign = other.sign;
    numerator = other.numerator;
    denominator = other.denominator;
    return *this;
}

//// Comparison operators

bool Rational::operator<(const Rational &other) const {
    if (sign == other.sign)
        return ((numerator * other.denominator) < (other.numerator * denominator)) == sign;
    else
        return !sign;
}

bool Rational::operator==(const Rational &other) const {
    return !(*this < other) && !(other < *this);
}

bool Rational::operator!=(const Rational &other) const {
    return *this < other || other < *this;
}

bool Rational::operator>=(const Rational &other) const {
    return !(*this < other);
}

bool Rational::operator>(const Rational &other) const {
    return other < *this;
}

bool Rational::operator<=(const Rational &other) const {
    return !(other < *this);
}

//// Arithmetical operators

Rational Rational::operator-() const {
    Rational a = *this;
    if (a.numerator != 0)
        a.sign = !sign;
    return a;
}

Rational &Rational::operator+=(const Rational &other) {
    if (other == *this) {
        *this *= 2;
        cut_back();
        return *this;
    }

    if (sign == other.sign)
        numerator = numerator * other.denominator + denominator * other.numerator;
    else if (sign)
        numerator = numerator * other.denominator - denominator * other.numerator;
    else
        numerator = other.denominator * numerator - other.numerator * denominator;

    if (numerator < 0) {
        numerator *= -1;
        sign = !sign;
    }
    denominator *= other.denominator;
    cut_back();
    return *this;
}

Rational &Rational::operator-=(const Rational &other) {
    if (other == *this) {
        *this = 0;
        return *this;
    }
    return *this += -other;
}

Rational operator+(const Rational &a, const Rational &b) {
    Rational c = a;
    c += b;
    return c;
}

Rational operator-(const Rational &a, const Rational &b) {
    Rational c = a;
    c -= b;
    return c;
}

Rational &Rational::operator*=(const Rational &other) {
    sign = (sign == other.sign);
    numerator *= other.numerator;
    denominator *= other.denominator;
    cut_back();
    return *this;
}

Rational operator*(const Rational &a, const Rational &b) {
    Rational c = a;
    c *= b;
    return c;
}

Rational &Rational::operator/=(const Rational &other) {
    if (other.numerator == 0)
        assert(false);
    if (*this == other) {
        *this = 1;
        return *this;
    }

    sign = (sign == other.sign);
    numerator *= other.denominator;
    denominator *= other.numerator;
    cut_back();
    return *this;
}

Rational operator/(const Rational &a, const Rational &b) {
    Rational c = a;
    c /= b;
    return c;
}

std::string Rational::toString() const {
    std::string ans;
    if (!sign)
        ans += '-';
    ans += numerator.toString();
    if (denominator != 1) {
        ans += '/';
        ans += denominator.toString();
    }
    return ans;
}

std::string Rational::asDecimal(size_t precision) const {
    std::string ans;
    if (!sign)
        ans += '-';
    BigInteger integer_part = numerator / denominator;
    BigInteger fractional_part = numerator % denominator;

    ans += integer_part.toString();

    if (precision > 0) {
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

Rational::operator double() const {
    double ans = 0;
    int num = int(numerator), dem = int(denominator);
    double num_db = num, dem_db = dem;
    ans = num_db / dem_db;
    if (!sign)
        ans *= -1;
    return ans;
}

#endif //BIGINTEGER_BIGINTEGER_H
