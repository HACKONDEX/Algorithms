#include <iostream>
#include "biginteger.cpp"
int main() {



//    std::cout << ( a != e ) << std::endl;
//    std::cout << ( c != b ) << std::endl;
//    std::cout << ( a < c ) << std::endl;
//    std::cout << ( c < d ) << std::endl;
//    std::cout << ( d < a ) << std::endl;
//    std::cout << ( c < b ) << std::endl;
//    std::cout << ( e < c ) << std::endl;
//    std::cout << ( c < e ) << std::endl;
//    std::cout << ( f < b ) << std::endl

    BigInteger h;
    BigInteger j;

    std::cin >> j >> h;

//    std::cout << " j = " << j << std::endl;
//    std::cout << " h = " << h << std::endl;

//    BigInteger l = j * h;
//    std::cout << l << std::endl;

    j %= h;

    std::cout << " j = " << j << std::endl;
    std::cout << " h = " << h << std::endl;






    return 0;
}