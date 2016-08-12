
// clang++ -o factorial fast_factorial.cpp -std=c++11 -O3

#include <chrono>
#include "fast_factorial.hpp"

int main() {
    try {
        const uint32_t number = 1000000;
        auto s = std::chrono::high_resolution_clock::now();
        math::BigInt res = math::Factorial(number);
        auto e = std::chrono::high_resolution_clock::now();
        auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(e - s);
        std::cerr << "Time taken by Factorial(" << number << ") in nano seconds " << d.count() << " or " << d.count()/double(1000000000) << " seconds\n";

        /* Uncomment the below code to print factorial value onto screen. */

        /*
           s = std::chrono::high_resolution_clock::now();
           std::cout << res << '\n';
           e = std::chrono::high_resolution_clock::now();
           d = std::chrono::duration_cast<std::chrono::microseconds>(e - s);
           std::cerr << "Time taken to print " << d.count() << '\n';
        */

    } catch (std::exception &e) {
        std::cout << "Exception: " << e.what() << '\n';
    }

    return 0;
}
