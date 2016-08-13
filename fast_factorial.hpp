
#ifndef FAST_FACTORIAL_HPP
#define FAST_FACTORIAL_HPP

#include <vector>
#include <sstream>
#include <cstdint>
#include <string>
#include <cassert>
#include <algorithm>

#include <iostream>
#include <iterator>

/* Implementation of factorial in the Python3 way. Factorial is asymptotic by expressing the factorial result
   as products of prime powers. However, integer prime factorization isn't easy job. There are efficient
   algorithm to evaluate factorial. Refer the following link author website for details of such algorithms,
   
   http://www.luschny.de/math/factorial/binarysplitfact.html
   
   The current code is an implementation of above mentioned idea.
 
   To understand the code flow, start at Factorial() function.
   
   ++++++++++++Divide-Conquer Multiplication Algorithm++++++++++++

   We use polynomial representation of integral numbers in large base, preferbly power of 2. We are considering
   64-bit system with base of 2 ^ 31, where 31 is exponent (scale). The number can be expressed as,

   N = An * 2^(scale*n) + An-1 * 2^(scale*(n-1)) + ... + A0 * 2^(scale*0)
   
   i.e. 'scale' is representation weight and A0 ... An are co-efficients.

   These co-efficients can be stored in a collection. In our implementation, we store A0 at 0-index.

   We can think the number as slices of 32 integers in a collection. Index of collection corresponds to
   the weight in terms of (2^scale*index). We represent the least significant digit of number in 0 location
   followed by higher co-efficients on higher side of collection.
   
   We have used vector as our collection data structure.

   The beauty in using exact power of 2 simplifies multiplication as left shift. I didn't consider all 32
   bits of uint32_t to ignore carry (borrow) during intermediate addition (subtraction). We can use all 32 bits
   to represent each co-efficient, yet it complicates addition & subtraction during addition of intermediate
   results.

   Infact, BINARY_BASE can anything below 32. The algorithms used holds good.
   
   ++++++++++++Competetive Programming Case++++++++++++

   Most of the CP cases require truncating N! with a remainder operation over largest prime. The current idea
   is adopted from Python3 implementation, can be tweaked for CP needs.
*/

namespace math {

#if defined (__x86_64__)
    /* word_t is alias of uint32_t to hold maximum of 2^31-1 as co-efficient value.
     * dword_t is used during multiplication of two word_t data.
     */
    typedef uint32_t  word_t;
    typedef uint64_t dword_t;

    const static int32_t RECURSIVE_MUL_THRESHOLD = 120; // A trial & error limit, varies across platfroms.

    // Base - 2**31, leaving one bit for carry manipulation during positional addition & subtraction
    const static uint32_t BINARY_SHIFT  = 31;
    const static uint32_t BINARY_BASE   = word_t(1) << BINARY_SHIFT;
    const static uint32_t BINARY_MASK   = BINARY_BASE - 1;

    const static uint32_t DECIMAL_SHIFT = 9;
    const static uint32_t DECIMAL_BASE  = 1000 * 1000 * 1000;
#else
    #error "Undefined word_t on 32 bit or lower size platforms"
#endif

/* C++ vector along with move semantics are used to hold co-efficients. */
typedef std::vector<word_t> value_container_t;

/* Products of odd numbers in N! to be multiplied with 2^pexp, where 'pexp' is exponent of 2 in prime
   factors of N! Interestingly the following relation holds between population count (# set bits) and N.
   Using GCC intrinisics to found set bits at the expense of one instruction.
*/
inline
int32_t two_exponent(uint32_t n) {
    return n - __builtin_popcount(n);
}

/* A helper function to drop zero valued higher co-efficient(s) after every intermediate result. */
inline
void drop_leading_zeros(value_container_t &val) {
    while (!val.empty() && 0 == val.back())
        val.pop_back();
}

/* Splits input co-efficient vector into two parts @shift. */
inline
void split(const value_container_t &v, uint32_t shift, value_container_t &vh, value_container_t &vl) {
    const auto split_at = std::min<uint32_t>(shift, v.size());
    std::copy(v.begin(), v.begin()+split_at, std::back_inserter(vl));
    std::copy(v.begin()+split_at, v.end(), std::back_inserter(vh));
}

/* word_t used in both positional operations as we don't overflow (underflow) due to one bit left in 32-bits

   In the expression - using 3 multiplications and 4 additions,

   a*b = ah * bh * 2^2*scale + ((ah+al)*(bh+bl) - ah*bh - al*bl)2^scale + al*bl

   Since the weight is exact power of 2, co-efficients can be manipulated at their position of 'scale' digits
   or '2*scale' digits from the begining of final result. Each term can be placed at their position by doing
   positional addition or subtraction. Selection of 31 as shift (scale) helps in ignoring carry or borrow.

   We add (subtract) at position depending on the term being considered into result. Precalculate the result
   size in 'ra', place ah*bh at higher side of 'result', place al*bl at lower side of result and do positional
   add (subtract) the term (ah+al)*(bh_bl) or ah*bh or al*bl at 'scale' shift.
*/
void positional_add(value_container_t &ra, uint32_t position, uint32_t no_digits, value_container_t &rb) {

    word_t carry = 0; // not dword_t as we don't care of carry bit
    uint32_t i = 0;

    for (; i < rb.size(); i++) {
        carry += ra[position+i] + rb[i];
        ra[position+i] = word_t(carry & BINARY_MASK);
        carry >>= BINARY_SHIFT;
        carry &= 1;
    }

    for (; carry && i < no_digits; i++) {
        carry += ra[position+i];
        ra[position+i] = word_t(carry & BINARY_MASK);
        carry >>= BINARY_SHIFT;
        carry &= 1;
    }
}

void positional_sub(value_container_t &ra, uint32_t position, uint32_t no_digits, value_container_t &rb) {
    word_t borrow = 0;
    uint32_t i = 0;

    for (; i < rb.size(); i++) {
        borrow = ra[position+i] - rb[i] - borrow;
        ra[position+i] = word_t(borrow & BINARY_MASK);
        borrow >>= BINARY_SHIFT;
        borrow &= 1;
    }

    for (; borrow && i < no_digits; i++) {
        borrow = ra[position+i] - borrow;
        ra[position+i] = word_t(borrow & BINARY_MASK);
        borrow >>= BINARY_SHIFT;
        borrow &= 1;
    }
}

/* Adds two numbers reprented in polynomial notation. Since the weight is 2^31 result vector size can't 
   overflow as we have guard bit. We drop off leading zero co-efficients if needed.
*/
value_container_t add(value_container_t const *pa, value_container_t const *pb) {

    uint32_t asize = pa->size();
    uint32_t bsize = pb->size();

    // Assumption: asize <= bsize
    if (asize > bsize) {
        std::swap(pa, pb);
        std::swap(asize, bsize);
    }

    auto const &ra = *pa;
    auto const &rb = *pb;

    value_container_t result(bsize+1, 0);

    uint32_t i = 0;
    word_t carry = 0;
    for (; i < asize; i++) {
        carry += ra[i] + rb[i];
        result[i] = carry & BINARY_MASK;
        carry >>= BINARY_SHIFT;
    }

    for (; i < bsize; i++) {
        carry += rb[i];
        result[i] = carry & BINARY_MASK;
        carry >>= BINARY_SHIFT;
    }

    if (carry) {
        assert(i < result.size());
        result[i] = carry;
    }

    drop_leading_zeros(result);
    return std::move(result);
}

/* Old school algorithm:
   Calculates product in O(m * n) --- m, n are number of digits in a, b respectively.
*/
value_container_t traditional_multiplication(value_container_t const *pa, value_container_t const *pb) {
    // Assumption asize <= bsize
    // Keeping outer loop smaller results in less computation cost over the other way
    uint32_t asize = pa->size();
    uint32_t bsize = pb->size();

    if (asize > bsize) {
        std::swap(pa, pb);
        std::swap(asize, bsize);
    }

    /* Didn't find a better trick. */
    auto const &ra = *pa;
    auto const &rb = *pb;

    /* k is partial result index
     * 'ai' should be dword_t, to save precision in multiplication
     *
     * number of digits in result = log(a) + log(b) to the base 2. It can't exceed.
    */
    value_container_t result(asize+bsize, 0);
    for (uint32_t i = 0; i < asize; i++) {
        dword_t prod = 0;
        dword_t   ai = ra[i];
        uint32_t   k = i;

        for (uint32_t j = 0; j < bsize; j++, k++) {
            prod      +=    (result[k] + rb[j] * ai);
            result[k] =     word_t(prod & BINARY_MASK);
            prod      >>=   BINARY_SHIFT;
        }

        if (prod) {
            assert(k < result.size());
            result[k] += word_t(prod & BINARY_MASK);
        }
    }

    drop_leading_zeros(result);
    return std::move(result);
}

value_container_t recursive_multiplication(value_container_t const *pa, value_container_t const *pb) {

    /* a x b = ah x bh x 2^2xscale - ((ah+al)x(bh+bl) - ahxbh - alxbl) x 2^scale + al x bl
     *
     * x1 = ah x bh
     *
     * x2 = al x bl
     *
     * x3 = (ah+al)x(bh+bl) - ahxbh - alxbl
    */

    auto asize = pa->size();
    auto bsize = pb->size();

    if (asize > bsize) {
        std::swap(pa, pb);
        std::swap(asize, bsize);
    }

    if (asize <= RECURSIVE_MUL_THRESHOLD)
        return traditional_multiplication(pa, pb);

    auto const &ra = *pa;
    auto const &rb = *pb;

    uint32_t shift = bsize/2;
    value_container_t ah, al;
    split(ra, shift, ah, al);

    value_container_t bh, bl;
    split(rb, shift, bh, bl);

    value_container_t result(asize + bsize, 0);

    auto x1(recursive_multiplication(&ah, &bh));

    // Place (ah x bh) on the upper words of 'result'
    for (uint32_t i = 2*shift, j = 0; j < x1.size(); i++, j++)
        result[i] = x1[j];

    // Place (al x bl) on the lower words of 'result'
    auto x2(recursive_multiplication(&al, &bl));
    for (uint32_t i = 0, j = 0; j < x2.size(); i++, j++)
        result[i] = x2[j];

    // ['result' - (al x bl) * 2^scale] - can be done by shifted arithmeric.
    // Since we are leaving one bit, carry or borrow can be ignored.
    positional_sub(result, shift, result.size()-shift, x2);
    positional_sub(result, shift, result.size()-shift, x1);

    x1 = add(&ah, &al);
    x2 = add(&bh, &bl);
    auto x3 = recursive_multiplication(&x1, &x2);

    // ['result' + (ah+al)x(bh+bl) * 2^scale] - can be done by shifted arithmeric.
    // Refer function comments.
    positional_add(result, shift, result.size()-shift, x3);

    drop_leading_zeros(result);
    return std::move(result);
}

/* The divide-and-conquer algorithm and grade school algorithms to multiply two polynomial number are doing the
   actual work. Below class BigInt being used to abstract arbitary precision integer. A few needed constructors,
   conversion operators are added.
*/

struct BigInt {
    value_container_t value_;

    // Constructures
    BigInt() {}
    BigInt(value_container_t const &value) {
        value_.clear();
        std::copy(value.begin(), value.end(), std::back_inserter(value_));
    }

    // C++11 move ctor
    BigInt(BigInt&& rhs) noexcept : value_(std::move(rhs.value_)) { }

    // Copy ctor
    BigInt(BigInt const &rhs) {
        value_.clear();
        std::copy(rhs.value_.begin(), rhs.value_.end(), std::back_inserter(value_));
    }

    // Copy Assignment operator
    // Supports expressions like .. BigInt a = c ... where c is also BigInt
    BigInt& operator = (const BigInt &rhs) {
        if (this != &rhs) {
            value_.clear();
            std::copy(rhs.value_.begin(), rhs.value_.end(), std::back_inserter(value_));
        }

        return *this;
    }

    /* Conversion ctor-s
       supports expressions like BigInt(1654600160), BigInt(1654600165416L), etc...
       'const ref' to bind numeric literals.
       Added only required constructors to support only numeric types, template 'from_number'
       is used for code reuse.
    */
    BigInt(const int32_t &rhs) {
        from_number(rhs);
    }

    BigInt(const uint32_t &rhs) {
        from_number(rhs);
    }

    BigInt(const int64_t &rhs) {
        from_number(rhs);
    }

    BigInt(const uint64_t &rhs) {
        from_number(rhs);
    }

    // Conversion assignment to support expression like, BigInt i = 10
    void operator = (const int64_t &rhs) {
        value_.clear();
        from_number(rhs);
    }

    // Global operator overloading - stream printing
    friend std::ostream& operator << (std::ostream &os, const BigInt &v) {
        os << v.to_sting();
        return os;
    }

    // Integer arithmetic - only multiplication is implemented
    const BigInt operator * (const BigInt &rhs) const {
        if (value_.size() == 0 || rhs.value_.size() == 0)
            return BigInt();

        return BigInt(recursive_multiplication(&value_, &(rhs.value_)));
    }

    /* Multiplies given number by 2^shift. We should have implemented it as saperate function.
     *
     * We shift all bits by 'shift' positions, and adjust bits if shift is not multible of 31.
     */
    BigInt left_shift(uint32_t shift) {
        uint32_t words_to_shift = shift/BINARY_SHIFT;
        uint32_t residue_bits   = shift - words_to_shift*BINARY_SHIFT;

        const auto &ra = value_;
        value_container_t result(ra.size() + words_to_shift + (residue_bits != 0), 0);

        dword_t shifter = 0;
        for (uint32_t i = words_to_shift, j = 0; j < ra.size(); i++, j++) {
            shifter |= (dword_t(ra[j]) << residue_bits);
            result[i] = shifter & BINARY_MASK;
            shifter >>= BINARY_SHIFT;
        }

        if (residue_bits)
            result[result.size()-1] = word_t(shifter);
        drop_leading_zeros(result);

        value_.clear();
        std::copy(result.begin(), result.end(), std::back_inserter(value_));

        return *this;
    }

private:
    template<typename T>
    void from_number(const T &rhs) {
        if (rhs == 0)
            return;

        if (rhs < 0) {
            std::cout << "Negative numbers are not supported\n";
            exit(1);
        }

        T abs = rhs;
        while (abs) {
            value_.push_back(abs & BINARY_MASK);
            abs >>= BINARY_SHIFT;
        }
    }

    /* Convert binary representation of BigInt to decimal base with weight of DECIMAL_BASE.
     * We can represent inplace value_ in DECIMAL_BASE. To support const BigInt, returning modified container.
     * We are not utilizing full bandwidth of uint32_t in decimal form for simplicity.
     *
     * Algorithm - Usual grade school methond of converting one base to another base.
    */
    value_container_t convert_to_decimal(void) const {
        int dsize = 1 + int(value_.size() * BINARY_SHIFT / (3 * DECIMAL_SHIFT));

        value_container_t decimal(dsize, 0);
        dsize = 0;
        for (int i = int(value_.size())-1; i >= 0; i--) {
            word_t hi_word = value_[i];

            for (int j = 0; j < dsize; j++) {
                dword_t x = (dword_t(decimal[j]) << BINARY_SHIFT) | hi_word;
                hi_word = word_t(x / DECIMAL_BASE);
                decimal[j] = word_t(x % DECIMAL_BASE); // % seems to be costly
            }

            while (hi_word) {
                decimal[dsize++] = hi_word % DECIMAL_BASE;
                hi_word /= DECIMAL_BASE;
            }
        }

        if (decimal.size() > 1)
            drop_leading_zeros(decimal);
        return std::move(decimal);
    }

    std::string to_sting() const {

        // Convert binary representation to decimal
        auto decimal(convert_to_decimal());

        // Print decimal representation
        assert(decimal.size() >= 1);
        std::string s(std::to_string(decimal[decimal.size()-1]));

        for (int i = int(decimal.size())-2; i >= 0; i--) {
            std::string number = std::to_string(decimal[i]);
            std::string zeros(DECIMAL_SHIFT-number.size(), '0');
            s += (zeros + number); // A number can be of [24 04 17] - 170424
        }

        return std::move(s);
    }



}; // Class definition

/* Generates product of all odd numbers in the range [A ... B], where A <= B + 1.
 * Forms a tree structure of evaluation, for e.g. B = 31 and A = 17, the recursion tree
 *
 *                         P(17, 31)
 *                     /     m=23      \
 *             P(17, 23)              P(25, 31)
 *            /  m=19  \             /  m=27  \
 *      P(17, 19)   P(21, 23)   P(25, 27)   P(29, 31)
 *       (17*19)     (21*23)     (25*27)     (29*31)
 *
 * i.e. we get product of 31,29,27,...19,17 as result.
 *
 * NOTE: The type of l and u should fit in word_t and their product should fit in dword_t.
 *       Since (l, u) <= N in N!, uint32_t is sufficient enough.
*/

BigInt partial_product(uint32_t l, uint32_t u) {
    if (u <= (l + 1))
        return BigInt(uint64_t(l));

    if (u == (l + 2))
        return BigInt(uint64_t(l) * uint64_t(u));

    // Find mid (left) odd element - we are spliting the range here.
    const static uint32_t no_branching[] = {1, 0};
    uint64_t m = l + ((u - l) >> 1);
    m = m - no_branching[m & 1];
    auto a = partial_product(l, m);
    auto b = partial_product(m+2, u);

    return a * b;
}

/* Recursively calculate all odd products, each call narrows the range.
   top most call - product of odd numbers in the range [N ... N/2)
   next call - product of odd numbers in the range [N/2 ... N/4)
   ... so on
   Until N > 3
*/
void odd_product_in_range(int n, BigInt &odds_product, BigInt &result) {
    if (n <= 2)
        return;

    odd_product_in_range(n/2, odds_product, result);
    odds_product = odds_product * partial_product(n / 2 + 1 + ((n / 2) & 1), n - 1 + (n & 1));
    result = result * odds_product;
}

/*  Algorithm credit to the following link

    http://www.luschny.de/math/factorial/binarysplitfact.html

    After few pages of note book work, the algorithm working became life long experience.

*/
BigInt Factorial(int32_t n) {
    if (n < 0) {
        std::cout << "Negative input to factorial\n";
        exit(1);
    }

    if (n == 0)
        return BigInt(1);

    BigInt p(1);
    BigInt r(1);
    odd_product_in_range(n, p, r);
    return r.left_shift(two_exponent(n));
}

} // math namespace

#endif // FAST_FACTORIAL_HPP

