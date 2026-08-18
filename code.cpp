#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"

namespace sjtu {
class int2048 {
  static const int base = 10000;
  std::vector<int> digit;
  bool negative;

  void normalize();
  int compare_abs(const int2048 &) const;
  void add_abs(const int2048 &);
  void sub_abs(const int2048 &);
  void multiply_uint(int);
  static int2048 multiply_abs(const int2048 &, const int2048 &);
  static void divide_abs(const int2048 &, const int2048 &, int2048 &, int2048 &);
  static void divide_floor(const int2048 &, const int2048 &, int2048 &, int2048 &);
public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif

namespace sjtu {

namespace {

long long mod_pow(long long value, long long power, int mod) {
  long long result = 1;
  while (power) {
    if (power & 1) result = result * value % mod;
    value = value * value % mod;
    power >>= 1;
  }
  return result;
}

void ntt(std::vector<int> &a, bool invert, int mod, int primitive_root) {
  const int n = static_cast<int>(a.size());
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      int temp = a[i]; a[i] = a[j]; a[j] = temp;
    }
  }
  for (int length = 2; length <= n; length <<= 1) {
    long long root = mod_pow(primitive_root, (mod - 1) / length, mod);
    if (invert) root = mod_pow(root, mod - 2, mod);
    for (int start = 0; start < n; start += length) {
      long long w = 1;
      const int half = length >> 1;
      for (int j = 0; j < half; ++j) {
        int left = a[start + j];
        int right = static_cast<int>(a[start + j + half] * w % mod);
        int sum = left + right;
        if (sum >= mod) sum -= mod;
        int difference = left - right;
        if (difference < 0) difference += mod;
        a[start + j] = sum;
        a[start + j + half] = difference;
        w = w * root % mod;
      }
    }
  }
  if (invert) {
    const long long inverse_n = mod_pow(n, mod - 2, mod);
    for (int i = 0; i < n; ++i) a[i] = static_cast<int>(a[i] * inverse_n % mod);
  }
}

void convolution_mod(const std::vector<int> &left, const std::vector<int> &right,
                     std::vector<int> &out, int mod, int root) {
  int n = 1;
  const int needed = static_cast<int>(left.size() + right.size() - 1);
  while (n < needed) n <<= 1;
  std::vector<int> a(n, 0), b(n, 0);
  for (int i = 0; i < static_cast<int>(left.size()); ++i) a[i] = left[i];
  for (int i = 0; i < static_cast<int>(right.size()); ++i) b[i] = right[i];
  ntt(a, false, mod, root);
  ntt(b, false, mod, root);
  for (int i = 0; i < n; ++i) a[i] = static_cast<int>(static_cast<long long>(a[i]) * b[i] % mod);
  ntt(a, true, mod, root);
  out.swap(a);
}

} // namespace

int2048::int2048() : digit(1, 0), negative(false) {}

int2048::int2048(long long value) : negative(value < 0) {
  unsigned long long magnitude;
  if (value < 0) magnitude = static_cast<unsigned long long>(-(value + 1)) + 1;
  else magnitude = static_cast<unsigned long long>(value);
  do {
    digit.push_back(static_cast<int>(magnitude % base));
    magnitude /= base;
  } while (magnitude);
}

int2048::int2048(const std::string &text) : digit(1, 0), negative(false) { read(text); }

int2048::int2048(const int2048 &other) : digit(other.digit), negative(other.negative) {}

void int2048::normalize() {
  while (digit.size() > 1 && digit.back() == 0) digit.pop_back();
  if (digit.size() == 1 && digit[0] == 0) negative = false;
}

int int2048::compare_abs(const int2048 &other) const {
  if (digit.size() != other.digit.size()) return digit.size() < other.digit.size() ? -1 : 1;
  for (int i = static_cast<int>(digit.size()) - 1; i >= 0; --i) {
    if (digit[i] != other.digit[i]) return digit[i] < other.digit[i] ? -1 : 1;
  }
  return 0;
}

void int2048::add_abs(const int2048 &other) {
  const int old_size = static_cast<int>(digit.size());
  if (digit.size() < other.digit.size()) digit.resize(other.digit.size(), 0);
  int carry = 0;
  for (int i = 0; i < static_cast<int>(digit.size()); ++i) {
    int sum = digit[i] + (i < static_cast<int>(other.digit.size()) ? other.digit[i] : 0) + carry;
    if (sum >= base) { sum -= base; carry = 1; } else carry = 0;
    digit[i] = sum;
  }
  if (carry) digit.push_back(carry);
  (void)old_size;
}

void int2048::sub_abs(const int2048 &other) { // requires |*this| >= |other|
  int borrow = 0;
  for (int i = 0; i < static_cast<int>(digit.size()); ++i) {
    int value = digit[i] - (i < static_cast<int>(other.digit.size()) ? other.digit[i] : 0) - borrow;
    if (value < 0) { value += base; borrow = 1; } else borrow = 0;
    digit[i] = value;
  }
  normalize();
}

void int2048::multiply_uint(int factor) {
  if (factor == 0) { digit.assign(1, 0); negative = false; return; }
  long long carry = 0;
  for (int i = 0; i < static_cast<int>(digit.size()); ++i) {
    long long value = static_cast<long long>(digit[i]) * factor + carry;
    digit[i] = static_cast<int>(value % base);
    carry = value / base;
  }
  while (carry) { digit.push_back(static_cast<int>(carry % base)); carry /= base; }
}

void int2048::read(const std::string &text) {
  digit.assign(1, 0); negative = false;
  int start = 0;
  if (!text.empty() && (text[0] == '-' || text[0] == '+')) { negative = text[0] == '-'; start = 1; }
  digit.clear();
  for (int end = static_cast<int>(text.size()); end > start; end -= 4) {
    int begin = end - 4 < start ? start : end - 4;
    int value = 0;
    for (int i = begin; i < end; ++i) value = value * 10 + text[i] - '0';
    digit.push_back(value);
  }
  if (digit.empty()) digit.push_back(0);
  normalize();
}

void int2048::print() { std::cout << *this; }

int2048 &int2048::add(const int2048 &other) { return *this += other; }
int2048 add(int2048 left, const int2048 &right) { return left += right; }
int2048 &int2048::minus(const int2048 &other) { return *this -= other; }
int2048 minus(int2048 left, const int2048 &right) { return left -= right; }

int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const { int2048 result(*this); if (!(result.digit.size() == 1 && result.digit[0] == 0)) result.negative = !result.negative; return result; }
int2048 &int2048::operator=(const int2048 &other) { if (this != &other) { digit = other.digit; negative = other.negative; } return *this; }

int2048 &int2048::operator+=(const int2048 &other) {
  if (negative == other.negative) add_abs(other);
  else {
    int comparison = compare_abs(other);
    if (comparison >= 0) sub_abs(other);
    else { int2048 result(other); result.sub_abs(*this); *this = result; }
  }
  return *this;
}
int2048 operator+(int2048 left, const int2048 &right) { return left += right; }

int2048 &int2048::operator-=(const int2048 &other) { return *this += -other; }
int2048 operator-(int2048 left, const int2048 &right) { return left -= right; }

int2048 int2048::multiply_abs(const int2048 &left, const int2048 &right) {
  if ((left.digit.size() == 1 && left.digit[0] == 0) || (right.digit.size() == 1 && right.digit[0] == 0)) return int2048();
  int2048 result;
  const int n = static_cast<int>(left.digit.size()), m = static_cast<int>(right.digit.size());
  if (static_cast<long long>(n) * m <= 4096) {
    result.digit.assign(n + m, 0);
    for (int i = 0; i < n; ++i) {
      long long carry = 0;
      for (int j = 0; j < m || carry; ++j) {
        long long current = result.digit[i + j] + carry + (j < m ? static_cast<long long>(left.digit[i]) * right.digit[j] : 0);
        result.digit[i + j] = static_cast<int>(current % base);
        carry = current / base;
      }
    }
  } else {
    std::vector<int> first, second;
    convolution_mod(left.digit, right.digit, first, 998244353, 3);
    convolution_mod(left.digit, right.digit, second, 1004535809, 3);
    const long long inverse = mod_pow(998244353, 1004535809 - 2, 1004535809);
    result.digit.assign(n + m + 1, 0);
    long long carry = 0;
    for (int i = 0; i < n + m - 1 || carry; ++i) {
      long long coefficient = carry;
      if (i < n + m - 1) {
        long long delta = second[i] - first[i];
        if (delta < 0) delta += 1004535809;
        long long multiplier = delta * inverse % 1004535809;
        coefficient += first[i] + 998244353LL * multiplier;
      }
      result.digit[i] = static_cast<int>(coefficient % base);
      carry = coefficient / base;
    }
  }
  result.normalize();
  return result;
}

int2048 &int2048::operator*=(const int2048 &other) {
  const bool sign = negative != other.negative;
  *this = multiply_abs(*this, other);
  if (!(digit.size() == 1 && digit[0] == 0)) negative = sign;
  return *this;
}
int2048 operator*(int2048 left, const int2048 &right) { return left *= right; }

void int2048::divide_abs(const int2048 &left, const int2048 &right, int2048 &quotient, int2048 &remainder) {
  if (left.compare_abs(right) < 0) { quotient = int2048(); remainder = left; remainder.negative = false; return; }
  if (right.digit.size() == 1) {
    quotient = left; quotient.negative = false;
    long long rem = 0;
    for (int i = static_cast<int>(quotient.digit.size()) - 1; i >= 0; --i) {
      long long value = rem * base + quotient.digit[i];
      quotient.digit[i] = static_cast<int>(value / right.digit[0]);
      rem = value % right.digit[0];
    }
    quotient.normalize(); remainder = int2048(rem); return;
  }
  const int norm = base / (right.digit.back() + 1);
  int2048 u(left), v(right);
  u.negative = v.negative = false;
  u.multiply_uint(norm); v.multiply_uint(norm);
  u.digit.push_back(0);
  const int n = static_cast<int>(v.digit.size());
  const int m = static_cast<int>(u.digit.size()) - n - 1;
  quotient.digit.assign(m + 1, 0); quotient.negative = false;
  for (int j = m; j >= 0; --j) {
    long long top = static_cast<long long>(u.digit[j + n]) * base + u.digit[j + n - 1];
    long long estimate = top / v.digit[n - 1];
    long long rem = top % v.digit[n - 1];
    if (estimate >= base) { estimate = base - 1; rem += v.digit[n - 1]; }
    while (estimate * v.digit[n - 2] > static_cast<long long>(base) * rem + u.digit[j + n - 2]) { --estimate; rem += v.digit[n - 1]; if (rem >= base) break; }
    long long carry = 0, borrow = 0;
    for (int i = 0; i < n; ++i) {
      long long product = estimate * v.digit[i] + carry;
      carry = product / base;
      long long value = u.digit[j + i] - (product % base) - borrow;
      if (value < 0) { value += base; borrow = 1; } else borrow = 0;
      u.digit[j + i] = static_cast<int>(value);
    }
    long long value = u.digit[j + n] - carry - borrow;
    if (value < 0) {
      --estimate;
      long long add_carry = 0;
      for (int i = 0; i < n; ++i) {
        long long sum = u.digit[j + i] + v.digit[i] + add_carry;
        if (sum >= base) { sum -= base; add_carry = 1; } else add_carry = 0;
        u.digit[j + i] = static_cast<int>(sum);
      }
      u.digit[j + n] = static_cast<int>(u.digit[j + n] + add_carry);
    } else u.digit[j + n] = static_cast<int>(value);
    quotient.digit[j] = static_cast<int>(estimate);
  }
  quotient.normalize();
  remainder.digit.assign(u.digit.begin(), u.digit.begin() + n);
  remainder.negative = false;
  long long rem = 0;
  for (int i = n - 1; i >= 0; --i) {
    long long value = rem * base + remainder.digit[i];
    remainder.digit[i] = static_cast<int>(value / norm);
    rem = value % norm;
  }
  remainder.normalize();
}

void int2048::divide_floor(const int2048 &left, const int2048 &right, int2048 &quotient, int2048 &remainder) {
  divide_abs(left, right, quotient, remainder);
  const bool opposite = left.negative != right.negative;
  const bool has_remainder = !(remainder.digit.size() == 1 && remainder.digit[0] == 0);
  if (opposite && has_remainder) {
    quotient.add_abs(int2048(1));
    remainder = [&]() { int2048 result(right); result.negative = false; result.sub_abs(remainder); return result; }();
  }
  if (!(quotient.digit.size() == 1 && quotient.digit[0] == 0)) quotient.negative = opposite;
  if (has_remainder) remainder.negative = opposite ? right.negative : left.negative;
  else remainder.negative = false;
}

int2048 &int2048::operator/=(const int2048 &other) { int2048 quotient, remainder; divide_floor(*this, other, quotient, remainder); *this = quotient; return *this; }
int2048 operator/(int2048 left, const int2048 &right) { return left /= right; }
int2048 &int2048::operator%=(const int2048 &other) { int2048 quotient, remainder; divide_floor(*this, other, quotient, remainder); *this = remainder; return *this; }
int2048 operator%(int2048 left, const int2048 &right) { return left %= right; }

std::istream &operator>>(std::istream &input, int2048 &value) { std::string text; input >> text; if (input) value.read(text); return input; }
std::ostream &operator<<(std::ostream &output, const int2048 &value) {
  if (value.negative) output << '-';
  output << value.digit.back();
  for (int i = static_cast<int>(value.digit.size()) - 2; i >= 0; --i) {
    int x = value.digit[i];
    if (x < 1000) output << '0';
    if (x < 100) output << '0';
    if (x < 10) output << '0';
    output << x;
  }
  return output;
}

bool operator==(const int2048 &left, const int2048 &right) { return left.negative == right.negative && left.digit == right.digit; }
bool operator!=(const int2048 &left, const int2048 &right) { return !(left == right); }
bool operator<(const int2048 &left, const int2048 &right) { if (left.negative != right.negative) return left.negative; int comparison = left.compare_abs(right); return left.negative ? comparison > 0 : comparison < 0; }
bool operator>(const int2048 &left, const int2048 &right) { return right < left; }
bool operator<=(const int2048 &left, const int2048 &right) { return !(right < left); }
bool operator>=(const int2048 &left, const int2048 &right) { return !(left < right); }

} // namespace sjtu
