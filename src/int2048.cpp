#include "int2048.h"

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
  // Long division is quadratic.  Use two stored limbs at once here: 10^8 is
  // still safe for every intermediate product in a signed long long, while
  // reducing the number of inner-loop iterations by almost a factor of four.
  const long long division_base = 100000000LL;
  std::vector<long long> u, v;
  for (int i = 0; i < static_cast<int>(left.digit.size()); i += 2)
    u.push_back(left.digit[i] + (i + 1 < static_cast<int>(left.digit.size()) ? static_cast<long long>(base) * left.digit[i + 1] : 0));
  for (int i = 0; i < static_cast<int>(right.digit.size()); i += 2)
    v.push_back(right.digit[i] + (i + 1 < static_cast<int>(right.digit.size()) ? static_cast<long long>(base) * right.digit[i + 1] : 0));
  const long long norm = division_base / (v.back() + 1);
  long long carry = 0;
  for (int i = 0; i < static_cast<int>(u.size()); ++i) {
    long long value = u[i] * norm + carry;
    u[i] = value % division_base; carry = value / division_base;
  }
  if (carry) u.push_back(carry);
  carry = 0;
  for (int i = 0; i < static_cast<int>(v.size()); ++i) {
    long long value = v[i] * norm + carry;
    v[i] = value % division_base; carry = value / division_base;
  }
  if (carry) v.push_back(carry);
  u.push_back(0);
  const int n = static_cast<int>(v.size());
  const int m = static_cast<int>(u.size()) - n - 1;
  std::vector<long long> q(m + 1, 0);
  for (int j = m; j >= 0; --j) {
    long long top = u[j + n] * division_base + u[j + n - 1];
    long long estimate = top / v[n - 1];
    long long rem = top % v[n - 1];
    if (estimate >= division_base) { estimate = division_base - 1; rem += v[n - 1]; }
    while (estimate * v[n - 2] > division_base * rem + u[j + n - 2]) { --estimate; rem += v[n - 1]; if (rem >= division_base) break; }
    long long product_carry = 0, borrow = 0;
    for (int i = 0; i < n; ++i) {
      long long product = estimate * v[i] + product_carry;
      product_carry = product / division_base;
      long long value = u[j + i] - (product % division_base) - borrow;
      if (value < 0) { value += division_base; borrow = 1; } else borrow = 0;
      u[j + i] = value;
    }
    long long value = u[j + n] - product_carry - borrow;
    if (value < 0) {
      --estimate;
      long long add_carry = 0;
      for (int i = 0; i < n; ++i) {
        long long sum = u[j + i] + v[i] + add_carry;
        if (sum >= division_base) { sum -= division_base; add_carry = 1; } else add_carry = 0;
        u[j + i] = sum;
      }
      u[j + n] += add_carry;
    } else u[j + n] = value;
    q[j] = estimate;
  }
  quotient.digit.clear();
  for (int i = 0; i < static_cast<int>(q.size()); ++i) {
    quotient.digit.push_back(static_cast<int>(q[i] % base));
    quotient.digit.push_back(static_cast<int>(q[i] / base));
  }
  quotient.negative = false; quotient.normalize();
  remainder.digit.assign(n * 2, 0);
  long long rem = 0;
  for (int i = n - 1; i >= 0; --i) {
    long long value = rem * division_base + u[i];
    long long denormalized = value / norm;
    rem = value % norm;
    remainder.digit[2 * i] = static_cast<int>(denormalized % base);
    remainder.digit[2 * i + 1] = static_cast<int>(denormalized / base);
  }
  remainder.negative = false; remainder.normalize();
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
