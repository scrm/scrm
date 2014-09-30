#ifndef scrm_macros
#define scrm_macros

// Unless compiled with options NDEBUG, we will produce a debug output using 
// 'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG
#define dout std::cout
#else
#define dout 0 && std::cout
#endif

// from Knuths "The art of computer programming"
inline bool areSame(const double a, const double b, 
                    const double epsilon = std::numeric_limits<double>::epsilon()) {
  return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

#endif
