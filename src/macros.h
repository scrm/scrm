#ifndef scrm_macros
#define scrm_macros

// Unless compiled with options NDEBUG, we will produce a debug output using 
// 'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG
#define dout std::cout
#else
#define dout 0 && std::cout
#endif

#endif
