#include <util/TIMER.h>

#if _WIN32
bool TIMER::_firstCall = true;
double TIMER::_divisor;
#endif
