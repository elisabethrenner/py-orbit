#ifndef WRAP_ERROR_FUNCTIONS_H
#define WRAP_ERROR_FUNCTIONS_H


#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_errorfunctions{
    void initErrorfunctions(const char* errorfunctions_name);
  }
	
#ifdef __cplusplus
}
#endif // __cplusplus

#endif 
