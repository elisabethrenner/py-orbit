#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_utils.hh"

#include <iostream>
#include <string>
#include <cmath>
#include <cfloat>
#include <complex>

#include "ErrorFunctions.hh"

using namespace OrbitUtils;
using namespace wrap_orbit_utils;

namespace wrap_errorfunctions{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	static PyObject* errorfunctions_cerrf(PyObject* self, PyObject* args)
	{
		std::complex<double> z;
		if(!PyArg_ParseTuple(args,"D:cerrf",&z)){
			error("cerrf(x) - parameter is needed.");
		}
		std::complex<double> w = cerrf(z);
		return Py_BuildValue("D", &w);	
	}

	
	static PyMethodDef ErrorfunctionsModuleMethods[] = { 
		{"cerrf", errorfunctions_cerrf , METH_VARARGS, "cerrf(x) - calculated the complex error function."},
		{NULL}
	};
  
	//--------------------------------------------------
	//Initialization functions of the errorfunctions module
	//--------------------------------------------------
	void initErrorfunctions(const char* errorfunctions_name){
	//create errorfunctions module
	PyObject* module_ef = Py_InitModule(errorfunctions_name,ErrorfunctionsModuleMethods);			
	}

#ifdef __cplusplus
}
#endif

}
