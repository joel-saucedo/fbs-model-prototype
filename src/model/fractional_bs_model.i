%module fractional_bs_model

%{
// Include necessary headers
#include "fractional_bs_model.h"
#include <vector>
#include <string>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>  // NumPy C API header

// Function to initialize NumPy's C API
int init_numpy() {
    if (_import_array() < 0) {
        PyErr_Print();  // Print the error to help diagnose the issue
        return -1;  // Indicate failure
    }
    return 0;  // Indicate success
}
%}

%include <std_vector.i>  // SWIG support for std::vector
%include <std_string.i>  // SWIG support for std::string
%include "numpy.i"// SWIG support for NumPy arrays


// Enable std::vector handling for Python
%template(VectorDouble) std::vector<double>;
%template(VectorString) std::vector<std::string>;

// Expose functions to Python
%include "fractional_bs_model.h"

// Wrap the fractional Black-Scholes PDE function
%inline %{
std::vector<double> fractional_black_scholes_alpha(const std::vector<double>& S, double alpha, double sigma, double r,
                                                   double dx, int N, int L, double K, double T, const std::string& option_type) {
    return ::fractional_black_scholes_alpha(S, alpha, sigma, r, dx, N, L, K, T, option_type);
}

// Wrap the calibration objective function
double calibration_objective(double alpha, const std::vector<std::vector<double>>& options_data,
                             const std::string& option_type, const std::vector<double>& S, double sigma,
                             double r, double dx, int N, int L, double S0) {
    return ::calibration_objective(alpha, options_data, option_type, S, sigma, r, dx, N, L, S0);
}
%}

// Add NumPy support for std::vector<double>
%typemap(in) std::vector<double> (PyObject* py_obj) {
    if (!PyArray_Check(py_obj)) {
        PyErr_SetString(PyExc_TypeError, "Expected a numpy array");
        return nullptr;
    }
    PyArrayObject* array = (PyArrayObject*)py_obj;
    if (PyArray_TYPE(array) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "Expected a numpy array of type double");
        return nullptr;
    }
    int size = PyArray_SIZE(array);
    double* data = (double*)PyArray_DATA(array);
    $1 = std::vector<double>(data, data + size);
}

%typemap(out) std::vector<double> {
    npy_intp dims[] = { (npy_intp)$1.size() };
    PyObject* result = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (!result) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to create numpy array");
        SWIG_fail;
    }
    memcpy(PyArray_DATA((PyArrayObject*)result), $1.data(), $1.size() * sizeof(double));
    $result = result;
}
