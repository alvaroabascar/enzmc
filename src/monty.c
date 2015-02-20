#include <Python.h>

static PyObject *monty_wrapper(PyObject *self, PyObject *args)
{
  char *modelname;   // name of the model
  PyObject *params;  // array with the "real" value of the parameters
  PyObject *guess;   // array with the initial guess of the parameters
  PyObject *fit;     // array indicating which parameters to fit
  PyObject *xi;      // matrix of independent variable(s) values
  double dev;        // standard deviation of the gaussian error
  int nsims;         // number of simulations
  if (!PyArg_ParseTuple(args, "sOOOOdi", &modelname, &params, &guess, &fit, &xi,
                        &dev, &nsims))
    return NULL;
  printf("string = %s\n", modelname);
  return Py_BuildValue("s", modelname);
}

static PyMethodDef montyMethods[] = {
  {"monty", monty_wrapper, METH_VARARGS, "Run simulations."},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef monty = {
  PyModuleDef_HEAD_INIT,
  "monty",
  NULL,
  -1,
  montyMethods
};

PyMODINIT_FUNC PyInit_monty(void)
{
  return PyModule_Create(&monty);
}
