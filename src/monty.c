#include <Python.h>
#include <montecarlo.h>
#include <models.h>

#define NSIMS 10000

/* prototypes */
int PyList_to_array_double(PyObject *list, double **array_double);
int PyList_to_array_int(PyObject *list, int **array_int);
int fill_data_matrix(PyObject *py_matrix, double ***matrix);
void fill_py_list_double(PyObject *list, double *array, int n);
int array_int_count(int *array, int len, int value);
void free_matrix(double **matrix, int nrows);

static PyObject *monty_wrapper(PyObject *self, PyObject *args)
{
  // name of the model
  char *modelname;
  struct model *model;
  // array with the "real" value of the parameters
  PyObject *params;
  double *cparams;
  int nparams;
  // array with the initial guess of the parameters
  PyObject *guess;
  double *cguess;
  // array indicating which parameters to fitd
  PyObject *fit;
  int *cfit;
  int nfit; // number of elements of cfit whic hare 1
  // matrix of independent variable(s) values
  PyObject *xi;
  double **cxi;
  int npoints; // number of data points
  // standard deviation of the gaussian error
  double dev;
  // number of simulations
  int nsims = NSIMS;

  /* OUTPUT from montecarlo: */
  int nsuccess;
  PyObject *py_nsuccess; // this will be used as part of the output

  /* Parse the arguments, on failure propagate exception */
  if (!PyArg_ParseTuple(args, "sOOOOdi", &modelname, &params, &guess, &fit, &xi,
                        &dev, &nsims))
    return NULL;

  /* Get the model. return with error if the name is unknown */
  if (!(model = get_model(modelname))) {
      return NULL;
  }
  /********** parameters **********/
  /* Put parameter values in cparams */
  nparams = PyList_to_array_double(params, &cparams);
  /* Check if the number of parameters is correct */
  if (nparams != model->nparams) {
    free(cparams);
    return NULL;
  }

  /********** guess of the parameters **********/
  /* Put guess values in cguess */
  nparams = PyList_to_array_double(guess, &cguess);
  /* Check if the number of parameters is correct */
  if (nparams != model->nparams) {
    free(cparams);
    free(cguess);
    return NULL;
  }

  /********** indicator of params to fit **********/
  nparams = PyList_to_array_int(fit, &cfit);
  /* Check if the number of parameters is correct */
  if (nparams != model->nparams) {
    free(cparams);
    free(cguess);
    free(cfit);
    return NULL;
  }
  nfit = array_int_count(cfit, model->nparams, 1);
  /* nfit must contain ONLY ones and zeros! */
  if ((nfit + array_int_count(cfit, model->nparams, 0)) != model->nparams) {
    free(cparams);
    free(cguess);
    free(cfit);
    return NULL;
  }

  /********* matrix of experimental points **********/
  if ((npoints = fill_data_matrix(xi, &cxi)) == -1) {
    free(cparams);
    free(cguess);
    free(cfit);
    free_matrix(cxi, npoints);
    return NULL;
  }

  /* now call montecarlo!!!! */
  printf("calling montecarlo with %d points and error %f\n", npoints, dev);
  double params_mean[model->nparams];
  double params_variance[model->nparams];
  nsuccess = montecarlo(model->function, cparams, cguess, dev, nsims, npoints,
                        model->nvars, model->nparams, nfit, cfit, cxi,
                        params_mean, params_variance, NULL);


  /* reuse params and guess to store the means and variances of the parameters,
   * respectively
   */
  fill_py_list_double(params, params_mean, model->nparams);
  fill_py_list_double(guess, params_variance, model->nparams);
  free(cparams);
  free(cguess);
  free(cfit);
  free_matrix(cxi, npoints);
  Py_DECREF(fit);
  py_nsuccess = Py_BuildValue("i", nsuccess);
  Py_INCREF(py_nsuccess);
  PyObject *out = PyList_New(3);
  out = Py_BuildValue("[OOO]", params, guess, py_nsuccess);
  Py_INCREF(params);
  Py_INCREF(guess);
  Py_INCREF(out);
  return out;
}

/* Take a PyObject (list) and an array of doubles, and place "n" doubles from 
 * the array into the list.
 */
void fill_py_list_double(PyObject *list, double *array, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    printf("%f\n", array[i]);
    PyList_SetItem(list, i, Py_BuildValue("f", array[i]));
  }
}

int array_int_count(int *array, int len, int item)
{
  int count = 0;
  int i;
  for (i = 0; i < len; i++) {
    if (array[i] == item)
      count++;
  }
  return count;
}

/* Take a list of python floats and place them as doubles in a dinamically
 * allocated array. Return the number of elements in the list.
 *
 * Note: the array of doubles must be manually freed!
 */
int PyList_to_array_double(PyObject *list, double **array_double)
{
  int i;
  int len = PyList_Size(list);
  PyObject *tmp_pyobject;
  *array_double = malloc(len * sizeof(double));
  for (i = 0; i < len; i++) {
    tmp_pyobject = PyList_GetItem(list, i);
    (*array_double)[i] = PyFloat_AsDouble(tmp_pyobject);
  }
  return len;
}

/* Take a list of lists of floats (list[nrows][ncols]), and place them in a C
 * matrix (an array of arrays) dinamically allocated. Return the number of
 * rows.
 *
 * Note: the matrix must be manually freed calling free_matrix(matrix, nrows)
 */
int fill_data_matrix(PyObject *py_matrix, double ***double_matrix)
{
  int i, j;
  int nrows = PyList_Size(py_matrix);
  int ncols = -1;
  int ncols_last = -1;
  PyObject *row;
  PyObject *point;
  // allocate memory for nrows pointers to double arrays
  *double_matrix = malloc(nrows * sizeof(double *));
  for (i = 0; i < nrows; i++) {
    // save the row
    row = PyList_GetItem(py_matrix, i);
    ncols = PyList_Size(row);
    /* raise Exception: all rows must have an equal number of columns */
    if ((ncols_last != -1) && (ncols != ncols_last)) {
      return -1;
    }
    // allocate space for ncols doubles
    (*double_matrix)[i] = malloc(ncols * sizeof(double));
    for (j = 0; j < ncols; j++) {
      point = PyList_GetItem(row, j);
      (*double_matrix)[i][j] = PyFloat_AsDouble(point);
      Py_DECREF(point);
    }
    Py_DECREF(row);
  }
  return nrows;
}

void free_matrix(double **matrix, int nrows)
{
  int i;
  for (i = 0; i < nrows; i++) {
    free(matrix[i]);
  }
  free(matrix);
}

int PyList_to_array_int(PyObject *list, int **array_int)
{
  int i;
  int len = PyList_Size(list);
  PyObject *tmp_pyobject;
  *array_int = malloc(len * sizeof(int));
  for (i = 0; i < len; i++) {
    tmp_pyobject = PyList_GetItem(list, i);
    (*array_int)[i] = (int) PyLong_AsLong(tmp_pyobject);
    Py_DECREF(tmp_pyobject);
  }
  return len;
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
