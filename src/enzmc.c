/* Coolest CLI ever */
#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <regex.h>
#include "include/montecarlo.h"
#include "include/enzyme.h"
#include "include/enzmc.h"
#include "include/models.h"

#include "models/models.c"

/* MAX_INDEP, MAX_DATA_CHARS, MAX_PARAMS in models.h */
#define NREPS 10000

const char *argp_program_bug_address = "alvaroabascar@gmail.com";
const char *argp_program_version = "version 1.0";

static int parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *args = state->input;
    static nargs = 0;
    switch(key) {
        case 'f':
            args->mode = FILE_MODE;
            args->fileinput = arg;
            break;
        case 't':
            args->mode = TEMPLATE_MODE;
            args->model = arg;
            break;
        case 'v':
            args->verbose = 1;
            break;
        case 444: /* model */
            args->model = arg;
            break;
        case 555: /* parameters */
            args->params = arg;
            break;
        case 666: /* error */
            args->error = arg;
            break;
        case 777: /* data points */
            args->data = arg;
            break;
        case 999: /* parameters to fix */
            args->fixed_params = arg;
            break;
        case ARGP_KEY_ARG:
            args->fileoutput = arg;
            nargs++;
        case ARGP_KEY_END:
            if (args->mode == NORMAL_MODE) {
              if (!(args->model && args->params && args->error &&
                    args->data))
                argp_failure(state, 1, 0, ERROR_LACK_OPTS);
            }
            if (nargs > 1) {
              argp_failure(state, 1, 0, ERROR_TOO_MANY_ARGS);
            }
            if (args->mode == TEMPLATE_MODE)
              if (nargs < 1)
                argp_failure(state, 1, 0, ERROR_NO_FILENAME);
            break;
    }
    return 0;
}

int main(int argc, char *argv[])
{
  int result;
  struct argp_option options[] = {
    {0, 0, 0, 0, "Other program modes:", 1},
    {"template", 't', "model", 0, "Create a template for the specified model"},
    {"file", 'f', "input file", 0, "Get input from file"},
    {0, 0, 0, 0, "Mandatory parameters:", 2},
    {"model", 444, "\"model name\"", 0, "Choose a model"},
    {"params", 555, "\"Vm=5 Km=17\"", 0, "Set the parameters of the model"},
    {"error", 666, "error", 0, "Set the estimated measurement error (absolute value)"},
    {"data", 777, "\"X1=[0.3,0.45,0.6,...] X2=[0.1,0.2,0.4...]\"", 0, "Set the values of the independent variables"},
    {0, 0, 0, 0, "Optional parameters:", 3},
    {"fixed", 999, "\"Vm Kd\"", 0, "Indicate what parameters are fixed"},
    {0, 0, 0, 0, "Informational options:", -1},
    {"verbose", 'v', 0, 0, "Display arguments on output"},
    {0}};
  struct argp argp = {options, parse_opt, "OUTPUT_FILE"};
  struct arguments args = {
      .mode = NORMAL_MODE,
      .model = NULL,
      .params = NULL,
      .fixed_params = NULL,
      .error = NULL,
      .data = NULL,
      .fileoutput = NULL,
      .fileinput = NULL,
      .verbose = 0
    };
    argp_parse(&argp, argc, argv, 0, 0, &args);

    /********************************************
     * create a simple API to contain the modes?
     * *******************************************/
    switch(args.mode) {
    case NORMAL_MODE:
        result = run_cli_mode(&args);
        break;
    case FILE_MODE:
        result = run_file_mode(&args);
        break;
    case TEMPLATE_MODE:
        result = create_template(args.model, args.fileoutput);
        break;
    }
    return result;
}

/* Given the user input (as a set of arguments stored in a
 * struct arguments, call montecarlo to do the simulations and
 * print the output to the user
 */
int run_cli_mode(struct arguments *args)
{
  int ret_code; // return code of each of the functions called */
  /* Retrieve the data related to this model, necessary to parse the input */
  struct model *model = get_model(args->model);
  if (!model) {
    fprintf(stderr, "Unrecognized model name: %s\n", args->model);
    return -1;
  }
  /* Array to store the arrays of values for the independent variables. nvars
   * is the number of independent variables of the model.
   */
  double *data[model->nvars];
  /* Array to hold the values of the parameters. */
  double params[model->nparams];
  /* Array to indicate the fixed/adjustable parameters. */
  int fixed_params[model->nparams];
  int *fixed_params_ptr = fixed_params;
  /* Standard deviation of the gaussian error. */
  double error;
  /* Number of points (num of values of the each variable). */
  int npoints;
  /* Number of parameters to fit */
  int mfit;
  /* Means and variacnes of the adjusted parameters */
  double means[model->nparams];
  double variances[model->nparams];

  /* parse data passed in the --data option (values of the indep variables), and
   * put it into the matrix "data"
   */
  if (ret_code = get_indep_vars(model, args->data, data)) {
    return ret_code;
  }
  /* parse parameters and place them in the array params */
  if (ret_code = get_params(model, args->params, params)) {
    return ret_code;
  }
  /* parse error term */
  if (ret_code = get_error(model, args->error, &error)) {
    return ret_code;
  }
  /* parse the fixed parameters */
  /* get_fixed_params(model, args.fixed_params, fixed_params_ptr); */
  printf("%s\n", model->name);
  /*
  montecarlo(model.function, params, params, error, NREPS, npoints,
             model->nvars, model->nparams, mfit, fixed_params_ptr, data,
             means, variances, NULL);
   */
  free_array_double(data, model->nvars);
  return -1;
}

int run_file_mode(struct arguments *args)
{

}

int create_template(char *modelname, char *fileout)
{

}

/* Given the name of a model, return a ptr to the struct model with all the
 * information necessary to use the model (see models.h).
 */
struct model *get_model(char *modelname)
{
  int i;
  /* models is an array of struct model, defined in models.c
   * end of models is defined by a struct model with all its elements set to
   * zero / NULL
   */
  for (i = 0; models[i].function != NULL; i++) {
    if (!strcmp(models[i].name, modelname)) {
      return &(models[i]);
    }
  }
  return NULL;
}

/* Given a struct model and raw data as a string, parse this string and place
 * the values of the independent variables in a matrix. A reference to this
 * matrix must be passed as third element.
 * 
 * eg. of raw data: a string with the format   "S=[1,2,3,4,5,6] I=[2,2,4,4,5,5]"
 * eg. of output:   array of arrays type double { {1,2,3,4,5,6}, {2,2,4,4,5,5} }
 *
 * returns: 0 on success
 *         -1 on failure (not all indep variables of the model have been found
 *                        in the data)
 */
int get_indep_vars(struct model *model, char *raw_data,
                   double *data[model->nvars])
{
  int i;
  /* simple regexp which will match a string of type:
   * foo = [1, 2, 3, 4, 5, 6, 7], accepting commas, spaces and numbers either
   * in decimal or exponential notation
   */
  char *regexp_base = "%s[ ]*=[ ]*[[eE., [:digit:]-]+]";
  char regexp_complete[40];
  char match_str[200]; // will store the matched string
  /* for each independent variable in the model... */
  for (i = 0; i < model->nvars; i++) {
    /* build the full regex to match the ith indep. variable of the model */
    sprintf(regexp_complete, regexp_base, model->indep_vars[i]);
    if (extract_str(raw_data, match_str, regexp_complete)) {
      fprintf(stderr, "Error: variable %s is required by the model %s, but it was not found.\n", model->indep_vars[i], model->name);
      return -1;
    }
    data[i] = parse_array_double(match_str);
  }
  return 0;
}

/* Given a struct model and raw data as a string, parse this string to find the
 * values of the parameters of the model. An array of type double must be passed
 * as third argument, in which the values of the parameters will be placed.
 * 
 * eg. of raw data: a string with the format   "Vm=[2,3,4,5,6] Km=[2,4,4,5,5]"
 * eg. of output:   array of arrays type double { {2,3,4,5,6}, {2,4,4,5,5} }
 *
 * returns: 0 on success
 *         -1 on failure (not all parameters of the model have been found
 *                        in the data)
 */
int get_params(struct model *model, char *raw_data, double *params)
{
  /* simple regexp which will match a string of type foo = [1,3,4,5,6],
   * accepting commas or spaces as separators, and numbers either in decimal or
   * exponential notation
   */
  char *regexp_base = "%s[ ]*=[ ]*[eE.[:digit:]-]+";
  char regexp_complete[40];
  char match_str[200];
  int i;
  /* for each parameter in the model... */
  for (i = 0; i < model->nparams; i++) {
    /* build the full regex to match the ith parameter of the model */
    sprintf(regexp_complete, regexp_base, model->params[i]);
    if (extract_str(raw_data, match_str, regexp_complete)) {
      fprintf(stderr, "Error: parameter %s is required by the model %s, but it was not found.\n", model->params[i], model->name);
      return -1;
    }
    /* place the number in the correct place of the array params */
    sscanf(match_str, "%*[^0-9]%lf", &params[i]);
  }
  return 0;
}

/* Given a model, and a string containing the error, find the value
 * and store it in error (double)
 */
int get_error(struct model *model, char *raw_data, double *error)
{
  int d = sscanf(raw_data, "%lf", error);
  if (d != 1) {
    return -1;
  }
  return 0;
}

/* Given a model, and a string containing a list of parameters, find them and
 * fill the array pointed to by ptr in the following way:
 *
 * *ptr is an int array with the length of model->nparams (one entry for each
 * parameter). If a given entry is 1, the parameter will be fitted; if it is 0,
 * it will be kept fixed.
 */
/* int get_fixed_params(struct model *model, char *raw_data, int *ptr[]) */
/* { */
/*   char str[10]; // useless but required by extract_str */
/*   int i; */
/*   #<{(| for each parameter of the model... |)}># */
/*   for (i = 0; i < model->nparams; i++) { */
/*     #<{(| can we match params[i]? if so, then params[i] is in the string raw_data, */
/*      * which means that it should be fixed */
/*      |)}># */
/*     if (extract_str(raw_data, str, model->params[i])) */
/*       ptr[i] = 0; */
/*     else */
/*       ptr[i] = 1; */
/*   } */
/*  */
/* } */

/* Given a source and destiny string, and a regular expression, match the
 * regular expression in the source string and save the matched region in the
 * destiny string.
 */
int extract_str(char *src, char *dst, char *regexp_str)
{
  size_t nchars; // to be used later
  regex_t regex;
  regmatch_t pmatch[1];
  regcomp(&regex, regexp_str, REG_EXTENDED);
  /* if didn't match... */
  if (regexec(&regex, src, 1, pmatch, 0)) {
    return -1;
  }
  /* copy the matched string to the destiny string. pmatch[0].rm_so contains the
   * starting position of the match, pmatch[0].rm_eo the ending position. So we
   * copy (rm_eo - rm_so) chars from the string starting at rm_so.
   */
  nchars = pmatch[0].rm_eo - pmatch[0].rm_so;
  strncpy(dst, &src[pmatch[0].rm_so], nchars);
  dst[nchars] = '\0'; // end string!
  return 0;
}

/* Given a string in format [1,2,3,3,5.6, 3.2], return a pointer to an array
 * with this numbers as doubles. Numbers might be separated by commas or spaces,
 * and exponential notation is accepted as well as decimal notation.
 *
 * The array might be preceded or followed by something else (eg. "A = [1,2,3]")
 * In this case the function will consider as array the first match of a set of
 * numbers between brackets.
 */
double *parse_array_double(char *array_raw)
{
  int i;
  char *token;
  double *result;
  /* set result_tmp to maximum possible size */
  double *result_tmp = malloc(MAX_VALUES * sizeof(double));
  /* first step: clean the array, leave only the digits + brackets */
  char array_clean[100];
  if (extract_str(array_raw, array_clean, "\[[eE., [:digit:]-]+]")) {
    return NULL;
  }
  /* part 2: turn it into double array */
  /* tokenize by commas, spaces and brackets (we remove brackets this way) */
  token = strtok(array_clean, ", []");
  for (i = 0; token; i++) {
    result_tmp[i] = atof(token);
    token = strtok(NULL, ", []");
  }

  /* now copy the array to result, using the minimum required space */
  result = malloc(i * sizeof(double));
  copy_array_double(result, result_tmp, i);
  /* now free result_tmp, which was (probably) excesively large */
  free(result_tmp);
  return result;
}

void free_array_double(double *array[], int n)
{
  int i;
  for (i = 0; i < n; i++) {
    free(array[i]);
  }
}

/* Copy n elements from src to dst */
void copy_array_double(double *dst, double *src, int n)
{
  while (--n >= 0)
    dst[n] = src[n];
}
