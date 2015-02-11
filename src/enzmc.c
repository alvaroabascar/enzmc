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
        case 888: /* guess of the parameters */
            args->guess = arg;
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
    {"guess", 888, "\"Vm=4 Km=5\"", OPTION_HIDDEN, "Set the guess of the parameters"},
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
      .guess = NULL,
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

  printf("%s\n", model->name);
  /*
  montecarlo(model.function, params, params, error, NREPS, npoints, model.nvars,
             model.nparams, mfit, fixed_params_ptr, data, means, variances,
             NULL);
   */
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
