#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include <math.h>
#include <regex.h>
#include "include/montecarlo.h"
#include "include/enzyme.h"

/* modes */
#define NORMAL_MODE 52
#define INTERACTIVE_MODE 53
#define GRAPHICAL_MODE 54
#define FILE_MODE 55
#define TEMPLATE_MODE 56

/* maximum number of independent variables */
#define MAX_INDEP 10
/* maximum number of chars in the --data field */
#define MAX_DATA_CHARS 1000
/* maximum number of parameters */
#define MAX_PARAMS 10

#define NREPS 10000
const char *_models_[20] = {
    "michaelis",
    "alberty",
    "pingpong",
    "mixed",
    "competitive",
    "uncompetitive",
    "noncompetitive",
    "ph",
    "michaelistemp",
    "michaelisinactiv"};

const char *argp_program_bug_address = "alvaroabascar@gmail.com";
const char *argp_program_version = "version 0.5";

struct arguments {
    short mode;
    short verbose;
    char *model;
    char *params;
    char *fixed_params;
    char *guess;
    char *data;
    char *error;
    char *fileoutput;
    char *fileinput;
};

int run_interactive_mode(struct arguments *);
int run_graphical_mode(struct arguments *);
int run_normal_mode(struct arguments *);
int parse_array(double *S_arr, int n, char S_str[]);
int strcopy(char *dst, char *src, int start, int end);

static int parse_opt(int key, char *arg, struct argp_state *state)
{
    struct arguments *args = state->input;
    static nargs = 0;
    switch(key) {
        case 'I':
            args->mode = INTERACTIVE_MODE;
            break;
        case 'G':
            args->mode = GRAPHICAL_MODE;
            break;
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
                    argp_failure(state, 1, 0, "lacking options (model, params, data and error are mandatory)");
            }
            if (nargs > 1) {
                argp_failure(state, 1, 0, "too many arguments (see --usage)");
            }
            if (args->mode == TEMPLATE_MODE)
                if (nargs < 1)
                    argp_failure(state, 1, 0, "you must specify a file name:\n./enzmc --template model filename");
            break;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    int result;
    struct argp_option options[] = {
        {0, 0, 0, 0, "Other program modes:", 1},
        {"interactive", 'I', 0, 0, "Enter interactive mode"},
        {"template", 't', "model", 0, "Create a template for the specified model"},
        {"file", 'f', "input file", 0, "Get input from file"},
        {"graphical", 'G', 0, OPTION_HIDDEN, "Enter graphical mode"},
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
    struct arguments args;
    args.mode = NORMAL_MODE;
    args.model = NULL;
    args.params = NULL;
    args.fixed_params = NULL;
    args.guess = NULL;
    args.error = NULL;
    args.data = NULL;
    args.fileoutput = NULL;
    args.fileinput = NULL;
    args.verbose = 0;
    argp_parse(&argp, argc, argv, 0, 0, &args);
    switch(args.mode) {
    case INTERACTIVE_MODE:
        result = run_interactive_mode(&args);
        break;
    case GRAPHICAL_MODE:
        result = run_graphical_mode(&args);
        break;
    case NORMAL_MODE:
        result = run_normal_mode(&args);
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

int create_template(char *model, char *filename)
{
    FILE *fp;
    char *params, *indep_vars;
    if (!(fp = fopen(filename, "w"))) {
        fprintf(stderr, "Error: could not open %s.\n", filename);
        return -1;
    }
    /* michaelis */
    if (strcmp(model, _models_[0]) == 0) {
        params = "Vmax=\nKm=";
        indep_vars = "S=[  ]";
    /* alberty */
    } else if (strcmp(model, _models_[1]) == 0) {
        params = "Vmax=\nKmA=\nKmB=\nKsA=\n";
        indep_vars = "A=[  ]\nB=[  ]";
    /* pingpong */
    } else if (strcmp(model, _models_[2]) == 0) {
        params = "Vmax=\nKmA=\nKmB=";
        indep_vars = "A=[  ]\nB=[  ]";
    /* mixed */
    } else if (strcmp(model, _models_[3]) == 0) {
        params = "Vmax=\nKm=\nKIa=\nKIb=";
        indep_vars = "S=[  ]\nI=[  ]";
    /* competitive */
    } else if (strcmp(model, _models_[4]) == 0) {
        params = "Vmax=  Km=  KIa=";
        indep_vars = "S=[  ]\nI=[  ]";
    /* uncompetitive */
    } else if (strcmp(model, _models_[5]) == 0) {
        params = "Vmax=\nKm=\nKIb=";
        indep_vars = "S=[  ]\nI=[  ]";
    /* noncompetitive */
    } else if (strcmp(model, _models_[6]) == 0) {
        params = "Vmax=\nKm=\nKIb=";
        indep_vars = "S=[  ]\nI=[  ]";
    /* ph */
    } else if (strcmp(model, _models_[7]) == 0) {
        params = "Vmax=\nKm=\nKa1=\nKa2=\nKa3=\nKa4=";
        indep_vars = "S=[  ]\nH=[  ]";
    /* michaelis + temperature */
    } else if (strcmp(model, _models_[8]) == 0) {
        params = "Vmax=\nEa=\nKm=\nT1=\n";
        indep_vars = "S=[ ]\nT=[ ]";
    /* michaelis + thermal inactivation */
    } else if (strcmp(model, _models_[9]) == 0) {
        params = "Vmax=\nKm=\nkt=\n";
        indep_vars = "S=[ ]\nt=[ ]";
    } else {
        fprintf(stderr, "Error: model \"%s\" is not recognized.\n", model);
        return -1;
    }

    fprintf(fp, "- Model: %s\n\n", model);
    fprintf(fp, "- Independent variables:\n\n%s\n\n", indep_vars);
    fprintf(fp, "- Parameters:\n\n%s\n\n", params);
    fprintf(fp, "- Fixed parameters: None\n\n");
    fprintf(fp, "- Error (absolute):");
}


int run_interactive_mode(struct arguments *args)
{
    char tmp[10];
    struct arguments newargs;
    newargs.model = malloc(100 * sizeof(char));
    newargs.data = malloc(MAX_DATA_CHARS * sizeof(char));
    newargs.params = malloc(100 * sizeof(char));
    newargs.guess = malloc(100 * sizeof(char));
    newargs.error = malloc(100 * sizeof(char));
    newargs.fixed_params = malloc(100 * sizeof(char));
    if (args->model == NULL) {
        printf("Enter the name of the model (e.g. \"michaelis\")\n>> ");
        fgets(newargs.model, 100, stdin);
    } else {
        newargs.model = args->model;
    }
    if (args->data == NULL) {
        printf("Enter the independent variables (e.g. S=[1,2,3] I=[1.5,2,4])\n>> ");
        fgets(newargs.data, MAX_DATA_CHARS, stdin);
    } else {
        newargs.data = args->data;
    }

    if (args->params == NULL) {
        printf("Enter the parameters (e.g. Km=5.5 Vmax=32.1)\n>> ");
        fgets(newargs.params, 100, stdin);
    } else {
        newargs.params = args->params;
    }

    if (args->error == NULL) {
        printf("Enter an estimation of the error (absolute value)\n>> ");
        fgets(newargs.error, 100, stdin);
    } else {
        newargs.error = args->error;
    }

    if (args->guess == NULL) {
        printf("Do you want to use the value of the parameters as the initial guess (P), or to provide a new guess (N)?\n>> [P/n] ");
        fgets(tmp, 10, stdin);
        if (!(strcmp(tmp, "n\n") || strcmp(tmp, "N\n"))) {
            printf("Enter the guess (e.g. Km=5 Vmax=20\n)>> ");
            fgets(newargs.guess, 100, stdin);
        } else {
            newargs.guess = newargs.params;
        }
    } else {
        newargs.guess = args->guess;
    }

    if (args->fixed_params == NULL) {
        printf("Enter the parameters that you wish to keep fixed, or nothing if you want to fit all of them\n>> ");
        fgets(newargs.fixed_params, 100, stdin);
    } else {
        newargs.fixed_params = args->fixed_params;
    }
    return run_normal_mode(&newargs);
}

int run_graphical_mode(struct arguments *args)
{
    printf("Oops! I hope this was already available\n");
    return -1;
}

int run_file_mode(struct arguments *args)
{
    char model[20], data[MAX_DATA_CHARS], params[100], fixed[100], error[10];
    char *line = malloc(sizeof(char) * MAX_DATA_CHARS);
    char *strtmp = malloc(sizeof(char) * MAX_DATA_CHARS);
    size_t nbytes = MAX_DATA_CHARS;
    regex_t regex;
    regmatch_t pmatch[1];
    FILE *fp = fopen(args->fileinput, "r");

    if (fp == NULL) {
        fprintf(stderr, "Error: could not open input file %s.\n", args->fileinput);
        return -1;
    }

    /* model name */
    getline(&line, &nbytes, fp);
    regcomp(&regex, "[a-zA-Z]*\n", REG_EXTENDED);
    if (regexec(&regex, line, 1, pmatch, 0)) {
        fprintf(stderr, "Error: Could not find model name.\n");
        exit(1);
    }
    strcopy(model, line, pmatch[0].rm_so, pmatch[0].rm_eo - 1);
    strtmp[0] = '\0';
    /* independent variables */
    /* jump to location of independent variables */
    do {
        getline(&line, &nbytes, fp);
    } while (strcmp(line, "\n") == 0 || strcmp(line, "- Independent variables:\n") == 0);

    regcomp(&regex, "[a-zA-Z]*[ ]*=[ ]*[[0-9eE., -]*]\n", REG_EXTENDED);
    while (!regexec(&regex, line, 1, pmatch, 0)) {
        sprintf(data, "%s%s", strtmp, line);
        strcpy(strtmp, data);
        if (!getline(&line, &nbytes, fp))
            break;
    }
    /* parameters */
    /* jump to parameters */
    do {
        getline(&line, &nbytes, fp);
    } while (strcmp(line, "\n") == 0 || strcmp(line, "- Parameters:\n") == 0);

    regcomp(&regex, "[a-zA-Z0-9]*[ ]*=[ ]*[0-9eE.-]*\n", REG_EXTENDED);
    strtmp[0] = '\0';
    while (!regexec(&regex, line, 1, pmatch, 0)) {
        sprintf(params, "%s%s", strtmp, line);
        strcpy(strtmp, params);
        if (!getline(&line, &nbytes, fp))
            break;
    }
    /* fixed parameters */
    do {
        getline(&line, &nbytes, fp);
    } while (strcmp(line, "\n") == 0);
    
    regcomp(&regex, "[0-9a-zA-Z ]*\n", REG_EXTENDED);
    regexec(&regex, line, 1, pmatch, 0);
    strcopy(fixed, line, pmatch[0].rm_so, pmatch[0].rm_eo - 1);
    /* error */
    do {
        getline(&line, &nbytes, fp);
    } while (strcmp(line, "\n") == 0);
    getline(&line, &nbytes, fp);

    regcomp(&regex, "[0-9eE. -]*\n", REG_EXTENDED);
    regexec(&regex, line, 1, pmatch, 0);
    strcopy(error, line, pmatch[0].rm_so, pmatch[0].rm_eo - 1);
    free(line);
    free(strtmp);
    args->model = model;
    args->data = data;
    args->fixed_params = (strcmp(fixed, "") == 0) ? NULL : fixed;
    args->error = error;
    args->params = params;
    return run_normal_mode(args);
}

int run_normal_mode(struct arguments *args)
{
    int i, j, m, mfit, n, nvars, *fitparams;
    char *data, *data2;
    char *params_name[MAX_PARAMS];
    double error;
    double (*modelfunc)(double X[], double p[]);
    double *S[MAX_INDEP];
    double *params, *guess;
    FILE *fp = fopen(args->fileoutput, "w");
    params = malloc(MAX_PARAMS * sizeof(double));
    guess = malloc(MAX_PARAMS * sizeof(double));
    fitparams = malloc(MAX_PARAMS * sizeof(double));
    for (i = 0; i < MAX_INDEP; i++)
        S[i] = NULL;
    if (args->verbose) {
        printf("Model: %s\n", args->model);
        printf("Parameters: %s\n", args->params);
        printf("Fixed parameters: %s\n", args->fixed_params);
        printf("Guess: %s\n", args->guess ? args->guess : args->params);
        printf("Error: %s\n", args->error);
        printf("Data: %s\n", args->data);
        printf("\n");
    } else {
        printf("Data: %s\n", args->data);
    }

    if ((error = atof(args->error)) == 0) {
        fprintf(stderr, "Error: \"error\" must be a numeric value > 0\n");
        return -1;
    }
    /* S will hold the values of all independent variables at each
     * measured point. MAX_INDEP is the maximum number of independent variables
     */
    n = parse_data(args, S, &nvars, &m, params, guess, fitparams, &modelfunc, params_name);
    /* n -> number of points
     * S -> values of the independent variables
     * nvars -> number of independent variables
     * m -> number of parameters
     * params -> values of the parameters
     * guess -> guess of the parameters
     * fitparams -> parameters to fit (1) or keep fixed (0)
     * modelfunc -> the model function  
     */
     if (n < 0) {
         return n;
     }

    /* build final matrix of independent variables values. Note that
     * S is nvars by n, while the matrix must by n by nvars*/
    double X[n][nvars];
    for (i = 0; i < nvars; i++) {
        for (j = 0; j < n; j++) {
            X[j][i] = S[i][j];
        }
        free(S[i]);
    }
    /* Copy parameters and guesses in arrays of the right size and free
     * memory */
    double p[m]; /* parameters */
    double p_guess[m]; /* guess */
    int fit[m]; /* array indicating what parameters to fit/keep fixed */
    int *fitp[m]; /* pointer to fit, to be passed to montecarlo */
    *fitp = fit;
    mfit = 0;
    for (i = 0; i < m; i++) {
        p[i] = params[i];
        p_guess[i] = guess[i];
        if ((fit[i] = fitparams[i]) == 1)
            mfit++;
    }
    free(params);
    free(guess);
    free(fitparams);
    double params_mean[m];
    double params_variance[m];
    if (fp != NULL) {
        fprintf(fp, "Model: %s\nParams:\n%s\nPoints:\n%s\n", args->model,
                args->params, args->data);
    }
    n = montecarlo(modelfunc, p, p_guess, error, NREPS, n, nvars, m, mfit,
        fitp, X, params_mean, params_variance, fp);
    /* --- output --- */
    printf("\n Parameter     Mean      Standard Dev       CV(%%)\n");
    printf("----------------------------------------------------\n");
    for (i = 0; i < m; i++) {
        if (params_variance[i] > 0.001) {
            printf("%d | %-6s %10.6f   %10.4f    %10.2f%%\n", i, params_name[i],
            params_mean[i], sqrt(params_variance[i]), 100*sqrt(params_variance[i])/params_mean[i]);
        } else {
            printf("%d | %-6s %10.6e   %10.4e    %10.2f%%\n", i, params_name[i],
            params_mean[i], sqrt(params_variance[i]), 100*sqrt(params_variance[i])/params_mean[i]);
        }
        printf("----------------------------------------------------\n");
    }
    printf("Number of succesful adjustments: %d of %d (%.2f%%)\n", n, NREPS, 100*(float)n/(float)NREPS);
    return 0;
}

/* parses the data points, returns the number of values if this is
 * equal for all the independent variables, or -1 if they do not.
 * Fills S[] with "nvars" pointers to double arrays, being "nvars" the
 * number of independent variables of the model
 */
int parse_data(struct arguments *args, double *S[], int *nvars, int *m,
               double params[], double guess[], int fit[],
               double (**f)(double X[], double p[]),
               char *params_name[MAX_PARAMS])
{
    int i, j, n[MAX_INDEP];
    char data[MAX_DATA_CHARS], data2[MAX_DATA_CHARS];
    char *pattern_data[MAX_INDEP];
    char *pattern_params[MAX_PARAMS];
    regex_t regex;
    regmatch_t pmatch[1];
    strtok(args->model, "\n");
    /* michaelis-menten */
    if (strcmp(args->model, "michaelis") == 0) {
        (*f) = michaelis;
        *m = 2; /* num of parameters */
        *nvars = 1; /* num of independent variables */
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_params[0] = "Vmax[ ]*=[ ]*[.0-9eE-]*";
        pattern_params[1] = "Km[ ]*=[ ]*[.0-9eE-]*";
        params_name[0] = "Vmax";
        params_name[1] = "Km";

    /* Alberty's equation (multisubstrate general equation) */
    } else if (strcmp(args->model, "alberty") == 0) {
        (*f) = alberty;
        *m = 4;
        *nvars = 2;
        pattern_data[0] = "A[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "B[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "KmA";
        params_name[2] = "KmB";
        params_name[3] = "KsA";
    /* Double substitution (ping-pong) */
    } else if (strcmp(args->model, "pingpong") == 0) {
        (*f) = pingpong;
        *m = 3;
        *nvars = 2;
        pattern_data[0] = "A[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "B[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "KmA";
        params_name[2] = "KmB";

    /* mixed inhibition */
    } else if (strcmp(args->model, "mixed") == 0) {
        (*f) = mixed;
        *m = 4;
        *nvars = 2;
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "I[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "Km"; 
        params_name[2] = "KIa";
        params_name[3] = "KIb";
    /* competitive inhibition */
    } else if (strcmp(args->model, "competitive") == 0) {
        (*f) = competitive;
        *m = 3;
        *nvars = 2;
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "I[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "Km"; 
        params_name[2] = "KIa";
    /* uncompetitive inhibition */
    } else if (strcmp(args->model, "uncompetitive") == 0) {
        (*f) = uncompetitive;
        *m = 3;
        *nvars = 2;
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "I[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "Km"; 
        params_name[2] = "KIb";
    /* noncompetitive inhibition */
    } else if (strcmp(args->model, "noncompetitive") == 0) {
        (*f) = noncompetitive;
        *m = 3;
        *nvars = 2;
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "I[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "Km"; 
        params_name[2] = "KIb";
    /* effect of the pH */
    } else if (strcmp(args->model, "ph") == 0) {
        (*f) = ph;
        *m = 6; 
        *nvars = 2;
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "H[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "Km";
        params_name[2] = "Ka1";
        params_name[3] = "Ka2";
        params_name[4] = "Ka3";
        params_name[5] = "Ka4";
    /* michaelis + temperature */
    } else if (strcmp(args->model, "michaelistemp") == 0) {
        (*f) = michaelistemp;
        *m = 4;
        *nvars = 2;
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "T[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "Km";
        params_name[2] = "Ea";
        params_name[3] = "T1";
    } else if (strcmp(args->model, "michaelisinactiv") == 0) {
        (*f) = michaelisinactiv;
        *m = 3;
        *nvars = 2;
        pattern_data[0] = "S[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        pattern_data[1] = "t[ ]*=[ ]*(\\[[ ,.0-9eE-]*\\])";
        params_name[0] = "Vmax";
        params_name[1] = "Km";
        params_name[2] = "kt";
    } else {
        fprintf(stderr, "Error: unrecognized model name: \"%s\"\n", args->model);
        return -1;
    }
    /* Create pattern (pattern_params) to search for every parameter */
    for (i = 0; i < *m; i++)
        asprintf(&pattern_params[i], "%s%s", params_name[i], "[ ]*=[ ]*[.0-9eE-]*");

    /* build matrix of independent variables values */
    for (i = 0; i < *nvars; i++) {
        if (regcomp(&regex, pattern_data[i], REG_EXTENDED)) {
            fprintf(stderr, "Error: Could not compile regex\n");
            exit(1);
        }
        if (regexec(&regex, args->data, 1, pmatch, 0)) {
            fprintf(stderr, "Error: Could not find all the independent variables.\n");
            exit(1);
        }
        strcopy(data, args->data, pmatch[0].rm_so, pmatch[0].rm_eo);
        regfree(&regex);
        if((n[i] = array_length(data)) < 0)
            return -1;

        S[i] = malloc((n[i]+1) * sizeof(double));
        if (parse_array(S[i], n[i], data) < 0)
            return -1;
    }
    i = n[0];
    for (j = 1; j < *nvars; j++) {
        if (n[j] != i) {
            fprintf(stderr, "Error: the number of values is not equal for all the independent variables\n");
            return -1;
        }
    }
    /* obtain parameters, and the guess of the parameters */
    if (args->guess == NULL) {
        args->guess = args->params;
    }
    for (i = 0; i < *m; i++) {
        /* Compile regular expression */
        if (regcomp(&regex, pattern_params[i], REG_EXTENDED)) {
            fprintf(stderr, "Error: Could not compile regex\n");
            return -1;
        }
        /* Extract first parameter, save it into "data" */
        if (regexec(&regex, args->params, 1, pmatch, 0)) {
            fprintf(stderr, "Error: Could not find all the parameters.\nThis model uses the following parameters:  ");
            for (i = 0; i < *m; i++)
                printf("%s%s", params_name[i], i == *m - 1 ? "\n" : ", ");
            return -1;
        }
        strcopy(data, args->params, pmatch[0].rm_so, pmatch[0].rm_eo);
        /* Extract first parameter guess, save it into "data" */
        if (regexec(&regex, args->guess, 1, pmatch, 0)) {
            fprintf(stderr, "Error: Could not find all the parameters guesses.\nThis model uses the following parameters:");
            for (i = 0; i < *m; i++)
                printf("%s%s", params_name[i], i == *m - 1 ? "\n" : ", ");
            return -1;
        }
        strcopy(data2, args->guess, pmatch[0].rm_so, pmatch[0].rm_eo);
        regfree(&regex);
        /* Check for fixed parameters */
        if (args->fixed_params != NULL) {
            regcomp(&regex, params_name[i], REG_EXTENDED);
            if (regexec(&regex, args->fixed_params, 0, 0, 0))
                fit[i] = 1;
            else
                fit[i] = 0;
            regfree(&regex);
        } else {
            fit[i] = 1;
        }
        /* Extract the value of the parameter */
        if (regcomp(&regex, "[0-9.]*$", REG_EXTENDED)) {
            fprintf(stderr, "Error: Could not compile regex\n");
            return -1;
        }
        if (regexec(&regex, data, 1, pmatch, 0)) {
            fprintf(stderr, "Error: Could not find the value of all the parameters.");
            return -1;
        }
        strcopy(data, data, pmatch[0].rm_so, pmatch[0].rm_eo);
        /* Extract the value of the parameter guess */
        if (regexec(&regex, data2, 1, pmatch, 0)) {
            fprintf(stderr, "Error: Could not find the value of all the parameters guesses");
            return -1;
        }
        strcopy(data2, data2, pmatch[0].rm_so, pmatch[0].rm_eo);
        regfree(&regex);
        /* Save value of first parameter and first guess*/
        sscanf(data, "%lf", &params[i]);
        sscanf(data2, "%lf", &guess[i]);
    }
    return n[0];

}

/* checks how many elements has the array represented by the string S */
int array_length(char *S)
{
    char c, *S2 = S;
    double w;
    int i = 0, j = 0;
    /* first number */
    while (!sscanf(S2, "[%e", &w)) {
        if (S2[0] == '\0') {
            fprintf(stderr, "Error in array of substrate concentrations\n");
            return -1;
        }
        S2++;
    }
    i++;
    /* jump comma */
    while (S2[0] != '\0') {
        for (c = *(S+j); c != ',' && c != '\0' && c != ' '; c = *(S+j)) {
            j++;
        }
        if (c == '\0')
            break;
        j++;
        S2 = S+j;
        if (sscanf(S2, "%e,", &w))
            i++;
    }
    /* last number */
    if (sscanf(S, "%e]", &w))
        i++;
    return i;
}

/* parses the "nvars" strings in ... and saves the values in S_arr */
int parse_array(double *S_arr, int n, char *S)
{
    char c, *S2;
    int i = 0, j = 0;
    S2 = S;
    while (!sscanf(S2, "[%lf", &S_arr[i])) {
        if (S2[0] == '\0') {
            fprintf(stderr, "Error in array of substrate concentrations\n");
            return -1;
        }
        S2++;
    }
    i++;
    /* jump comma */
    while (S2[0] != '\0') {
        for (c = *(S+j); c != ',' && c != '\0' && c != ' '; c = *(S+j)) {
            j++;
        }
        if (c == '\0')
            break;
        j++;
        S2 = S+j;
        if (sscanf(S2, "%lf,", &S_arr[i]))
            i++;
    }
    /* last number */
    if (sscanf(S, "%lf]", &S_arr[i]))
        i++;
    return i;
}

int strcopy(char *dst, char *src, int off, int end)
{
    int i;
    for (i = 0; i + off < end; i++) {
        dst[i] = src[i+off];
    }
    dst[i] = '\0';
    return i;
}
