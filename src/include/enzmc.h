#ifndef _ENZMC_

#define NORMAL_MODE 52
#define INTERACTIVE_MODE 53
#define GRAPHICAL_MODE 54
#define FILE_MODE 55
#define TEMPLATE_MODE 56

/* error messages */
#define ERROR_LACK_OPTS "lacking options (model, params, data and error are mandatory"
#define ERROR_TOO_MANY_ARGS "too many arguments (see --usage)"
#define ERROR_NO_FILENAME "you must specify a file name:\n./enzmc --template model filename"

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

/* Program modes */
int run_interactive_mode(struct arguments *);
int run_graphical_mode(struct arguments *);
int run_normal_mode(struct arguments *);

/* Internal stuff */
int parse_array(double *S_arr, int n, char S_str[]);
int strcopy(char *dst, char *src, int start, int end);


#endif
