#ifndef __ENZMC_H__
#define __ENZMC_H__

/* Stuff directly related to the implementation of the cli */

/* program modes */
#define NORMAL_MODE 52
#define FILE_MODE 55
#define TEMPLATE_MODE 56

/* Error messages */
#define ERROR_LACK_OPTS "lacking options (model, params, data and error are mandatory"
#define ERROR_TOO_MANY_ARGS "too many arguments (see --usage)"
#define ERROR_NO_FILENAME "you must specify a file name:\n./enzmc --template model filename"

/* Arguments which must be provided to do the simulation. */
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
int run_cli_mode(struct arguments *);
int run_file_mode(struct arguments *);
int create_template(char *modelname, char *fileout);

/* Internal functions */
struct model *get_model(char *modelname);

#endif /* __ENZMC_H__ */
