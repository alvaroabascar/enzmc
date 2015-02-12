#ifndef __ENZMC_H__
#define __ENZMC_H__

#include <models.h>

/* Stuff directly related to the implementation of the cli */

#define MAX_VALUES 200

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
int get_indep_vars(struct model *model, char *raw_data,
                   double *data[model->nvars]);
int get_params(struct model *model, char *raw_data, double *params);
int get_error(struct model *model, char *raw_data, double *error);
int get_fixed_params(struct model *model, char *raw_data, int *fixed_ptr);

int extract_str(char *src, char *dst, char *regexp_str);
int parse_array_double(char *array_str, double **dst);

void free_matrix_double(double *array[], int n);
void copy_array_double(double *dst, double *src, int n);
#endif /* __ENZMC_H__ */
