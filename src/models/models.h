#ifndef __MODELS_H__
#define __MODELS_H__

/* maximum number of independent variables in a model*/
#define MAX_INDEP 10
/* max num of chars in the --data option */
#define MAX_DATA_CHARS 1000
/* max number of parameters of a model */
#define MAX_PARAMS 10

struct model {
  char *name;                                  // name of the model
  double (*function) (double X[], double p[]); // function (ex. see enzyme.h)
  int nparams;                                 // number of parameters
  int nvars;                                   // number of indep vars
  char *params[MAX_PARAMS];                    // names of the parameters
  char *indep_vars[MAX_INDEP];                 // names of the indep vars
};

#endif /* __MODELS_H__ */
