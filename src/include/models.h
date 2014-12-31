#ifndef _ENZYME_
#include "enzyme.h"
#endif

#ifndef _MODELS_

/* maximum number of independent variables */
#define MAX_INDEP 10
/* maximum number of chars in the --data field */
#define MAX_DATA_CHARS 1000
/* maximum number of parameters */
#define MAX_PARAMS 10

struct model {
  char *name;                                  // name of the model
  double (*function) (double X[], double p[]); // function (ex. see enzyme.h)
  int nparams;                                 // number of parameters
  int nvars;                                   // number of indep vars
  char *params[MAX_PARAMS];                    // names of the parameters
  char *indep_vars[MAX_INDEP];                 // names of the indep vars
};

#endif
