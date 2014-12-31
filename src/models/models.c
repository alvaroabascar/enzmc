#include <stdlib.h>
#ifndef _ENZYME_
#include "enzyme.h"
#endif

#ifndef _MODELS_
#include "models.h"
#endif

struct model {
  char *name;                                  // name of the model
  double (*function) (double X[], double p[]); // function (ex. see enzyme.h)
  int nparams;                                 // number of parameters
  int nvars;                                   // number of indep vars
  char *params[MAX_PARAMS];                    // names of the parameters
  char *indep_vars[MAX_INDEP];                 // names of the indep vars
};

struct model models[10] = {
  {
    "michaelis",
    michaelis,
    2,
    1,
    {"Vmax", "Km"},
    {"S"}
  },
  {
    "alberty",
    alberty,
    3,
    2,
    {"Vmax", "KmA", "KmB"},
    {"A", "B"}
  },
  {
    "pingpong",
    pingpong,
    3,
    2,
    {"Vmax", "KmA", "KmB"},
    {"A", "B"}
  },
  {
    "mixed",
    mixed,
    4,
    2,
    {"Vmax", "Km", "KIa", "KIb"},
    {"S", "I"}
  },
  {
    "competitive",
    competitive,
    3,
    2,
    {"Vmax", "Km", "KIa"},
    {"S", "I"}
  },
  {
    "uncompetitive",
    uncompetitive,
    3,
    2,
    {"Vmax", "Km", "KIb"},
    {"S", "I"}
  },
  {
    "noncompetitive",
    noncompetitive,
    3,
    2,
    {"Vmax", "Km", "KIb"},
    {"S", "I"}
  },
  {
    "ph",
    ph,
    6,
    2,
    {"Vmax", "Km", "Ka1", "Ka2", "Ka3", "Ka4"},
    {"S", "H"}
  },
  {
    "michaelistemp",
    michaelistemp,
    3,
    2,
    {"Ea", "Km", "T1"},
    {"S", "T"}
  }
};
