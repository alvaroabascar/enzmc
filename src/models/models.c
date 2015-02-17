#include "models.h"
#include "enzyme.h"
#include "enzyme.c"

/* How to implement a model?
 *
 * 1. In enzyme.c, define your function
 * 2. In enzyme.h, add the function prototype
 * 3. Here, add a new entry defining the name, name of the function
 * (as in enzyme.c), num of parameters and indep variables, names of the
 * parameters and name of the independent variables.
 *
 * NOTE: enzyme.c and enzyme.h contain enzymatic models. If you want to include
 * models from a field not related to enzymology, you might want to create a
 * couple of files: e.g. phisics.h and phisics.c and place your functions there.
 * Then just include them at the top of this file as follows:
 *
 * #include "phisics.h"
 * #include "phisics.c"
 *
 * IMPORTANT: leave the last "model" as is: it is used to recognized the end of
 * the array of models. Modifying it will break the program
 */

struct model models[] = {
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
  },
  {
    "",
    NULL,
    0,
    0,
    {},
    {}
  }
};
