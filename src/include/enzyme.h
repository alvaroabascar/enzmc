#ifndef __MICHAELIS__

double michaelis(double S[1], double p[2]);
double michaelis_dyda(double S[1], double p[2], int k);
double alberty(double S[2], double p[4]);
double pingpong(double S[2], double p[3]);
double mixed(double S[2], double p[4]);
double competitive(double S[2], double p[3]);
double uncompetitive(double S[2], double p[3]);
double noncompetitive(double S[2], double p[3]);
double ph(double S[2], double p[6]);
double michaelistemp(double S[2], double p[4]);
double michaelisinactiv(double S[2], double p[3]);

#endif
