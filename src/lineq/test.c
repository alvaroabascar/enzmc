#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 3
#define M 3
int main()
{
   int n = N;
   int m = M;
   int i, j, k;
   static double A[N][N];
   static double B[N][M];
   double id;
   srand(time(NULL));
   for (i = 0; i < n; i++) {
       B[i][0] = 0;
       for (j = 0; j < n; j++) {
           A[i][j] = rand() % 100;
           for (k = 0; k < M; k++)
               B[i][k] += A[i][j] * (j+1);
       }
   }
   printf("Modify N and M in test.c to obtain systems of different size\n");
   printf("Every column of B must be a vector [1, 2, 3, ... N]\n\n");
   printf("Before gaussj:\n");
   printf("A:\n");
   mprint(n, n, A);
   printf("B:\n");
   mprint(n, m, B);
   if (gaussj(N, M, A, B) != 0) {
       printf("singular matrix\n");
       return -1;
   }
   printf("After gaussj:\n");
   printf("A:\n");
   mprint(n, n, A);
   printf("B:\n");
   mprint(n, m, B);
   return 0;
}
