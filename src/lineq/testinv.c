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
   static double A[N][N] = {{2,3,4},
                            {3,2,1},
                            {2,1,2}};
   static double B[N][M] = {{1,0,0},
                            {0,1,0},
                            {0,0,1}};
   double id;
   /*srand(time(NULL));
   for (i = 0; i < n; i++) {
       B[i][0] = 0;
       for (j = 0; j < n; j++) {
           A[i][j] = rand() % 100;
           for (k = 0; k < M; k++)
               B[i][k] += A[i][j] * (j+1);
       }
   }*/
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
