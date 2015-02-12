#include "gaussjbs.h"
#include "../misc/matrix.h"

#ifndef abs
#define abs(x) ((x>=0) ? x : -1*x)
#endif

/* gaussj:
 * given an n-by-n matrix of coefficients A, and a matrix of m sets
 * of right values B (n-by-m), performs Gauss-Jordan elimination
 * with backsubstitution of A and returns the solutions (X) in B.
 *
 * Note: this version uses full, explicit pivoting.
 */
int gaussj(int n, int m, double A[][n], double B[][m])
{
    int i, j, k, row, max_row, max_col, itmp;
    int piv[n]; /* keeps track of column pivoting */
    double tmp; /* used to switch rows/cols */
    double max;
    /* set colpiv to initial state: each column (index) stores the original
    column (value) */
    for (i = 0; i < n; i++) {
        piv[i] = i;
    }

    /* outer loop: iterate over each diagonal element */
    for (row = 0; row < n; row++) {
        /* find best pivot in row (at right of A[row][row]): */
        max = 0;
        for (i = row; i < n; i++) {
            if (abs(A[row][i]) > abs(max)) {
                max = A[row][i];
                max_row = row;
                max_col = i;
            }
        }
        /* find best pivot in column (below A[row][row]): */
        for (i = row; i < n; i++) {
            if (abs(A[i][row]) > abs(max)) {
                max = A[i][row];
                max_row = i;
                max_col = row;
            }
        }
        /* find best pivot in row (at right of A[row][row]): */
        /* if best pivot == 0 -> singular matrix */
        if (abs(max) == 0) {
            return -1;
        }
        /* pivot in same column, switch rows */
        if (max_col == row && max_row != row) {
            /* switch rows "row" and "max_row" in A and B*/
            switchrows(n, n, A, row, max_row);
            switchrows(n, m, B, row, max_row);
        /* pivot in same row, switch columns */
        } else if (max_col != row && max_row == row) {
            switchcols(n, n, A, row, max_col);
            /* keep track of column switchs */
            i = piv[row];
            piv[row] = piv[max_col];
            piv[max_col] = i;
        }
        /* pivoting has been completed, start operating:
        max = 0;
        1. divide row by pivot */
        for (i = 0; i < n; i++) {
            A[row][i] /= max;
        }
        /* same for B: */
        for (i = 0; i < m; i++) {
            B[row][i] /= max;
        }
        /* 2. for each element below in the column, substract to
        its row this row multiplied by that element */
        for (i = row + 1; i < n; i++) {
            tmp = A[i][row]; /* multiplier */
            for (j = row; j < n; j++) {
                A[i][j] = A[i][j] - tmp * A[row][j];
            }
            for (j = 0; j < m; j++) {
                B[i][j] -= tmp * B[row][j];
            }
        }
    }
    /* back-substitution: */
    for (k = 0; k < m; k++) {
        for (i = n - 1; i >= 0; i--) {
            tmp = 0;
            for (j = n - 1; j > i; j--) {
                tmp += B[j][k]*A[i][j];
            }
            B[i][k] -= tmp;
        }
    }
    /* reverse column changes to obtain B in the correct order */
    for (i = 0; i < n; i++) {
        while (piv[i] != i) {
            switchrows(n, m, B, i, piv[i]);
            /* same switch in piv */
            itmp = piv[i];
            piv[i] = piv[itmp];
            piv[itmp] = itmp;
        }
    }
    return 0;
}
