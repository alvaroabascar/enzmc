
/*
 **********************[ IMPORTANT NOTE ]***********************
 * This version of GE does NOT build the inverse matrix. Until *
 * this is fixed, it makes no sense to use this version.       *
 * The version contained in gaussjbs (gauss-jordan elimination *
 * with backsubstitution) should be used instead.              *
 ***************************************************************
 */

int gaussj(int n, int m, double [n][n], double [m][n]);

/* This code is my own implementation of Gauss-Jordan, according
 * to the instructions (but not to de code) from numerical recipes
 * in Pascal */

/* gausj:
 * given an n-by-n matrix of coefficients A, and a matrix of m sets
 * of right values B (n-by-m), performs Gauss-Jordan elimination
 * (without backsubstitution) of A and returns the solutions (X) in B.
 * A does NOT contain its inverse at the end of the process.
 *
 * RETURNS: 0 on success, -1 on failure (singular matrix)
 *
 * Note: this version uses full, explicit pivoting.
 */
int gaussj(int n, int m, double A[][n], double B[][m])
{
    int i, j, k, row, col, max_row, max_col, itmp;
    int piv[n]; /* keeps track of column pivoting */
    double tmp; /* used to interchange rows/cols */
    double max;

    /* set colpiv to initial state: each column (index) stores the original
    column (value) */
    for (i = 0; i < n; i++)
        piv[i] = i;

    /* outer loop: iterate over each diagonal element */
    for (row = 0; row < n; row++) {
        /* find best pivot in row (right of A[row][row]): */
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
                max_row = row;
                max_row = i;
                max_col = row;
            }
        }
        /* if best pivot == 0 -> singular matrix */
        if (abs(max) == 0) {
            return -1;
        }
        /* pivot in same column, interchange rows */
        if (max_col == row && max_row != row) {
            for (i = 0; i < n; i++) { /* interchange in A */
                tmp = A[row][i];
                A[row][i] = A[max_row][i];
                A[max_row][i] = tmp;
            }
            for (i = 0; i < m; i++) { /* interchange in B */
                tmp = B[row][i];
                B[row][i] = B[max_row][i];
                B[max_row][i] = tmp;
            }
        /* pivot in same row, interchange columns */
        } else if (max_row == row && max_col != row) {
            for (i = 0; i < n; i++) {
                tmp = A[i][row];
                A[i][row] = A[i][max_col];
                A[i][max_col] = tmp;
            }
        }
        /* keep track of column interchanges */
        i = piv[row];
        piv[row] = piv[max_col];
        piv[max_col] = i;

        /* pivoting has been completed, start operating:
        1. divide row by pivot */
        for (i = 0; i < n; i++) {
            A[row][i] /= max;
        }
        /* same for B: */
        for (i = 0; i < m; i++) {
            B[row][i] /= max;
        }
        /* 2. for each element in the column, substract to
        its row this row multiplied by that element */
        for (i = 0; i < n; i++) {
            if (i == row)
                continue; /* skip current row */
            tmp = A[i][row]; /* multiplier */
            for (j = row; j < n; j++) {
                A[i][j] = A[i][j] - tmp * A[row][j];
            }
            for (j = 0; j < m; j++) {
                B[i][j] -= tmp * B[row][j];
            }
        }
    }

    /* reverse column changes to obtain B in the correct order */
    for (i = 0; i < n; i++) {
        while (piv[i] != i) {
            for (j = 0; j < m; j++)  {
                tmp = B[i][j];
                B[i][j] = B[piv[i]][j];
                B[piv[i]][j] = tmp;
            }
            itmp = piv[i];
            piv[i] = piv[itmp];
            piv[itmp] = itmp;
        }
    }
    return 0;
}
