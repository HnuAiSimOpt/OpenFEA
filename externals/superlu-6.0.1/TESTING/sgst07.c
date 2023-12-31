/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

/*! \file
 * SGST07 tests the error bounds from iterative refinement for the
 * computed solution to a system of equations op(A)*X = B, where A is a
 * general n by n matrix and op(A) = A or A**T, depending on TRANS.
 *
 * \ingroup TestingS
 */

#include <math.h>
#include "slu_sdefs.h"

/*!
 * SGST07 tests the error bounds from iterative refinement for the
 * computed solution to a system of equations op(A)*X = B, where A is a
 * general n by n matrix and op(A) = A or A**T, depending on TRANS.
 *
 * RESLTS(1) = test of the error bound
 *           = norm(X - XACT) / ( norm(X) * FERR )
 * A large value is returned if this ratio is not less than one.
 *
 * RESLTS(2) = residual from the iterative refinement routine
 *           = the maximum of BERR / ( (n+1)*EPS + (*) ), where
 *             (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )
 *
 * \param[in] trans   Specifies the form of the system of equations.
 *                    = NOTRANS:  A *x = b
 *                    = TRANS  :  A'*x = b, where A' is the transpose of A
 *                    = CONJ   :  A'*x = b, where A' is the transpose of A
 * \param[in] n       The number of rows of the matrices X and XACT.  N >= 0.
 * \param[in] nrhs    The number of columns of the matrices X and XACT.  NRHS >= 0.
 * \param[in] A       The original n by n matrix A, dimension(A->nrow, A->ncol).
 * \param[in] b       The right hand side vectors for the system of linear
 *                    equations, dimension (LDB,NRHS).
 * \param[in] ldb     The leading dimension of the array B.  LDB >= max(1,N).
 * \param[in] x       The computed solution vectors, dimension(LDX,NRHS).
 *                    Each vector is stored as a column of the matrix X.
 * \param[in] ldx     The leading dimension of the array X.  LDX >= max(1,N).
 * \param[in] xact    The exact solution vectors, dimension(LDX,NRHS).
 *                    Each vector is stored as a column of the matrix XACT.
 * \param[in] ldxact  The leading dimension of the array XACT.  LDXACT >= max(1,N).
 * \param[in] ferr    The estimated forward error bounds for each solution
 *                    vector X, dimension(NRHS).
 *                    If XTRUE is the true solution, FERR bounds the magnitude of
 *                    the largest entry in (X - XTRUE) divided by the magnitude
 *                    of the largest entry in X.
 * \param[in] berr    The componentwise relative backward error of each solution
 *                    vector (i.e., the smallest relative change in any entry of A
 *                    or B that makes X an exact solution), dimension(NRHS).
 * \param[out] reslts The maximum over the NRHS solution vectors of the ratios, dimension (2):
 *                    RESLTS(1) = norm(X - XACT) / ( norm(X) * FERR )
 *                    RESLTS(2) = BERR / ( (n+1)*EPS + (*) )
 */
int sgst07(trans_t trans, int n, int nrhs, SuperMatrix *A, float *b,
           int ldb, float *x, int ldx, float *xact,
           int ldxact, float *ferr, float *berr, float *reslts)
{
    /* Table of constant values */
    int c__1 = 1;

    /* System generated locals */
    float d__1, d__2;

    /* Local variables */
    float diff, axbi;
    int    imax, irow, n__1;
    int    i, j, k;
    float unfl, ovfl;
    float xnorm;
    float errbnd;
    int    notran;
    float eps, tmp;
    float *rwork;
    float *Aval;
    NCformat *Astore;

    /* Function prototypes */
    extern int    isamax_(int *, float *, int *);


    /* Quick exit if N = 0 or NRHS = 0. */
    if ( n <= 0 || nrhs <= 0 ) {
	reslts[0] = 0.;
	reslts[1] = 0.;
	return 0;
    }

    eps = smach("Epsilon");
    unfl = smach("Safe minimum");
    ovfl   = 1. / unfl;
    notran = (trans == NOTRANS);

    rwork  = (float *) SUPERLU_MALLOC(n*sizeof(float));
    if ( !rwork ) ABORT("SUPERLU_MALLOC fails for rwork");
    Astore = A->Store;
    Aval   = (float *) Astore->nzval;
    
    /* Test 1:  Compute the maximum of   
       norm(X - XACT) / ( norm(X) * FERR )   
       over all the vectors X and XACT using the infinity-norm. */

    errbnd = 0.;
    for (j = 0; j < nrhs; ++j) {
	n__1 = n;
	imax = isamax_(&n__1, &x[j*ldx], &c__1);
	d__1 = fabs(x[imax-1 + j*ldx]);
	xnorm = SUPERLU_MAX(d__1,unfl);
	diff = 0.;
	for (i = 0; i < n; ++i) {
	    d__1 = fabs(x[i+j*ldx] - xact[i+j*ldxact]);
	    diff = SUPERLU_MAX(diff, d__1);
	}

	if (xnorm > 1.) {
	    goto L20;
	} else if (diff <= ovfl * xnorm) {
	    goto L20;
	} else {
	    errbnd = 1. / eps;
	    goto L30;
	}

L20:
#if 0	
	if (diff / xnorm <= ferr[j]) {
	    d__1 = diff / xnorm / ferr[j];
	    errbnd = SUPERLU_MAX(errbnd,d__1);
	} else {
	    errbnd = 1. / eps;
	}
#endif
	d__1 = diff / xnorm / ferr[j];
	errbnd = SUPERLU_MAX(errbnd,d__1);
	/*printf("Ferr: %f\n", errbnd);*/
L30:
	;
    }
    reslts[0] = errbnd;

    /* Test 2: Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where 
       (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) + abs(b))_i ) */

    for (k = 0; k < nrhs; ++k) {
	for (i = 0; i < n; ++i) 
            rwork[i] = fabs( b[i + k*ldb] );
	if ( notran ) {
	    for (j = 0; j < n; ++j) {
		tmp = fabs( x[j + k*ldx] );
		for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
		    rwork[Astore->rowind[i]] += fabs(Aval[i]) * tmp;
                }
	    }
	} else {
	    for (j = 0; j < n; ++j) {
		tmp = 0.;
		for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
		    irow = Astore->rowind[i];
		    d__1 = fabs( x[irow + k*ldx] );
		    tmp += fabs(Aval[i]) * d__1;
		}
		rwork[j] += tmp;
	    }
	}

	axbi = rwork[0];
	for (i = 1; i < n; ++i) axbi = SUPERLU_MIN(axbi, rwork[i]);
	
	/* Computing MAX */
	d__1 = axbi, d__2 = (n + 1) * unfl;
	tmp = berr[k] / ((n + 1) * eps + (n + 1) * unfl / SUPERLU_MAX(d__1,d__2));
	
	if (k == 0) {
	    reslts[1] = tmp;
	} else {
	    reslts[1] = SUPERLU_MAX(reslts[1],tmp);
	}
    }

    SUPERLU_FREE(rwork);
    return 0;

} /* sgst07 */
