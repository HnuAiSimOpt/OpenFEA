/*******************************************************************************
* Copyright 2017 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
!  Content:
!      Intel(R) Math Kernel Library (Intel(R) MKL) interfaces for Compact format
!******************************************************************************/
#ifndef _MKL_COMPACT_H
#define _MKL_COMPACT_H

#include "mkl_types.h"

#ifndef mkl_compact_complex_float
#define mkl_compact_complex_float MKL_Complex8
#endif

#ifndef mkl_compact_complex_double
#define mkl_compact_complex_double MKL_Complex16
#endif

#ifdef __cplusplus
extern "C" {            /* Assume C declarations for C++ */
#endif /* __cplusplus */

MKL_COMPACT_PACK mkl_get_format_compact( void );

MKL_INT mkl_sget_size_compact( const MKL_INT LD, const MKL_INT SD,
                               const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_sgepack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          const float * const *A, const MKL_INT LDA, float *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_sgeunpack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          float * const *A, const MKL_INT LDA, const float *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_sgemm_compact( const  MKL_LAYOUT Layout, const  MKL_TRANSPOSE TransA,
                         const  MKL_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                         const MKL_INT K, const float alpha, const float *A,
                         const MKL_INT a_stride, const float *B, const MKL_INT b_stride,
                         const float beta, float *C, const MKL_INT c_stride,
                         const MKL_COMPACT_PACK format, const MKL_INT NM );
void mkl_strsm_compact( const  MKL_LAYOUT Layout, const  MKL_SIDE Side,
                         const  MKL_UPLO Uplo, const  MKL_TRANSPOSE TransA,
                         const  MKL_DIAG Diag, const MKL_INT M, const MKL_INT N,
                         const float alpha, const float *A, const MKL_INT a_stride,
                         float *B, const MKL_INT b_stride,
                         const MKL_COMPACT_PACK format, const MKL_INT NM );

MKL_INT mkl_dget_size_compact( const MKL_INT LD, const MKL_INT SD,
                               const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_dgepack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          const double * const *A, const MKL_INT LDA, double *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_dgeunpack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          double * const *A, const MKL_INT LDA, const double *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_dgemm_compact( const  MKL_LAYOUT Layout, const  MKL_TRANSPOSE TransA,
                     const  MKL_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                     const MKL_INT K, const double alpha, const double *A,
                     const MKL_INT a_stride, const double *B, const MKL_INT b_stride,
                     const double beta, double *C, const MKL_INT c_stride,
                     const MKL_COMPACT_PACK format, const MKL_INT NM );
void mkl_dtrsm_compact( const  MKL_LAYOUT Layout, const  MKL_SIDE Side,
                         const  MKL_UPLO Uplo, const  MKL_TRANSPOSE TransA,
                         const  MKL_DIAG Diag, const MKL_INT M, const MKL_INT N,
                         const double alpha, const double *A, const MKL_INT a_stride,
                         double *B, const MKL_INT b_stride,
                         const MKL_COMPACT_PACK format, const MKL_INT NM );

MKL_INT mkl_cget_size_compact( const MKL_INT LD, const MKL_INT SD,
                               const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_cgepack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          const mkl_compact_complex_float * const *A, const MKL_INT LDA, float *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_cgeunpack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          mkl_compact_complex_float * const *A, const MKL_INT LDA, const float *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_cgemm_compact( const MKL_LAYOUT Layout, const  MKL_TRANSPOSE TransA,
                         const MKL_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                         const MKL_INT K, const mkl_compact_complex_float *alpha, const float *A,
                         const MKL_INT a_stride, const float *B, const MKL_INT b_stride,
                         const mkl_compact_complex_float *beta, float *C, const MKL_INT c_stride,
                         const MKL_COMPACT_PACK format, const MKL_INT NM );
void mkl_ctrsm_compact( const  MKL_LAYOUT Layout, const  MKL_SIDE Side,
                         const  MKL_UPLO Uplo, const  MKL_TRANSPOSE TransA,
                         const  MKL_DIAG Diag, const MKL_INT M, const MKL_INT N,
                         const mkl_compact_complex_float *alpha, const float *A, const MKL_INT a_stride,
                         float *B, const MKL_INT b_stride,
                         const MKL_COMPACT_PACK format, const MKL_INT NM );

MKL_INT mkl_zget_size_compact( const MKL_INT LD, const MKL_INT SD,
                               const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_zgepack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          const mkl_compact_complex_double * const *A, const MKL_INT LDA, double *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_zgeunpack_compact( const MKL_LAYOUT Layout, const MKL_INT ROWS, const MKL_INT COLUMNS,
                          mkl_compact_complex_double * const *A, const MKL_INT LDA, const double *Ap, const MKL_INT LDAP,
                          const MKL_COMPACT_PACK FORMAT, const MKL_INT NM );
void mkl_zgemm_compact( const  MKL_LAYOUT Layout, const  MKL_TRANSPOSE TransA,
                         const  MKL_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
                         const MKL_INT K, const mkl_compact_complex_double *alpha, const double *A,
                         const MKL_INT a_stride, const double *B, const MKL_INT b_stride,
                         const mkl_compact_complex_double *beta, double *C, const MKL_INT c_stride,
                         const MKL_COMPACT_PACK format, const MKL_INT NM );
void mkl_ztrsm_compact( const  MKL_LAYOUT Layout, const  MKL_SIDE Side,
                         const  MKL_UPLO Uplo, const  MKL_TRANSPOSE TransA,
                         const  MKL_DIAG Diag, const MKL_INT M, const MKL_INT N,
                         const mkl_compact_complex_double *alpha, const double *A, const MKL_INT a_stride,
                         double *B, const MKL_INT b_stride,
                         const MKL_COMPACT_PACK format, const MKL_INT NM );

/* LAPACK compact routines */

void mkl_cgetrinp_compact( MKL_LAYOUT layout, MKL_INT n, float* ap,
                           MKL_INT ldap, float* work, MKL_INT lwork,
                           MKL_INT* info, MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_dgetrinp_compact( MKL_LAYOUT layout, MKL_INT n, double* ap, MKL_INT ldap,
                           double* work, MKL_INT lwork, MKL_INT* info,
                           MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_sgetrinp_compact( MKL_LAYOUT layout, MKL_INT n, float* ap, MKL_INT ldap,
                           float* work, MKL_INT lwork, MKL_INT* info,
                           MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_zgetrinp_compact( MKL_LAYOUT layout, MKL_INT n, double* ap,
                           MKL_INT ldap, double* work, MKL_INT lwork,
                           MKL_INT* info, MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_cgetrfnp_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n,
                           float* ap, MKL_INT ldap, MKL_INT* info,
                           MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_dgetrfnp_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n, double* ap,
                           MKL_INT ldap, MKL_INT* info, MKL_COMPACT_PACK format,
                           MKL_INT nmat );

void mkl_sgetrfnp_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n, float* ap,
                           MKL_INT ldap, MKL_INT* info, MKL_COMPACT_PACK format,
                           MKL_INT nmat );

void mkl_zgetrfnp_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n,
                           double* ap, MKL_INT ldap, MKL_INT* info,
                           MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_cpotrf_compact( MKL_LAYOUT layout, MKL_UPLO uplo, MKL_INT n,
                         float* ap, MKL_INT ldap, MKL_INT* info,
                         MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_dpotrf_compact( MKL_LAYOUT layout, MKL_UPLO uplo, MKL_INT n, double* ap,
                         MKL_INT ldap, MKL_INT* info, MKL_COMPACT_PACK format,
                         MKL_INT nmat );

void mkl_spotrf_compact( MKL_LAYOUT layout, MKL_UPLO uplo, MKL_INT n, float* ap,
                         MKL_INT ldap, MKL_INT* info, MKL_COMPACT_PACK format,
                         MKL_INT nmat );

void mkl_zpotrf_compact( MKL_LAYOUT layout, MKL_UPLO uplo, MKL_INT n,
                         double* ap, MKL_INT ldap, MKL_INT* info,
                         MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_cgeqrf_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n, float* ap,
                         MKL_INT ldap, float* taup, float* work, MKL_INT lwork,
                         MKL_INT* info, MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_dgeqrf_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n, double* ap,
                         MKL_INT ldap, double* taup, double* work,
                         MKL_INT lwork, MKL_INT* info, MKL_COMPACT_PACK format,
                         MKL_INT nmat );

void mkl_sgeqrf_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n, float* ap,
                         MKL_INT ldap, float* taup, float* work, MKL_INT lwork,
                         MKL_INT* info, MKL_COMPACT_PACK format, MKL_INT nmat );

void mkl_zgeqrf_compact( MKL_LAYOUT layout, MKL_INT m, MKL_INT n, double* ap,
                         MKL_INT ldap, double* taup, double* work,
                         MKL_INT lwork, MKL_INT* info, MKL_COMPACT_PACK format,
                         MKL_INT nmat );

#ifdef __cplusplus
}
#endif    /* __cplusplus */

#endif /* _MKL_COMPACT_H */
