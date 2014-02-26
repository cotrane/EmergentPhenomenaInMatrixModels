#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>

/*
 * implementations of EV routines from LAPACK; they are used in EVonthefly, EVtoDist_Modulus_nH, EVtoDist_nonHermitian,   to extract eigenvalues
 * for distribution
 */

// computes eigenvalues of hermitian matrix
int findEigenvalues(doublereal *pVec, doublecomplex *m_pMatrix, integer m_nDim)
{
	//assert(m_nDim == vec.m_nDim);
	char jobz = 'N';	// specifies the desired output- 'N' only eigenvalues, 'V' eigenvalues and eigenvectors
	char stor = 'U';	// values 'U'/'L':  Upper/Lower triangle of the matrix is stored;

	integer lda;			// The leading dimension of the array A. LDA > = max(1,N) !?!?
	//__CLPK_real eig_vals[MATRIX_DIM];	// Array which contains the eigenvalues after the cheev_() cal
	doublecomplex wkopt;			// On exit, if INFO = 0, WORK(1) returns the optimal LDWORK. ?!?!?
	doublecomplex* work;			// we should allocate memory here depending on wkopt from the first call
	integer lwork = -1;		// the size of the workspace, if -1, then a workspace size query is assumed
	doublereal *rwork;			// another workspace ?!?!?
	integer info;			// output parameter: 0 on OK

	//static int bInitialized = 0;
	lda = m_nDim;

	// some initialization
	//if(!bInitialized)
	//{
		rwork = (doublereal*) malloc( (3 * m_nDim - 2) * sizeof(doublereal) );// don't know why, this is the adviced value
		// first call of cheev_- no calcs performed. Just fill lwork so we allocate enough memory in work pointer
		zheev_(&jobz, &stor, &m_nDim, m_pMatrix, &lda, pVec, &wkopt, &lwork, rwork, &info);
		// memory allocation
		lwork = (int) wkopt.r;
		work = (doublecomplex*) malloc( lwork * sizeof(doublecomplex) );
	//	bInitialized = 1;
	//}

	//printf("%d line work = 0x%x size = %d\n", __LINE__, (unsigned int) work, lwork * sizeof(__CLPK_complex));

	// 2nd call to cheev_- actual computing
	zheev_(&jobz, &stor, &m_nDim, m_pMatrix, &lda, pVec, work, &lwork, rwork, &info);

	// release allocated memory
	free( work );
	free( (void*) rwork);

	return info;
}

// computes eigenvalues AND eigenvectors of hermitian matrix
int findEigenvalues_vectors(doublereal *pVec, doublecomplex *m_pMatrix, integer m_nDim)
{
//	printmatEV(m_pMatrix, (int) m_nDim);

	//assert(m_nDim == vec.m_nDim);
	char jobz = 'V';	// specifies the desired output- 'N' only eigenvalues, 'V' eigenvalues and eigenvectors
	char uplo = 'L';	// values 'U'/'L':  Upper/Lower triangle of the matrix is stored;

	integer lda;			// The leading dimension of the array A. LDA > = max(1,N) !?!?
	//__CLPK_real eig_vals[MATRIX_DIM];	// Array which contains the eigenvalues after the cheev_() cal
	doublecomplex wkopt;			// On exit, if INFO = 0, WORK(1) returns the optimal LDWORK. ?!?!?
	doublecomplex *work;			// we should allocate memory here depending on wkopt from the first call
	integer lwork = -1;		// the size of the workspace, if -1, then a workspace size query is assumed
	doublereal *rwork;			// another workspace ?!?!?
	integer info;			// output parameter: 0 on OK

	lda = m_nDim;

	// some initialization
	rwork = (doublereal*) malloc( (3 * m_nDim - 2) * sizeof(doublereal) ); // don't know why, this is the adviced value
	// first call of cheev_- no calcs performed. Just fill lwork so we allocate enough memory in work pointer
	zheev_(&jobz, &uplo, &m_nDim, m_pMatrix, &lda, pVec, &wkopt, &lwork, rwork, &info);

	if(info!=0)
		printf("error in memory allocation info=%d\n", (int) info);

	// memory allocation
	lwork = (integer) wkopt.r;
	work = (doublecomplex*) malloc( lwork * sizeof(doublecomplex) );

	//printf("%d line work = 0x%x size = %d\n", __LINE__, (unsigned int) work, lwork * sizeof(__CLPK_complex));

	// 2nd call to cheev_- actual computing
	zheev_(&jobz, &uplo, &m_nDim, m_pMatrix, &lda, pVec, work, &lwork, rwork, &info);

//	printmatEV(m_pMatrix, (int) m_nDim);

	if(info!=0)
		printf("error in computation of eigenvalues/vectors! info=%d\n", (int) info);

	// release allocated memory
	free( work );
	free( (void*) rwork);

	return info;
}

// computes EV's of nonHermitian matrix
int findEigenvalues_NonHermitian(doublecomplex *pVec, doublecomplex *m_pMatrix, integer m_nDim)
{
	//assert(m_nDim == vec.m_nDim);
	char jobvs = 'N';		// specifies the desired output- 'N' only eigenvalues, 'V' eigenvalues and eigenvectors
	char sort = 'N';		// Specifies whether or not to order the eigenvalues on the diagonal of the Schur form.  = 'N': Eigenvalues are not ordered: = 'S': Eigenvalues are ordered (see SELECT).

	integer lda;			// The leading dimension of the array A. LDA > = max(1,N) !?!?
	integer sdim=m_nDim;	// is 0 if sort='N'; If SORT = 'S', SDIM =number of eigenvalues for which SELECT is true.
	doublecomplex *vs;		//(output) COMPLEX*16 array, dimension (LDVS,N); If JOBVS = 'V', VS contains the unitary
							//matrix Z of Schur vectors.  If JOBVS = 'N', VS is not referenced.
	integer ldvs;			//(input) INTEGER; The leading dimension of the array VS.  LDVS >= 1;if JOBVS = 'V', LDVS >= N.
	doublecomplex wkopt;
	doublecomplex *work;	//(workspace/output) COMPLEX*16 array, dimension (LWORK), On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	integer lwork = -1;		// the size of the workspace, if -1, then a workspace size query is assumed
	doublereal *rwork;		// another workspace
	logical *bwork;

	integer info;			// output parameter: 0 on OK

	//static int bInitialized = 0;
	lda = m_nDim; ldvs = m_nDim;

	// some initialization
	vs = (doublecomplex*) malloc(ldvs*sizeof(doublecomplex));
	rwork = (doublereal*) malloc( m_nDim*sizeof(doublereal) );// don't know why, this is the adviced value
	bwork = (logical*) malloc(m_nDim*sizeof(logical));

	// first call of cheev_- no calcs performed. Just fill lwork so we allocate enough memory in work pointer
	zgees_(&jobvs, &sort, NULL, &m_nDim, m_pMatrix, &lda, &sdim, pVec, vs, &ldvs, &wkopt, &lwork, rwork, bwork, &info);
//	zgees_(char *jobvs, char *sort, L_fp select, integer *n, doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w,
//		doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

	// memory allocation
	lwork = (int) wkopt.r;
	work = (doublecomplex*) malloc( lwork * sizeof(doublecomplex) );

	//printf("%d line work = 0x%x size = %d\n", __LINE__, (unsigned int) work, lwork * sizeof(__CLPK_complex));

	// 2nd call to cheev_- actual computing
	zgees_(&jobvs, &sort, NULL, &m_nDim, m_pMatrix, &lda, &sdim, pVec, vs, &ldvs, work, &lwork, rwork, bwork, &info);

	// release allocated memory
	free(work);
	free(vs);
	free( (void*) rwork);

	return info;
}
