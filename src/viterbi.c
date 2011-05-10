#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "Rinternals.h"

static void getIndexAndMaxVal(const double *pVec, const int len, double *pMaxVal, int *pMaxIdx)
{
  int i;
  *pMaxVal = *pVec;
  *pMaxIdx = 0;
  for (i=1; i<len; ++i)
  {
    if (pVec[i] > *pMaxVal)
    {
      *pMaxIdx = i;
      *pMaxVal = pVec[i];
    }
  }
}

static void getMatrixIndexAndMaxVal(const double *pMat, const int nCols, double *pMaxVal, int *pMaxIdx, const int nRows)
{
  int i;
  *pMaxVal = *pMat;
  *pMaxIdx = 0;
  for (i=1; i<nCols; ++i)
  {
    if (pMat[i*nRows] > *pMaxVal)
    {
      *pMaxIdx = i;
      *pMaxVal = pMat[i*nRows];
    }
  }
w}

/**
 * viterbi
 * \param pBeta - log emission probabilities SxT
 * \param initialP - log initial state probabilities - length S
 * \param tau - transition probabilities - original scale
 * \param pArm - indicator for chromosome arm
 * \param S - number of columns of beta (number of states)
 * \param T - number of rows of beta matrix
 * \param pQHat - output vector of integers
 * \param pDelta - output vector of doubles ...this is the forward variable
 * \param c1  Pr( Normal -> altered)
 * \param c2  Pr( altered -> normal)
 * \param c3  Pr( altered -> altered)
 * \param normalState  index of normal state
 * \param pAA
 */
void viterbi(double *pBeta, double *initialP, double *tau,
             int  *pArm, int *S, int *T, int *pQHat, double *pDelta,
             double *c1, double*c2, double *c3, int *normalState, double *pAA)
{
  /**  RS double *pDelta, *pAA, *pDeltaTempSum, Pstar; */
  /** double *pAA, *pDeltaTempSum, Pstar, *tp; */
  double *pDeltaTempSum, Pstar, *tp;
  int i,j,t;
  int nRows, nCols, *pPsi;
  int NS;


  NS = *normalState - 1;
  nRows = *T;
  nCols = *S;

  /** RS
      pDelta = (double *)R_alloc(sizeof(double), nRows * nCols); */
  pPsi = (int *)R_alloc(sizeof(int), nRows * nCols);
  /** transition probability matrix */
  /** pAA = (double *)R_alloc(sizeof(double), nCols * nCols); */
  pDeltaTempSum = (double *)R_alloc(sizeof(double), nCols);

  /** what is this notation? *(pDelta + nrows*j) C uses vectors and
      not matrices.  so the Nth row in pDelta is set to initialP + the
      emission probability for the Nth row*/
  for (j=0; j<nCols; ++j)
  {
    *(pDelta + nRows*j) = initialP[j] + *(pBeta + nRows*j);
    *(pPsi + nRows*j) = 0;
  }

  for (t=1; t<nRows; ++t)
    {
      /*      if (strcmp(*(ppArm + t), *(ppArm + t - 1)) != 0)*/
      if(*(pArm + t) != *(pArm + t - 1))
	{
	  for (j=0; j<nCols; ++j)
	    {
	      *(pDelta + j*nRows + t) = initialP[j] + *(pBeta + j*nRows + t);
	      *(pPsi + j*nRows + t) = 0;
	    }
	  continue;
	}
      for (i=0; i<nCols; ++i)
	{
	  for (j=0; j<nCols; ++j)
	    {
	      int offset;
	      offset = j * nCols + i;
	      /* if (i == j)*/
	      if(i == NS)
		{
		  if(i == j)  /* probability of staying in the normal state */
		    {
		      *(pAA + offset) = 1 - ((1-tau[t-1]) * (nCols - 1) * *c1);
		      /* printf(i, j); */
		      /* probability of staying in the same state */
		      /* *(pAA + offset) = tau[t-1]; */
		    }
		  else /* probability of leaving normal state */
		    {
		      *(pAA + offset) = *c1 * (1-tau[t-1]);
		    }
		}
	      else   /* transitioning from an altered state */
		{
		  if(i == j)  /* staying in the same altered state */
		    {
		      /* c2 = scalar for transitioning from normal to altered state */
		      *(pAA + offset) = 1 - (1 - tau[t-1]) * (*c2 + (nCols - 2) * *c3);
		      /* *(pAA + offset) = (1-tau[t-1])/(nCols-1); */
		    }
		  else /* leaving altered state */
		    {
		      if(j == NS) /* going back to normal state */
			{
			  *(pAA + offset) = *c2 * (1 - tau[t-1]);
			}
		      else  /* going to another altered state */
			{
			  *(pAA + offset) = *c3 * (1 - tau[t-1]);
			}
		    }
		}
	      /* *(pAA + offset) = log ( *(pAA + offset) * *(tau_scale + offset) );*/
	      *(pAA + offset) = log ( *(pAA + offset) );
	    }
	}
      for (j=0; j<nCols; ++j)
	{
	  double maxDeltaTempSum;
	  int maxDeltaSumIdx = 0;
	  /* Sum the jth column of AA and the (t-1)th row of delta.
             The jth column of AA occupies
	     consecutive memory locations starting at the memory location pAA + nrow(AA)*j.
	     Since AA is a square matrix, that starting address can be
	     expressed as pAA + nCols*j */
	  for (i=0; i<nCols; ++i)
	    {
	      pDeltaTempSum[i] = pAA[j * nCols + i] + pDelta[(t-1) + i * nRows];
	    }
	  /* Needs update */
	  getIndexAndMaxVal( (double *)(pDeltaTempSum), nCols, &maxDeltaTempSum, &maxDeltaSumIdx);
	  *(pPsi + j * nRows  + t) = maxDeltaSumIdx;
	  *(pDelta + j*nRows + t) = maxDeltaTempSum + *(pBeta + j*nRows + t);
	}
    }

  /* Needs update */
  getMatrixIndexAndMaxVal( (double *)(pDelta + nRows-1), nCols, &Pstar, (int *)(pQHat + nRows-1), nRows);
  for (t=nRows-2; t>= 0; --t)
    {
      /*if (strcmp(*(ppArm + t), *(ppArm + t + 1)) != 0)*/
      if(*(pArm + t) != *(pArm + t + 1))
	{
	  double maxVal;

	  /* Needs update */
	  getMatrixIndexAndMaxVal((double *)(pDelta + t), nCols, &maxVal, &pQHat[t], nRows);
	}
      else
	{
	  pQHat[t] = *(pPsi + pQHat[t+1] * nRows + (t+1));
	}
    }

  /* Array indices in R are 1-based so add one to these values. */
  for (i=0; i<nRows; ++i) {
    pQHat[i] += 1;
    if (i > 0)
      for (j=0; j<nCols; ++j)
        pPsi[j*nRows + i] += 1;}
}
