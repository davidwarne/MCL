/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2017  David J. Warne
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "mcl.h"
#include <time.h>

/**
 * @brief Generic Multi-level Monte Carlo estimator
 * @datails Computes an estimator of E[f(X)] using the mulitlevel Monte Carlo 
 * telesoping sum.  
 * bias against exact simulation as a tau-leaping with fine-grain tau.
 * This function computes for a tau leap approximation Z_L with tau = tau0*M^(-L),
 *      E[f(X)] ~ \mu[f(X_L)] = Q_0 + sum_l=1^L Q_l 
 * where
 *      Q_0 = 1/n0 sum_j=1^n0 f(X_0^j) with X_0^j the j-th realisation of X_0
 *      Q_l = 1/nl sum_j=1^nl f(X_l^j) - f(X_{l-1}^j where {X_l^j,X_{l-1}^j} are
 * a coupled pair of realisations at levels l and l-1.
 * @note This function computes the MLMC estimator for user provided Nl, bad 
 * choices in the nl's with result is poor perfromance or an un-reliable) 
 * estimator \see dmlmcnls to obtain empirical estimates of optimal nl if no
 * theoretical results are available.
 *
 * @param L Max level (L+1 levels 0,1,2,...,L)
 * @param nl array of sample numbers for each level 
 * @param params array of parameters for the stochastic process.
 * @param f function pointer to user provide simulation proces, it shall provide
 *          n realisations of the level l function f given parameters params. 
 *          The output is returned via the last argument which is expected to be
 *          a n by K.
 * @param K dimensionality of the function f ouput.
 * @param Ef an array of length K containing the MLMC estimate.
 * @param sigma2 estimate of the estiamtor variance. 
 *
 * @note the user define function must compute both the coarse samples and the
 *       bias differnce samples, as the difference is not explicitly computed. 
 * @todo extend to compute other moments (need to be a bit careful for K > 1)
 * @todo extend to include a convergence check.
 */
int 
dmlmcs(int L, int *nl, double *params, 
       int (*f)(int, int, int, double *, double*), int K, double *Ef, 
       double *sigma2)
{
    int i,j,k,l; 
    double Ql[K];
    double vl;
    double *fXl;
    int max_n;
    
    max_n = nl[0];
    for (l=1;l<=L;l++)
    {
        max_n = (max_n < nl[l]) ? nl[l] : max_n;
    }

    fXl = (double *)malloc(K*max_n*sizeof(double));
    *sigma2 = 0;

    /*compute coarsest level E_0*/
    (*f)(0,nl[0],K,params,fXl);
    dmcintv(nl[0],K,fXl,Ef,&vl);
    *sigma2 = vl/((double)nl[0]);
    printf("E_%d %f\n",0,Ef[0]);
    /* compute bias correction*/
    for (l=1;l<=L;l++)
    {
        (*f)(l,nl[l],K,params,fXl);
        dmcintv(nl[l],K,fXl,Ql,&vl);
        /*compute E_l = E_{l-1} + Q_l*/
        for (k=0;k<K;k++)
        {
            Ef[k] += Ql[k];
        }
        printf("E_%d %f\n",l,Ef[0]);
        *sigma2 += vl/((double)nl[l]);
    }

    return 0;
}
