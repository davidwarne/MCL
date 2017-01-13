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
 * @brief Emprical estimation of optimal sample numbers for Multi-level 
 * Monte Carlo estimator
 * @datails Uses a small fixed number of trails of the user process at each
 * level. Then contructs a sequence of nl's that is an optimal trade off of 
 * compute effort and accuracy give a target of order eps. 
 *
 * @param L Max level (L+1 levels 0,1,2,...,L)
 * @param n number of trial samples to perfrom at all levels. 
 * @param eps target esimator variance.
 * @param params array of parameters for the stochastic process.
 * @param f function pointer \see dmlmcs
 * @param K dimensionality of the function f ouput.
 * @param nl output array of length L+1 containing the optimal samples numbers.
 *
 * @note n is assumed to be small by comparison to the target accuracy 
 * requirements, but large enough to make a reasonable estimat of vl. As a rule
 * n = 100 is generally ok. If a large discrepancy is observed between the 
 * target estiamtor variance and the actual variance obtained by MLMC 
 * (\see dmlmcs) then n should be increased.
 */
int 
dmlmcnls(int L, int n, double eps, double *params, 
       int (*f)(int, int, int, double *, double*), int K, double *nl)
{
    int i,j,k,l; 
    double Ql[K];
    double fXl[K];
    double vl[L+1];
    double cl[L+1];
    double M;
    clock_t start_t, end_t;

    for (l=0;l<=L;l++)
    {
        /* run and time computation at this level*/
        start_t = clock();
        (*f)(l,n,K,params,fXl);
        dmcintv(n,K,fXl,Ql,&(vl[l]));
        end_t = clock();
        /*record time in seconds, the units actually cancel out, so it does not
         * really matter.*/
        cl[l] = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    }

    /*compute optimal nl according to the Langrange multiplier solution*/
    M = 0;
    for (l=0;l<=L;l++)
    {
        M += sqrt(vl[l]*cl[l]);
    }
    for (l=0;l<=L;l++)
    {
        nl[l] = (sqrt(vl[l]/cl[l])*M)/(eps*eps);
    }
    return 0;
}
