/* MCL: Monte Carlo Library.
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
/**
 * @brief Monte Carlo Integration for a random vector
 * @details Approximates E[X]  and V[X] = E[|X - E[X]|^2_\infty]
 *  X is a random vectors in R^k.
 *
 * @param N the number of iid samples.
 * @param k dimensionality of X random vector.
 * @param X an Nxk array of iid samples. 
 * @param E output of the expectation.
 * @param V output of the variance.
 */
int 
dmcintv(unsigned int N,unsigned int k, double *X, double *Y, double *params, 
        double (*f)(int, double *,double*), double *E, double *V)
{
    unsigned int i;
    double fX_i,fY_i;
    
    *E = 0;
    *V = 0;
    for (i=0;i<N;i++)
    {
        fX_i = (*f)(k,X + i*k,params);
        *E += fX_i - fY_i;
        *V += (fX_i - fY_i)*(fX_i - fY_i);
    }
    
    *E /= (double)N;
    *V /= (double)N;
    *V = ((double)N)*(*V - (*E)*(*E))/((double) N-1);
    //printf("%d %f %f\n",N,*E,*V);
}
