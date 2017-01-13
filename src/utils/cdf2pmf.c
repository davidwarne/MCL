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
 * Computes numerically the multivaruate PMF from the CDF using
 * finite differences.
 */
double 
cdf2pmf(CDF_estimate *M, double *p)
{
    int j,j_1,k,jp1,jm1;
    unsigned int x[255];    
    double *p_temp;
    double s;
    p_temp = (double*)malloc(M->G.numPoints*sizeof(double));
    memcpy(p_temp,M->F,M->G.numPoints*sizeof(double));
    for (k=0;k<M->G.dim;k++)
    {
        for (j=0;j<M->G.numPoints;j++)
        {
            IND2SUB(j,M->G.dim,M->G.D,x)
            if (x[k] > 0 && x[k] < M->G.D-1)
            {
                x[k] -= 1;
                SUB2IND(x,M->G.dim,M->G.D,jm1)
                x[k] +=2 ;
                SUB2IND(x,M->G.dim,M->G.D,jp1)
                /*central difference*/
                p[j] = (p_temp[jp1] - p_temp[jm1])/2.0;
            } 
            else if (x[k] == 0 )
            {
                /*forward difference*/
                x[k] += 1 ;
                SUB2IND(x,M->G.dim,M->G.D,jp1)
                p[j] = 0;
            } 
            else if (x[k] == M->G.D-1)
            {
                /*backward differnce*/
                x[k] -= 1 ;
                SUB2IND(x,M->G.dim,M->G.D,jm1)
                p[j] = 0;
            }
        }
        memcpy(p_temp,p,M->G.numPoints*sizeof(double));
    }

    /*test*/
    s = 0;
    for (j=0;j<M->G.numPoints;j++)
    {
        s += p[j];
    }
    fprintf(stderr,"%f\n",s);
    free(p_temp);
    return s;
}
