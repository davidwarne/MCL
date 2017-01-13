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
#include <time.h>
/**
 * @brief samples variance to esitmate number of required samples
 */
int 
dabcns(ABC_Parameters abc_p,int trials, double error, Dataset *data, 
       CDF_estimate *M, unsigned int *N)
{
    int l,j,i,k;
    double times;
    double *V;
    double sigma2;
    double *E;
    double *E2;
    double *theta;
    double c;
    clock_t start_t,end_t;

    abc_p.nacc = trials;

    if ((theta = (double*)malloc(trials*abc_p.k*sizeof(double))) == NULL) 
    {
        return 1;
    }

    /*begin timing trials*/
    start_t = clock();

    /* sample ABC */
    dabcrs(abc_p, data,theta, NULL);
    for (j=0;j<M->G.numPoints;j++)
    {
        dmcint(abc_p.nacc,abc_p.k,theta,M->G.coords +j*M->G.dim,M->g,M->F+j,
               M->V+j);
    }
    end_t = clock();
    sigma2 = M->V[0];
    for (j=0;j<M->G.numPoints;j++)
    {
        if (sigma2 < M->V[j])
        {
            sigma2 = M->V[j];
        }
    }
    times = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);

    *N = (unsigned int)ceil(log((double)M->G.numPoints)*sigma2/(error*error));
    free(theta);
    return 0;
}
