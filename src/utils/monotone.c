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
 * Replaces the estimate CDF with a monotonic version.
 * This cannot increase the L_infty error.
 *
 * Given continuous function F(x_1,x_2,...,x_n) we take
 *
 * Fm(s) = 1/2 (sup_{A_s} F + inf_{B_s} F)
 *
 * where A_s = {x : x_1 < s_1, x_2 < s_2, ... ,x_n < s_n}
 *
 */
void monotone(CDF_estimate *M, CDF_estimate *M_mon)
{
    int j,j_1,k;
    double *F_sup;
    double *F_inf;
    unsigned int x[255],ibits;    
    double F_min,F_max;
    j_1 = 1;
    /* j_1 = (((((D+1)*D + 1)*D + 1)*D + 1)....)*D + 1*/

    F_sup = (double*)malloc(M->G.numPoints*sizeof(double));
    F_inf = (double*)malloc(M->G.numPoints*sizeof(double));

    memcpy(F_sup,M->F,M->G.numPoints*sizeof(double));
    memcpy(F_inf,M->F,M->G.numPoints*sizeof(double));

    for (j=0;j<M->G.numPoints;j++)
    {
        for (ibits=1;ibits < (1 << M->G.dim);ibits++){
            IND2SUB(j,M->G.dim,M->G.D,x)
            unsigned int _ibits = ibits;
            for (k=0;k<M->G.dim;k++){
                x[k] -= (_ibits % 2)*(x[k] > 0);
                _ibits = _ibits >> 1;
            }
            SUB2IND(x,M->G.dim,M->G.D,j_1)
            F_sup[j] = (F_sup[j] < F_sup[j_1]) ? F_sup[j_1] : F_sup[j];
        }
    }
    for (j=M->G.numPoints-1;j>=0;j--)
    {
        for (ibits=1;ibits < (1 << M->G.dim);ibits++){
            unsigned int _ibits = ibits;
            IND2SUB(j,M->G.dim,M->G.D,x)
            for (k=0;k<M->G.dim;k++){
                x[k] += (_ibits % 2)*(x[k] < (M->G.D-1));
                _ibits = _ibits >> 1;
            }
            SUB2IND(x,M->G.dim,M->G.D,j_1)
            F_inf[j] = (F_inf[j] > F_inf[j_1]) ? F_inf[j_1] : F_inf[j];
        }
    }

    //F_min = 1;
    F_max = 0;
    for (j=0;j<M->G.numPoints;j++)
    {
        M_mon->F[j] = 0.5*F_sup[j] + 0.5*F_inf[j];
   //     M_mon->F[j] = F_sup[j] ;
//        F_min = (M_mon->F[j] < F_min) ? M_mon->F[j] : F_min;
    }
//    if (F_max > 1.0 && F_min < 0)
//    {
        for (j=0;j<M->G.numPoints;j++)
        {
            M_mon->F[j] = (M_mon->F[j] < 0) ? 0 : M_mon->F[j] ;
        }
 //   }
    for (j=0;j<M->G.numPoints;j++)
    {
        F_max = (M_mon->F[j] > F_max) ? M_mon->F[j] : F_max;
    }
    //if (F_max > 1.0)
   // {
        for (j=0;j<M->G.numPoints;j++)
        {
        
            M_mon->F[j] /= F_max;
        }
   // }

    F_max = 0;
    for (j=0;j<M->G.numPoints;j++)
    {
        F_max = (M_mon->F[j] > F_max) ? M_mon->F[j] : F_max;
    }

    free(F_sup);
    free(F_inf);
}
