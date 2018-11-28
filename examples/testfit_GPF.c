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
/**
 * @file TestPDE.c
 * @brief An example program using the PDE solver
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mcl.h"

Dataset data; 

double 
rho(Dataset *data, Dataset *data_s)
{
    double norm;
    double *c_d, *c_s;
    int Nt,Nx;
    norm = 0;

    c_d = (double *)data->fields[3].data_array;
    c_s = (double *)data_s->fields[3].data_array;
    
    Nt = data->fields[3].numRows;
    Nx = data->fields[3].numCols;
    
    for (int i=0;i<Nt*Nx;i++)
    {
        norm += (c_d[i] - c_s[i])*(c_d[i] - c_s[i]);
    }
    
    return norm;
}

/*diffusion function*/
void
D(double * c, int n, double * theta, double* Dc)
{
    double D_0;
    double K;
    double r;
    D_0 = theta[0];
    K = theta[2];
    r = theta[3];
    for (int i=0;i<n;i++)
    {
        Dc[i] = D_0*pow((c[i]/K),r);
    }
}



//#define K 1.7e-3
/*growth function*/
void
f(double *c, int n, double *theta, double *fc)
{
    double lambda;
    double K;
    lambda = theta[1];
    K = theta[2];
    for (int i=0;i<n;i++)
    {
        fc[i] = lambda*c[i]*(1.0 - c[i]/ K);
    }
}

unsigned int failure_counter;


/*simulation function... not stochastic at all, just solves the PDE*/
int
simulate(void *sim, double *theta, Dataset *data_s)
{
    int R = 4; /* the number of gridpoints to be placed between data points*/ 
    int Nx;
    double L,T; /*from data also?*/
    int Nt; /*currently predetermined*/
    double * C_s;
    double * C0;
    double * C_s_f;
    double * C0_f;
    double *Ti;
    double *Xi;
    int nt;
    int nx;

    /* observation times connected with real data.*/
    Ti = (double*)data.fields[0].data_array;
    nt = (int)data.fields[0].numRows;
    Xi = (double*)data.fields[1].data_array;
    nx = (int)data.fields[1].numCols;
    C0 = (double *)data.fields[2].data_array;
    C_s = (double *)data_s->fields[3].data_array;
   
    Nx = nx + R*(nx -1);
    /* need to choose Nt... we could use dtD_0 < dx^2 as a basis */
    /* for now we just pick something that is fine (possibly overkill)*/
    Nt = 32000;
    C0_f = (double *)malloc(Nx*sizeof(double));
    C_s_f = (double *)malloc(Nx*Nt*sizeof(double));

    /*For now a linear interpolation... but I'd like to improve this (e.g., C^1 or C^2)
     * as the initial condition should be differentiable*/
    for (int i=1;i<nx;i++)
    {
        /* Picture: */
        /* data point      i-1                           i        */
        /*            ...---O----------------------------O---...  */
        /* mesh points (i-1)*R  (i-1)R+1 ... (i-1)R+R-1  i*R      */
        /*            ...---O-------O---...------O-------O---...  */
        /* weights a_j      0      1/R        (R-1)/R    1        */

        double C_l,C_u;
        C_l = C0[i-1];
        C_u = C0[i];
        for (int j=0;j<=R+1;j++)
        {
            double alpha;
            alpha = ((double)j)/((double)(R + 1));
            C0_f[(i-1)*R + i-1 + j] = C_l*(1-alpha) + C_u*(alpha); 
        }
    }

    T = Ti[nt-1];
    L = Xi[nx-1]; 
    /*solve using BTCS with fixed-point iterations*/
    if(dbtcsfps(L,T,nt,Ti,theta,C0_f,&D,&f,Nx,Nt,1e-6,10000,C_s_f) < 0)
    {
        failure_counter++;/*just to query the mesh performance*/
    }

    /*print full simulation*/
    fprintf(stdout,"-1");
    for (int i=0;i<Nx;i++)
    {
        double alpha = ((double)i)/((double)Nx);
        fprintf(stdout,",%g",(1.0 - alpha)*Xi[0] + alpha*Xi[nx-1]);
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"0");
    for (int i=0;i<Nx;i++)
    {
        fprintf(stdout,",%g",C0_f[i]);
    }
    fprintf(stdout,"\n");
    for (int j=0;j<nt;j++)
    {
        fprintf(stdout,"%f",Ti[j]);
        for (int i=0;i<Nx;i++)
        {
            fprintf(stdout,",%g",C_s_f[j*Nx + i]);
        }
        fprintf(stdout,"\n");
    }

    /*by construction, data points are nodes in the mesh*/
    for (int j=0;j<nt;j++)
    {
        for(int i=0;i<nx;i++)
        {
            C_s[j*nx + i] = C_s_f[j*Nx + i*(1+R)];
        }
    }
    free(C0_f);
    free(C_s_f);
    sim_counter++;
    return 0;
}

/**
 * Import cell count data *.csv files
 */
void 
ImportCellCountData(char * filename, Dataset* data)
{
    FILE *fp;
    char chr;
    int numCols,numRows;
    double *xi;
    double *ti;
    double *c0;
    double *c;
    double *c_std;
    double *buffer;
    char temp[125];
    fp = fopen(filename,"r");
    numCols = 0;
    numRows = 0;
    while(!feof(fp))
    {
        chr = fgetc(fp);
        numCols += (chr == ',' && numRows == 0);
        numRows += (chr == '\n');
    }
    rewind(fp);
    buffer = (double*)malloc(numRows*(numCols+1)*sizeof(double));
    for (int j=0;j<numRows;j++)
    {
        fscanf(fp,"%lf",buffer + j*(numCols+1));
        for (int i=1;i<numCols+1;i++)
        {
            fscanf(fp,",%lf",buffer +j*(numCols+1) + i);
        }
        fscanf(fp,"\n");
    }
    fclose(fp);
    
    data->numFields = 5; /*time-axis, space-axis, initial condition, time-space counts + std*/
    data->fields = (field *)malloc(data->numFields*sizeof(field));
    /* time-axis*/
    data->fields[0].numRows = numRows-2;
    data->fields[0].numCols = 1;
    data->fields[0].type = REAL64_DATA;
    data->fields[0].numBytes = (data->fields[0].numRows)*sizeof(double);
    data->fields[0].data_array = malloc(data->fields[0].numBytes);
    /* space-axis*/
    data->fields[1].numRows = 1;
    data->fields[1].numCols = numCols/2;
    data->fields[1].type = REAL64_DATA;
    data->fields[1].numBytes = (data->fields[1].numCols)*sizeof(double);
    data->fields[1].data_array = malloc(data->fields[1].numBytes);
    /* initial condition */
    data->fields[2].numRows = 1;
    data->fields[2].numCols = numCols/2;
    data->fields[2].type = REAL64_DATA;
    data->fields[2].numBytes = (data->fields[2].numCols)*sizeof(double);
    data->fields[2].data_array = malloc(data->fields[2].numBytes);
    /* cell count data*/
    data->fields[3].numRows = numRows-2;
    data->fields[3].numCols = numCols/2;
    data->fields[3].type = REAL64_DATA;
    data->fields[3].numBytes = (data->fields[3].numRows)*(data->fields[3].numCols)*sizeof(double);
    data->fields[3].data_array = malloc(data->fields[3].numBytes);
    /*std-deviations (currently not used)*/
    data->fields[4].numRows = numRows - 1;
    data->fields[4].numCols = numCols/2;
    data->fields[4].type = REAL64_DATA;
    data->fields[4].numBytes = (data->fields[4].numRows)*(data->fields[4].numCols)*sizeof(double);
    data->fields[4].data_array = malloc(data->fields[4].numBytes);

    ti = (double*)data->fields[0].data_array;
    for (int i=0;i<data->fields[0].numRows;i++)
    {
        ti[i] = buffer[(i+2)*(numCols+1)];            
    }
    xi = (double*)data->fields[1].data_array;
    for (int i=0;i<data->fields[1].numCols;i++)
    {
        xi[i] = buffer[i + 1];            
    }
    c0 = (double*)data->fields[2].data_array;
    for (int i=0;i<data->fields[2].numCols;i++)
    {
        c0[i] = buffer[numCols + i + 2] ;
    }
    c = (double*)data->fields[3].data_array;
    for (int j=0;j<data->fields[3].numRows;j++)
    {
        for (int i=0;i<data->fields[3].numCols;i++)
        {
            c[j*(data->fields[3].numCols) + i] = buffer[(j+2)*(numCols +1) + i+1];
        }
    }
    free(buffer);
    /*shift spatial domain from [L,U] to [0,U-L]*/
    for (int i=1;i<data->fields[1].numCols;i++)
    {
        xi[i] = xi[i] - xi[0];
    }
    xi[0] = 0.0;
}


int 
main(int argc , char ** argv)
{
    double theta[4];
    double d;
    char* filename;
    double time;
    clock_t start_t,end_t;
    Dataset *data_s;

    if (argc < 6)
    {
        fprintf(stderr,"Usage: %s D lambda K datafile\n",argv[0]);
        exit(1);
    }
    else
    {
        theta[0] = (double)atof(argv[1]);
        theta[1] = (double)atof(argv[2]);
        theta[2] = (double)atof(argv[3]);
        theta[3] = (double)atof(argv[4]);
        filename = argv[5];
    }

    /*import scratch assay data*/
    ImportCellCountData(filename,&data);
    
    SSAL_Initialise(argc,argv);
    
    data_s = copyDataset(&data);
    simulate(NULL,theta,data_s);
    d = rho(&data,data_s);
    fprintf(stderr,"SSE = %g\n",d);
    exit(0); 
}
