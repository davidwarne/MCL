/* LIBABC: approximate Bayesian Computation
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

/* ABC rejection for inference of Fisher-KPP parameters 
 * using scratch assay data
 */
ABC_Parameters abc_p;
Dataset data_collection;

/*discrepency metric operating on data collections
 *  if \rho(D,D_s) is L2 norm (squared) then
 *   \rho(DC,DC_s) = \sum_j \rho(D_j,D_s,j)
 */
double 
rho(Dataset *data, Dataset *data_s)
{
    double norm;
    double *c_d, *c_s;
    int Nt,Nx,n_d;
    Dataset **datasets_s;
    Dataset **datasets_d;

    norm = 0;
    
    /* cast data arrays to access sub-sets */
    n_d = (data->fields[0].numCols)*(data->fields[0].numRows);
    datasets_s = (Dataset **)data_s->fields[0].data_array;
    datasets_d = (Dataset **)data->fields[0].data_array;

    /* accumulate L2 norm squared over all datasets*/
    for(int j=0;j<n_d;j++)
    {
        /* get data and dimensions*/
        c_d = (double *)datasets_d[j]->fields[3].data_array;
        c_s = (double *)datasets_s[j]->fields[3].data_array;
        Nt = datasets_d[j]->fields[3].numRows;
        Nx = datasets_s[j]->fields[3].numCols;
        /*compute L2 squared*/
        for (int i=0;i<Nt*Nx;i++)
        {
            norm += (c_d[i] - c_s[i])*(c_d[i] - c_s[i]);
        }
    }
    return norm;
}

/*prior sampling process*/
int 
prior(unsigned int dim, unsigned int nsamples, double *support, double *theta)
{
    /* we currently use uninformative priors*/
    for (int i=0;i<nsamples;i++)
    {
        for (int j=0;j<dim;j++)
        {
            theta[i*dim + j] = durngus(support[j],support[j+dim]);
        }
    }
    return 0;
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
        Dc[i] = D_0*pow(c[i]/K,r);
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

    int n_d;
    Dataset  **datasets_s;
    Dataset  **datasets_d;

    n_d = (data_collection.fields[0].numRows)*(data_collection.fields[0].numCols);
    datasets_d = (Dataset **)data_collection.fields[0].data_array;
    datasets_s = (Dataset **)data_s->fields[0].data_array;
    

    for (int j=0;j<n_d;j++)
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
        Ti = (double*)datasets_d[j]->fields[0].data_array;
        nt = (int)datasets_d[j]->fields[0].numRows;
        Xi = (double*)datasets_d[j]->fields[1].data_array;
        nx = (int)datasets_d[j]->fields[1].numCols;
        C0 = (double *)datasets_d[j]->fields[2].data_array;
        C_s = (double *)datasets_s[j]->fields[3].data_array;
   
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
        if(dbtcsfps(L,T,nt,Ti,theta,C0_f,&D,&f,Nx,Nt,1e-6,1000,C_s_f) < 0)
        {
            failure_counter++;/*just to query the mesh performance*/
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
    }
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

/**
 * Read All datasets
 */
void
ReadDatasets(int nfiles, char ** filenames, Dataset* data_collection)
{
    Dataset **datasets;
    /* data collection has a single field with which is a pointer array to
     * n datasets */
    
    /*set up the field*/
    data_collection->numFields = 1;
    data_collection->fields = (field*)malloc(sizeof(field));
    data_collection->fields[0].name[0] = '\0';
    data_collection->fields[0].numRows = 1;
    data_collection->fields[0].numCols = nfiles;
    data_collection->fields[0].numBytes = nfiles*sizeof(Dataset*);
    data_collection->fields[0].type = PTR_DATA;
    data_collection->fields[0].data_array = (void *)malloc(
                                           data_collection->fields[0].numBytes);
    
    /*import the datasets*/
    datasets = (Dataset **)data_collection->fields[0].data_array;
    for (int i=0;i<nfiles;i++)
    {
        datasets[i] = (Dataset *)malloc(sizeof(Dataset));
        ImportCellCountData(filenames[i], datasets[i]); 
    }

}


int 
main(int argc, char** argv)
{
    double *theta;
    double *d;
    clock_t start_t, end_t;
    double time;
    int numFiles;
    char **filenames;
    Dataset **datasets;
    
    /*read user args*/
    if (argc < 3)
    {
        fprintf(stderr,"Usage: %s Nf filename1 ... filenameNf N [sl1,..,slm,su1,...,sum]\n",argv[0]);
        exit(1);
    }
    else
    {
        /* construct ABC parameters*/
        numFiles = atoi(argv[1]);
        filenames = (char**) malloc(numFiles*sizeof(char*));
        //filename = argv[1];
        for (int i=0;i<numFiles;i++)
        {
            filenames[i] = argv[2+i];
        }
        abc_p.nacc = (unsigned int)atoi(argv[2 + numFiles]);
        //abc_p.eps = (double)atof(argv[3]);
        abc_p.nmax = (unsigned int)atoi(argv[2 + numFiles]);
        abc_p.k = 4; /*D, lambda, K*/
        abc_p.support = (double*)malloc(2*abc_p.k*sizeof(double));
        for (int k=0;k<abc_p.k*2;k++)
        {
            abc_p.support[k] = (double)atof(argv[3+numFiles+k]);
        }
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = NULL;
        abc_p.s = &simulate;
    }
    
    /*import scratch assay data*/
    ReadDatasets(numFiles, filenames, &data_collection);
    datasets = (Dataset **)data_collection.fields[0].data_array;
    {
        
        double *x,*t, *c, *c_std,*c0;
        int nx,nt;
        nt =  datasets[0]->fields[0].numRows;
        nx =  datasets[0]->fields[1].numCols;
        t =   (double *) datasets[0]->fields[0].data_array;
        x =   (double *) datasets[0]->fields[1].data_array;
        c0 =  (double *) datasets[0]->fields[2].data_array;
        c =   (double *) datasets[0]->fields[3].data_array;
        fprintf(stderr,"Dataset Information:\n");
        fprintf(stderr,"--------------------\n");
        fprintf(stderr,"Filename: %s\n",filenames[0]);
        fprintf(stderr,"# nodes: %d\n",nx);
        fprintf(stderr,"# timesteps: %d\n",nt);
        fprintf(stderr,"Spatial extent: %lf to %lf\n",x[0],x[nx-1]);
        fprintf(stderr,"Temporal extent: %f to %f\n",t[0],t[nt-1]);
        fprintf(stderr,"X = %f %f ... %f %f %f...%f %f\n",
                 x[0],x[1],x[nx/2 -1], x[nx/2], x[nx/2 +1], x[nx -2], x[nx -1]);
        fprintf(stderr,"T=00.000000 : %f %f ... %f %f %f...%f %f\n", 
                 c0[0], c0[1], c0[nx/2 -1], c0[nx/2], c0[nx/2 +1], c0[nx -2], 
                 c0[nx -1]);
        for (int i=0;i<nt;i++)
        {
            fprintf(stderr,"T=%f : %f %f ... %f %f %f...%f %f\n",
                   t[i],c[i*nx],c[i*nx + 1],c[i*nx + nx/2 -1], c[i*nx + nx/2],
                   c[i*nx + nx/2 +1],c[i*nx + nx -2],c[i*nx +nx -1]);
        }
    }

    /* initialise our RNG library*/
    SSAL_Initialise(argc,argv);

    theta = (double*)malloc(abc_p.nacc*abc_p.k*sizeof(double));
    d = (double*)malloc(abc_p.nacc*sizeof(double));
    sim_counter = 0;

    /*sample using ABC rejection*/
    start_t = clock();
    /*generate prior predictive samples using ABC*/
    dabcrs(abc_p,&data_collection,theta,d);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);

    fprintf(stdout,"\"Sample\",\"D\",\"lambda\",\"K\",\"r\",\"rho\",\"SEC\",\"NACC\",\"NSIMS\",\"NFAIL\"\n");
    for (int i=0;i<abc_p.nacc;i++)
    {
        fprintf(stdout,"%d,%g,%g,%g,%g,%g,%g,%d,%d,%d\n",i,theta[i*abc_p.k],theta[i*abc_p.k+1],theta[i*abc_p.k+2],theta[i*abc_p.k+3],d[i],time,abc_p.nacc,sim_counter,failure_counter);
        fflush(stdout);
    }

    /* clean up the memory mess we made...*/
    free(theta);
    free(d);
    for (int i=0;i<numFiles;i++)
    {
        for (int j=0;j<datasets[i]->numFields;j++)
        {
            free(datasets[i]->fields[j].data_array);
        }
        free(datasets[i]->fields);
        free(datasets[i]);
    }
    free(datasets);
    free(abc_p.support);
    exit(0);
}
