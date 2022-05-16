// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

// CONTAINS:
// 1. nonlinspace
// 2. binary_search
// 3. interp

////////////////////
// 1. NONLINSPACE //
////////////////////

// fills the x-vector with n unequally distributed points between x_min and x_max.
void nonlinespace(double *x, double low, double high, int n, double phi, int *err_flag)
{
    int i;

    x[0] = low;
    for(i = 1; i < n-1; i++){
        x[i] = x[i-1] + (high - x[i-1])/(pow((double)(n-i),phi));
    }
    x[n-1] = high; // numerically exact

        // OUTPUT CHECKS
        #if(CHECKS == 1)
        // monotonicity: x (asc)
        for(i = 1; i < n; i++){
            if(x[i] < x[i-1]){ *err_flag = 1; printf("nonlinspace: x[%d] < x[%d]\n", i, i-1);}
            if(*err_flag == 1) return;
        }
        #endif // CHECKS

    return;

} // nonlinespace

//////////////////////
// 2. BINARY SEARCH //
//////////////////////

// standard binary search algorithm.
// find imin such that x[imin] <= xi <= x[imin+1]
int binary_search(int imin, int num_x, double *x, double *xi, int *err_flag)
{

    int i;
    int imid;
    int imax = num_x - 1;

        // INPUT CHECKS
        #if CHECKS == 1

        // bounds
        if(imin >= imax){       *err_flag = 1; printf("binary_search: imin >= imax\n");}
        if(*xi < x[0]){         *err_flag = 1; printf("binary_search: *xi < x[0]\n");}
        if(*xi >= x[imax]){     *err_flag = 1; printf("binary_search: *xi >= x[imax]\n");}
        if(*err_flag == 1) return -1;

        // monotonicity: x (asc)
        for(i = 1; i < imax; i++){
            if(x[i] < x[i-1]){  *err_flag = 1; printf("binary_search: x[%d] < x[%d]\n",i, i-1);}
            if(*err_flag == 1) return -1;
        }

        if(*err_flag == 1) return -1;
        #endif // CHECKS

    while(x[imin+1] <= *xi){

        imid = (imax + imin)/2;

        if (x[imid] <= *xi){
            imin = imid;
        } else {
            imax = imid;
        }
    }

        // OUTPUT CHECKS
        #if CHECKS == 1
        // bounds
        if(*xi < x[imin]){        *err_flag = 1; printf("binary_search: *xi < x[min]\n");}
        if(*xi >= x[imin+1]){     *err_flag = 1; printf("binary_search: *xi >= x[imin+1]\n");}
        if(*err_flag == 1) return -1;
        #endif // CHECKS

    return imin;

} // binary search

///////////////
// 3. INTERP //
///////////////

// standard linear interpolation.
// extrapolation below: flat
// extrapolation below: linear
void interp(int num_xi, int num_x,
            double *xi, double *yi,
            double *x, double *y, int *err_flag)
{

    int i, imin=0, imax;
    double w;

        // INPUT CHECKS
        #if(CHECKS == 1)
        // bounds
        if(y[0] != 0){ *err_flag = 1; printf("interp: y[0] = %g != 0\n", y[0]);}
        if(*err_flag == 1) return;

        // monotonicity: xi (asc)
        for(i = 1; i < num_xi; i++){
            if(xi[i] < xi[i-1]){ *err_flag = 1; printf("interp: xi[%d] < xi[%d]\n", i, i-1);}
            if(*err_flag == 1) return;
        }

        // monotonicity: x and y (asc)
        for(i = 1; i < num_x; i++){
            if(x[i] < x[i-1]){ *err_flag = 1; printf("interp: x[%d] = %g < %g = x[%d]\n", i, x[i], x[i-1], i-1);}
            if(y[i] < y[i-1]){ *err_flag = 1; printf("interp: y[%d] = %g < %g = y[%d]\n", i, y[i], y[i-1], i-1);}
            if(*err_flag == 1) return;
        }
        #endif // CHECKS

    for(i = 0; i < num_xi; i++){

        // a-I. extrapolation: below (flat)
        if(xi[i] <= x[0]){

            yi[i] = y[0];
            continue; // jumps to the next i

        // a-II. extrapolation: above (linear)
        } else if(xi[i] >= x[num_x-2]) {

            imin = num_x-2;
            imax = num_x-1;

        // a-III. interpolation (linear)
        } else {

            imin = binary_search(0, num_x, &x[0], &xi[i], err_flag);
            imax = imin + 1;

                #if CHECKS == 1
                if(*err_flag == 1) return;
                #endif // CHECKS

        }

        // b. weight
        w = (xi[i] - x[imin]) / (x[imax] - x[imin]);

            #if CHECKS == 1
            if(w < 0){                   *err_flag = 1; printf("interp: w < 0\n");}
            if(w >= 1 & imax < num_x-1){ *err_flag = 1; printf("interp: w >= 1\n");}
            if(*err_flag == 1) return;
            #endif // CHECKS

        // c. result
        yi[i] = y[imin] + w * (y[imax] - y[imin]);

    }

    return;

} // interp
