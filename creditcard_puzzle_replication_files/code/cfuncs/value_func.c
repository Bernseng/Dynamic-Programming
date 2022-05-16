// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

// calculates the value of choice (d, c|n)
double value_func(int i_x_now, double d, double n, double c,
                  double *db_plus, double *nb_plus, int *i_nb_plus_f, int *Nnb_plus_db, double *vt_plus,
                  int *i_db_plus_exp, double *nb_plus_exp, double *vt_plus_exp,
                  double *psi, double *xit, double *w_scale,
                  double beta, double rho, double Gamma, double lambda, double ra, double rd,
                  int Nx, int Ndb, int Nnb_max, int Npsi, int  Nxit,
                  int *err_flag, long int *t_vec)
{

    int i_x, i_x_loop, i_psi, i_xit, i_db_plus;
    long int index_plus_w, index_plus_u, index_exp;
    double v, cont_v, u;
    clock_t start_val, end_val;

    //////////////////////////////
    // 1. next-period db and nb //
    //////////////////////////////

    i_db_plus_exp_func(d, &db_plus[0], &i_db_plus_exp[0],
                       &psi[0], Gamma, lambda,
                       Ndb, Npsi, err_flag);

    nb_plus_func(d, n, &nb_plus_exp[0], &psi[0], &xit[0], Npsi, Nxit, Gamma, ra, rd, err_flag);

        // shape(nb_plus_exp) = (Npsi, Nxit)

        #if CHECKS == 1
        if(*err_flag == 1) return 1.0;
        #endif // CHECKS

    ////////////////////////////////////
    // 2. continuation value (vector) //
    ////////////////////////////////////

    for(i_x_loop = 0; i_x_loop < Nx; i_x_loop++){
    for(i_psi = 0; i_psi < Npsi; i_psi++){

        if(d > 0 | i_x_now == 1){
            i_x = i_x_loop;
        } else {
            i_x = 0;
        }

        // a. next-period db (index)
        i_db_plus = i_db_plus_exp[i_psi];

            // working: next-period state index
            index_plus_w = 0*Nx*Ndb*Nnb_max
                           + i_x*Ndb*Nnb_max
                           + i_db_plus*Nnb_max
                           + i_nb_plus_f[0*Nx*Ndb + i_x*Ndb + i_db_plus];

            // unemployed: next-period state index
            index_plus_u = 1*Nx*Ndb*Nnb_max
                           + i_x*Ndb*Nnb_max
                           + i_db_plus*Nnb_max
                           + i_nb_plus_f[1*Nx*Ndb + i_x*Ndb + i_db_plus];

                    // shape(i_nb_plus_f) = (Nu,Nx,Ndb)

        // b. check for escape cases
        if(nb_plus_exp[i_psi*Nxit + 0] <= nb_plus[index_plus_u]){

            return -INFINITY;

        } else if(nb_plus_exp[i_psi*Nxit + 1] <= nb_plus[index_plus_w]){

            return -INFINITY;

        }

        // c. interpolation of vt

                #if TIMING == 1
                start_val = clock();
                #endif // TIMING

            index_exp = i_x_loop*Npsi*Nxit + i_psi*Nxit;

            // unemployed:
            interp(1, Nnb_plus_db[1*Nx*Ndb + i_x*Ndb + i_db_plus],
                   &nb_plus_exp[i_psi*Nxit + 0], &vt_plus_exp[index_exp + 0],
                   &nb_plus[index_plus_u], &vt_plus[index_plus_u], err_flag);

            // working:
            interp(Nxit-1, Nnb_plus_db[0*Nx*Ndb + i_x*Ndb + i_db_plus],
                   &nb_plus_exp[i_psi*Nxit + 1], &vt_plus_exp[index_exp + 1],
                   &nb_plus[index_plus_w], &vt_plus[index_plus_w], err_flag);

                        // shape(Nnb_plus_db) = (Nu,Nx,Ndb)
                        // shape(nb_plus_exp) = (Npsi,Nxit)
                        // shape(vt_plus_exp) = (Nx,Npsi,Nxit)
                        // shape(nb_plus) = shape(vt_plus) = (Nu,Nx,Ndb,Nnb_max)

                #if TIMING == 1
                end_val = clock();
                t_vec[6] += (end_val - start_val)*1000/CLOCKS_PER_SEC;
                #endif // TIMING

                #if CHECKS == 1
                if(*err_flag == 1) return 1.0;

                // bounds
                if(vt_plus_exp[index_exp] < 0){ *err_flag = 1; printf("vt_plus_exp[0] < 0\n");}

                // monotonicity
                for(i_xit = 1; i_xit < Nxit; i_xit++){
                    if(vt_plus_exp[index_exp+i_xit] <= vt_plus_exp[index_exp+i_xit-1]){
                        *err_flag = 1;
                        printf("vt_plus_exp[%d] = %g <= %g = vt_plus_exp[%d] (%g %g)\n",
                                index_exp+i_xit, vt_plus_exp[index_exp+i_xit],
                                vt_plus_exp[index_exp+i_xit-1], index_exp+i_xit-1,
                                nb_plus_exp[i_psi*Nxit + i_xit], db_plus[i_db_plus]);
                    }
                }
                if(*err_flag == 1) return 1.0;
                #endif // CHECKS

    } // psi
    } // x


    /////////////////////////////////
    // 3. continuation value (sum) //
    /////////////////////////////////

    cont_v = 0.0;
    for(i_x_loop = 0; i_x_loop < Nx; i_x_loop++){
    for(i_psi = 0; i_psi < Npsi; i_psi++){
    for(i_xit = 0; i_xit < Nxit; i_xit++){

        index_exp = i_x_loop*Npsi*Nxit + i_psi*Nxit + i_xit;

        cont_v += w_scale[index_exp]*(-1.0/vt_plus_exp[index_exp]);

            // shape(w_scale) = (Nx,Npsi,Nxit) (because of reduction in function call)
            // shape(vt_plus_exp) = (Nx,Npsi,Nxit)

    } } }


    ///////////////////////////////
    // 4. utility of consumption //
    ///////////////////////////////

    u = pow(c, 1.0-rho)/(1.0-rho);


    //////////////////////////////
    // 5. total value of choice //
    //////////////////////////////

    v = u + beta*cont_v;

        // OUTPUT CHECKS
        #if(CHECKS == 1)
        if(cont_v < -10000000000000000000000000000000000000000000000.0){
            *err_flag = 1;
            printf("(cont_v  too low, %g %g %g %g\n", d, n, c, cont_v);
        }
        if(cont_v >= 0){ *err_flag = 1; printf("cont_v >= 0\n");}
        if(isfinite(cont_v) == 0){ *err_flag = 1; printf("cont_v %g is not finite\n", cont_v);}

        if(v >= 0){ *err_flag = 1; printf("v >= 0\n");}
        if(isfinite(v) == 0){ *err_flag = 1; printf("v %g is not finite\n", v);}

        if(*err_flag == 1) return 1.0;
        #endif // CHECKS

    return v;

}
