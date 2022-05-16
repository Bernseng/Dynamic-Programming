// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

// fills c_vec and nb_vec (i.e. the unconstrained consumption function)
void ucon_c_func(int i_x_now, double d, double *n_vec, double *c_vec, double *nb_vec,
                 double *nb_plus, double *c_plus,
                 int *i_nb_plus_f, int *Nnb_plus_db,
                 int *i_db_plus_exp, double *nb_plus_exp, double *c_plus_exp,
                 double *psi, double *xit, double *w,
                 double beta, double rho, double Gamma, double ra, double rd,
                 int Nx, int Ndb, int Nnb_max, int Npsi, int  Nxit, int Nn,
                 int *err_flag, long int *t_vec)
{

    int i_n, i_x, i_x_loop, i_psi, i_xit, i_db_plus, i_rho;
    long int index_plus_w, index_plus_u, index_exp;
    double sum_temp, temp, temp_fac;
    clock_t start_val, end_val;

    ////////////////////////////////////
    // 1. c = 0 at state space border //
    ////////////////////////////////////

    c_vec[0] = 0;
    nb_vec[0] = n_vec[0];


    ///////////////////////////////////////////
    // 2. loop over n_vec -> c_vec -> nb_vec //
    ///////////////////////////////////////////

    for(i_n = 1; i_n < Nn; i_n++){

        ///////////////////////
        // 3. next-period nb //
        ///////////////////////

        nb_plus_func(d, n_vec[i_n], &nb_plus_exp[0], &psi[0], &xit[0], Npsi, Nxit, Gamma, ra, rd, err_flag);

            // shape(n_vec) = (Nn,)

            #if CHECKS == 1
            if(*err_flag == 1){return;}
            #endif // CHECKS


        ///////////////////////////////
        // 4. next-period c (vector) //
        ///////////////////////////////

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

            // b. interpolation of next-period c

                    #if TIMING == 1
                    start_val = clock();
                    #endif // TIMING

            index_exp = i_x_loop*Npsi*Nxit + i_psi*Nxit;

            // unemployed:
            interp(1, Nnb_plus_db[1*Nx*Ndb + i_x*Ndb + i_db_plus],
                   &nb_plus_exp[i_psi*Nxit + 0], &c_plus_exp[index_exp + 0],
                   &nb_plus[index_plus_u], &c_plus[index_plus_u], err_flag);

                    // shape(Nnb_plus_db) = (Nu,Nx,Ndb)
                    // shape(nb_plus_exp) = (Npsi,Nxit)
                    // shape(c_plus_exp) = (Nx,Npsi,Nxit)
                    // shape(nb_plus) = shape(c_plus) = (Nu,Nx,Ndb,Nnb_max)

            // working:
            interp(Nxit-1, Nnb_plus_db[0*Nx*Ndb + i_x*Ndb + i_db_plus],
                   &nb_plus_exp[i_psi*Nxit + 1], &c_plus_exp[index_exp + 1],
                   &nb_plus[index_plus_w], &c_plus[index_plus_w], err_flag);

                    // shape(Nnb_plus_db) = (Nu,Nx,Ndb)
                    // shape(nb_plus_exp) = (Npsi,Nxit)
                    // shape(c_plus_exp) = (Nx,Npsi,Nxit)
                    // shape(nb_plus) = shape(c_plus) = (Nu,Nx,Ndb,Nnb_max)

                    #if TIMING == 1
                    end_val = clock();
                    t_vec[2] += (end_val - start_val)*1000/CLOCKS_PER_SEC;
                    #endif // TIMING

                    #if CHECKS == 1
                    if(*err_flag == 1){return;}
                    #endif // CHECKS

        } // psi
        } // x


        ////////////////////////////
        // 5. next-period c (sum) //
        ////////////////////////////

            #if TIMING == 1
            start_val = clock();
            #endif // TIMING

        sum_temp = 0.0;
        for(i_x_loop = 0; i_x_loop < Nx; i_x_loop++){
        for(i_psi = 0; i_psi < Npsi; i_psi++){
        for(i_xit = 0; i_xit < Nxit; i_xit++){

            index_exp = i_x_loop*Npsi*Nxit + i_psi*Nxit + i_xit;

            /*
            // only works if rho is an integer:
            temp_fac = Gamma*psi[i_psi]*c_plus_exp[index_exp];
            temp = temp_fac;
            for(i_rho = 1; i_rho < (int)rho; i_rho++){
                temp *= temp_fac;
            }
            sum_temp += w[index_exp]/temp;
            */

            // is equivalent to:
            sum_temp += w[index_exp]*pow(Gamma*psi[i_psi]*c_plus_exp[index_exp], -rho);

                // shape(w) = (Nx,Npsi,Nxit) (because of reduction in function call)
                // shape(c_plus_exp) = (Nx,Npsi,Nxit)

        } } }

            #if TIMING == 1
            end_val = clock();
            t_vec[3] += (end_val-start_val)*1000/CLOCKS_PER_SEC;
            #endif // TIMING

        //////////////////////////////////
        // 6. construct c_vec -> nb_vec //
        //////////////////////////////////

        sum_temp    *= (1.0+ra)*beta;

        c_vec[i_n]  = pow(sum_temp,-1.0/rho);
        nb_vec[i_n] = n_vec[i_n] + c_vec[i_n];

            // shape(c_vec) = (Nn,)
            // shape(n_vec) = (Nn,)
            // shape(nb_vec) = (Nn,)

    } // n_vec

		////////////////////////////
		// IMPORTANT CHECK OF EGM //
		////////////////////////////

		for(i_n = 1; i_n < Nn; i_n++){
            if(nb_vec[i_n] <= nb_vec[i_n-1]){
                *err_flag = 1;
                printf("ucon_c_func: nb_vec[%d] = %g < %g = nb_vec[%d] (diff_nb = %g, diff_n = %g, diff_c = %g)\n", i_n, nb_vec[i_n], nb_vec[i_n-1], i_n-1, nb_vec[i_n]-nb_vec[i_n-1], n_vec[i_n]-n_vec[i_n-1], c_vec[i_n]-c_vec[i_n-1]);
            }
            if(*err_flag == 1) return;
        }

        // OUTPUT CHECKS
        #if(CHECKS == 1)
        // monotonicity: c_vec (asc)
        for(i_n = 1; i_n < Nn; i_n++){
            if(c_vec[i_n] < c_vec[i_n-1]){
                *err_flag = 1;
                printf("ucon_c_func: c_vec[%d] < c_vec[%d]\n", i_n, i_n-1);
            }
            if(*err_flag == 1) return;
        }
        #endif // CHECKS

    return;

} // ucon_c_func
