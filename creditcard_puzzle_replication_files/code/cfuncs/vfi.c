// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// number of threads in parallel region
#define MAXTHREADS MIN(24, omp_get_max_threads()-2)

#define CHECKSUM 1 // calculate a checksum
#define TIMING_TOT 1 // time overall
#define VTICHECK 0 // compare with value function iteration

// should not be used in parallel
#define CHECKS 0 // detailed checks
#define TIMING 0 // time in detail

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

// note: the order is important
#include "generic_funcs.c"
#include "base_funcs.c"
#include "ucon_c_func.c"
#include "value_func.c"
#include "grid_search.c"

void main(){}

void vfi(double *v_T, double *d_T, double *Delta_T, double *c_T, double *n_T, double *a_T, double *vt_T, double *ucon_c_d_T,
         double *db, double *nb_T, double *kappa, int *i_nb_f_T, int *Nnb_db_T,
         double *psi, double *xit, double *w, double *w_scale,
         double beta, double rho, double Gamma,
         double ra, double rd, double lambda, double eta, double varphi, double chi,
         double n_max, double phi_n, double eps,
         int T, int tcon, int Nu, int Nx, int Ndb, int Nnb_max, int Npsi, int  Nxit, int Nn, long int Nd, double d_stepsize)
{

    int i, *err_flag;
    long int t_tot, *t_vec;
    clock_t start, end, start_val, end_val, start_val_inner, end_val_inner;
    FILE *log_txt;
    double vfi_check[1];
    vfi_check[1] = 0.0;

    ////////////////////////////
    // allocate common memory //
    ////////////////////////////

    // error flag (remember to free)
    err_flag    = (int*) calloc(1, sizeof(int));
    err_flag[0] = 0;

        #if TIMING == 1
        t_vec = (long int*) calloc(20, sizeof(long int));
        for(i = 0; i < 20; i++){t_vec[i] = 0.0;}
        start = clock();
        #endif // TIMING

        #if TIMING_TOT == 1 && TIMING == 0
        start = clock();
        #endif // TIMING


    ////////////////////
    // begin parallel //
    ////////////////////

    #pragma omp parallel num_threads(MAXTHREADS)
    {

    // private variables for each thread
    int t, i_n, i_u, i_x, i_db, i_nb;
    int i_db_lower, i_db_upper, i_db_mid;
    long int index, index_db_min, index_db_max;
    long int *i_d_max, i_d, i_d_upper, i_d_lower, i_d_db_min;

    double *v, *d, *Delta, *c, *n, *a, *vt, *ucon_c_d, *nb;
    double *v_plus, *db_plus, *nb_plus, *c_plus, *vt_plus, *kappa_plus;
    int *i_nb_f, *Nnb_db, *i_nb_plus_f, *Nnb_plus_db;

    double *n_vec, *c_vec, *nb_vec;
    int *i_db_plus_exp;
    double *nb_plus_exp, *c_plus_exp, *vt_plus_exp;

    double d_now, varphi_now;
    double chekksum=0.0, max_diff, vfi_check_now, temp;
    FILE *log_txt;


    // allocate memory (remember to free)
    n_vec           = (double*) calloc(Nn, sizeof(double));
    c_vec           = (double*) calloc(Nn, sizeof(double));
    nb_vec          = (double*) calloc(Nn, sizeof(double));
    i_d_max         = (long int*) calloc(1, sizeof(long int));
    i_db_plus_exp   = (int*) calloc(Npsi, sizeof(int));
    nb_plus_exp     = (double*) calloc(Npsi*Nxit, sizeof(double));
    c_plus_exp      = (double*) calloc(Nx*Npsi*Nxit, sizeof(double));
    vt_plus_exp     = (double*) calloc(Nx*Npsi*Nxit, sizeof(double));


    ////////////////////////////////
    // 1. terminal value function //
    ////////////////////////////////

    #pragma omp single
    {

    for(i_u  = 0; i_u  < Nu;  i_u++){
    for(i_x  = 0; i_x  < Nx;  i_x++){
    for(i_db = 0; i_db < Ndb; i_db++){
    for(i_nb = 0; i_nb < Nnb_db_T[(T-1)*Nu*Nx*Ndb + i_u*Nx*Ndb + i_x*Ndb + i_db]; i_nb++){

        index = (T-1)*Nu*Nx*Ndb*Nnb_max
                + i_u*Nx*Ndb*Nnb_max
                + i_x*Ndb*Nnb_max
                + i_db*Nnb_max
                + i_nb_f_T[(T-1)*Nx*Nu*Ndb + i_u*Nx*Ndb + i_x*Ndb + i_db] + i_nb;

        c_T[index]     = nb_T[index];
        d_T[index]     = 0;
        Delta_T[index] = 0 - db[i_db];
        n_T[index]     = nb_T[index] - c_T[index];
        a_T[index]     = 0;

        temp = vt_T[index];
        if(i_nb == 0){

            v_T[index] = -INFINITY;
            vt_T[index] = 0.0;

        } else {

            v_T[index]     = pow(c_T[index], 1.0-rho)/(1.0-rho);
            vt_T[index]    = -1.0/v_T[index];
        }

    } } } }

    }
    #pragma omp barrier

    //////////////////
    // 2. time loop //
    //////////////////

    #pragma omp master
    {
    log_txt = fopen("log.txt", "w");
    if(log_txt == NULL){exit(-1); }
    fclose(log_txt);
    }

    for(t = T-2; t >= 0; t--){

        index = t*Nu*Nx*Ndb*Nnb_max;

            v       = &v_T[index];
            d       = &d_T[index];
            Delta   = &Delta_T[index];
            c       = &c_T[index];
            n       = &n_T[index];
            a       = &a_T[index];
            vt      = &vt_T[index];
            nb      = &nb_T[index];

            ucon_c_d = &ucon_c_d_T[0*Nu*Nx*Nd*Nnb_max]; // note: Nd instead of Ndb

        index = (t+1)*Nu*Nx*Ndb*Nnb_max; // note: t+1

            v_plus  = &v_T[index];
            db_plus = &db[0];
            nb_plus = &nb_T[index];
            vt_plus = &vt_T[index];
            c_plus  = &c_T[index];

        index = t*Nu*Nx*Ndb;

            i_nb_f = &i_nb_f_T[index];
            Nnb_db = &Nnb_db_T[index];

        index = (t+1)*Nu*Nx*Ndb; // note: t+1

            i_nb_plus_f = &i_nb_f_T[index];
            Nnb_plus_db = &Nnb_db_T[index];
            kappa_plus  = &kappa[index];


    //////////////////////////////////
    // 3. unconstrained consumption //
    //////////////////////////////////

        #if TIMING == 1
        start_val = clock();
        #endif // TIMING

    #pragma omp for schedule(static)
    for(i_d = 0; i_d < Nd; i_d++){

        // a. d choice
        d_now = (double)i_d*d_stepsize;

        // b. n vector and db plus index vector
        i_db_plus_exp_func(d_now, &db_plus[0], &i_db_plus_exp[0],
                           &psi[0], Gamma, lambda,
                           Ndb, Npsi, err_flag);

        n_vec_func(d_now, &n_vec[0], &kappa_plus[0], &i_db_plus_exp[0],
                   &psi[0], &xit[0], Gamma, ra, rd,
                   n_max, phi_n, eps,
                   Nu, Nx, Ndb, Npsi, Nn, err_flag);

            #if(CHECKS == 1)
            if(*err_flag == 1){ printf("\nEXIT_FAILURE (n_vec_func)\n\n"); exit(EXIT_FAILURE);}
            #endif

        // c. unconstrained consumption function
        for(i_u = 0; i_u < Nu; i_u++){
        for(i_x = 0; i_x < Nx; i_x++){

                #if TIMING == 1
                start_val_inner = clock();
                #endif // TIMING

            ucon_c_func(i_x, d_now, &n_vec[0], &c_vec[0], &nb_vec[0],
                        &nb_plus[0], &c_plus[0],
                        &i_nb_plus_f[0], &Nnb_plus_db[0],
                        &i_db_plus_exp[0], &nb_plus_exp[0], &c_plus_exp[0],
                        &psi[0], &xit[0], &w[i_w_scale_func(i_u, i_x, 0, 0, 0, Nu, Nx, Npsi, Nxit)],
                        beta, rho, Gamma, ra, rd,
                        Nx, Ndb, Nnb_max, Npsi, Nxit, Nn,
                        err_flag, &t_vec[0]);

                            // note: w is state dependent (on both u and x)
                            // shape(w) = ((Nu,Nx),(Nx,Npsi,Nxit))

                #if TIMING == 1
                end_val_inner = clock();
                t_vec[1] += (end_val_inner - start_val_inner)*1000/CLOCKS_PER_SEC;
                #endif // TIMING

                if(*err_flag == 1){ printf("\nEXIT_FAILURE (ucon_c_func)\n\n"); exit(EXIT_FAILURE);}

                interp(Nnb_db[i_u*Nx*Ndb + i_x*Ndb + Ndb-1], Nn,
                       &nb[i_u*Nx*Ndb*Nnb_max + i_x*Ndb*Nnb_max + (Ndb-1)*Nnb_max + 0],
                       &ucon_c_d[i_u*Nx*Nd*Nnb_max + i_x*Nd*Nnb_max + i_d*Nnb_max + 0],
                       &nb_vec[0], &c_vec[0], err_flag);

                    // shape(Nnb_db)   = (Nu, Nx, Ndb)
                    // shape(nb)       = (Nu, Nx, Ndb, Nnb_max)
                    // shape(ucon_c_d) = (Nu, Nx, Nd, Nnb_max)

        } // x
        } // u

    } // d (implicit barrier)

        #if TIMING == 1
        end_val = clock();
        t_vec[0] += (end_val - start_val)*1000/CLOCKS_PER_SEC;
        #endif // TIMING


    ///////////////////////////
    // 4. grid search over d //
    ///////////////////////////

        #if TIMING == 1
        start_val = clock();
        #endif // TIMING

    #pragma omp for schedule(dynamic) collapse(3) // dynamic is chosen due to unbalanced work load
    for(i_u  = 0; i_u  < Nu;      i_u++){
    for(i_x  = 0; i_x  < Nx;      i_x++){
    for(i_nb = 0; i_nb < Nnb_max; i_nb++){

        // continue if never in the state space (bad for work balancing)
        if(i_nb >= Nnb_db[i_u*Nx*Ndb + i_x*Ndb + (Ndb-1)]){ continue;}

            // shape(Nnb_db) = (Nu,Nx,Ndb)

        /////////////////////////////////
        // set up borrowing constraint //
        /////////////////////////////////

        if(i_x == 1){
            varphi_now = 0.0;
        } else if(i_u == 1){
            varphi_now = varphi*chi;
        } else {
            varphi_now = varphi;
        }


        ////////////////
        // minimum db //
        ////////////////

        i_db_lower = 0;

        // check if there is any border nodes
        if(i_nb <= i_nb_f[i_u*Nx*Ndb + i_x*Ndb + i_db_lower]){

                // shape(i_nb_f) = (Nu,Nx,Ndb)

            while(i_nb <= i_nb_f[i_u*Nx*Ndb + i_x*Ndb + i_db_lower] && i_db_lower < Ndb){

                // c = 0 at state space border
                if(i_nb == i_nb_f[i_u*Nx*Ndb + i_x*Ndb + i_db_lower]){

                    index = index_func(i_u, i_x, i_db_lower, i_nb, Nu, Nx, Ndb, Nnb_max);

                    d[index]       = 0;
                    Delta[index]   = 0;
                    c[index]       = 0;
                    n[index]       = 0;
                    a[index]       = 0;
                    vt[index]      = 0;

                    // shape(d) = ... = (Nu, Nx, Ndb, Nnb_max)

                }

                i_db_lower++;

            }

            // continue if all db-nodes have been visited
            if(i_db_lower == Ndb){ continue; }

        }


        ///////////////////////////////////////
        // grid search A: top->down approach //
        ///////////////////////////////////////

            #if TIMING == 1
            start_val_inner = clock();
            #endif // TIMING

        // a. find d-choice at highest db
        i_db = Ndb-1;
        grid_search(i_u, i_x, i_db, i_nb, 0, Nd-1, &i_d_max[0],
                    &v[0], &d[0], &Delta[0], &c[0], &n[0], &a[0], &vt[0], &ucon_c_d[0],
                    &db[0], &nb[0], &db_plus[0], &nb_plus[0], &vt_plus[0],
                    &i_nb_plus_f[0], &Nnb_plus_db[0],
                    &i_db_plus_exp[0], &nb_plus_exp[0], &vt_plus_exp[0],
                    &psi[0], &xit[0], &w_scale[i_w_scale_func(i_u, i_x, 0, 0, 0, Nu, Nx, Npsi, Nxit)],
                    beta, rho, Gamma,
                    ra, rd, lambda, eta, varphi_now,
                    Nu, Nx, Ndb, Nnb_max, Npsi, Nxit, Nd, d_stepsize,
                    err_flag, t_vec);

        index_db_max = index_func(i_u, i_x, i_db, i_nb, Nu, Nx, Ndb, Nnb_max);

        // b. copy to all db-nodes where the same choice is available (i.e only rejected choices removed)
        i_db_upper = Ndb-1;
        while(i_db_upper-1 >= i_db_lower && db[i_db_upper-1] >= d[index_db_max]){

            i_db_upper--;
            index = index_func(i_u, i_x, i_db_upper, i_nb, Nu, Nx, Ndb, Nnb_max);

            v[index]     = v[index_db_max];
            d[index]     = d[index_db_max];
            Delta[index] = Delta[index_db_max];
            c[index]     = c[index_db_max];
            n[index]     = n[index_db_max];
            a[index]     = a[index_db_max];
            vt[index]    = vt[index_db_max];

        }

            #if TIMING == 1
            end_val_inner = clock();
            t_vec[7] += (end_val_inner - start_val_inner)*1000/CLOCKS_PER_SEC;
            #endif // TIMING

        // c. continue if all db-nodes have been visited
        if(i_db_upper == i_db_lower){continue;}


        /////////////////////////////////////
        // grid search B1: bottom, initial //
        /////////////////////////////////////

            #if TIMING == 1
            start_val_inner = clock();
            #endif // TIMING

        // a. find d-choice at LOWEST db
        i_db = i_db_lower;
        grid_search(i_u, i_x, i_db, i_nb, 0, Nd-1, &i_d_max[0],
                    &v[0], &d[0], &Delta[0], &c[0], &n[0], &a[0], &vt[0], &ucon_c_d[0],
                    &db[0], &nb[0], &db_plus[0], &nb_plus[0], &vt_plus[0],
                    &i_nb_plus_f[0], &Nnb_plus_db[0],
                    &i_db_plus_exp[0], &nb_plus_exp[0], &vt_plus_exp[0],
                    &psi[0], &xit[0], &w_scale[i_w_scale_func(i_u, i_x, 0, 0, 0, Nu, Nx, Npsi, Nxit)],
                    beta, rho, Gamma,
                    ra, rd, lambda, eta, varphi_now,
                    Nu, Nx, Ndb, Nnb_max, Npsi, Nxit, Nd, d_stepsize,
                    err_flag, t_vec);

            i_d_db_min = i_d_max[0];
            index_db_min = index_func(i_u, i_x, i_db, i_nb, Nu, Nx, Ndb, Nnb_max);

            #if TIMING == 1
            end_val_inner = clock();
            t_vec[8] += (end_val_inner - start_val_inner)*1000/CLOCKS_PER_SEC;
            #endif // TIMING

        // b. continue if all db-nodes have been visited
        i_db_lower++;
        if(i_db_upper == i_db_lower){continue;}


        /////////////////////////////
        // grid search B2: mid->up //
        /////////////////////////////

            #if TIMING == 1
            start_val_inner = clock();
            #endif // TIMING

        // a. find i_db_mid
        //  req. 1: db[i_db_mid-1] > d[index_db_min]
        //  req. 2: db[i_db_mid] >= d[index_db_min]

        i_db_mid = i_db_lower;
        while(db[i_db_mid] < d[index_db_min]){i_db_mid++;}

        // b. upward search (bounded from below)
        i_d_lower = 0;
        for(i_db = i_db_mid; i_db < i_db_upper; i_db++){

            grid_search(i_u, i_x, i_db, i_nb, i_d_lower, Nd-1, &i_d_max[0],
                        &v[0], &d[0], &Delta[0], &c[0], &n[0], &a[0], &vt[0], &ucon_c_d[0],
                        &db[0], &nb[0], &db_plus[0], &nb_plus[0], &vt_plus[0],
                        &i_nb_plus_f[0], &Nnb_plus_db[0],
                        &i_db_plus_exp[0], &nb_plus_exp[0], &vt_plus_exp[0],
                        &psi[0], &xit[0], &w_scale[i_w_scale_func(i_u, i_x, 0, 0, 0, Nu, Nx, Npsi, Nxit)],
                        beta, rho, Gamma,
                        ra, rd, lambda, eta, varphi_now,
                        Nu, Nx, Ndb, Nnb_max, Npsi, Nxit, Nd, d_stepsize,
                        err_flag, t_vec);

            i_d_lower = i_d_max[0];

        } // db

            #if TIMING == 1
            end_val_inner = clock();
            t_vec[9] += (end_val_inner-start_val_inner)*1000/CLOCKS_PER_SEC;
            #endif // TIMING

        // c. continue if all db-nodes have been visited
        if(i_db_mid == i_db_lower){continue;}


        /////////////////////////////////
        // grid search B3: mid->bottom //
        /////////////////////////////////

            #if TIMING == 1
            start_val_inner = clock();
            #endif // TIMING

        // downward search (fixed upper bound with exit)
        for(i_db = i_db_mid-1; i_db >= i_db_lower; i_db--){

            grid_search(i_u, i_x, i_db, i_nb, 0, i_d_db_min, &i_d_max[0],
                        &v[0], &d[0], &Delta[0], &c[0], &n[0], &a[0], &vt[0], &ucon_c_d[0],
                        &db[0], &nb[0], &db_plus[0], &nb_plus[0], &vt_plus[0],
                        &i_nb_plus_f[0], &Nnb_plus_db[0],
                        &i_db_plus_exp[0], &nb_plus_exp[0], &vt_plus_exp[0],
                        &psi[0], &xit[0], &w_scale[i_w_scale_func(i_u, i_x, 0, 0, 0, Nu, Nx, Npsi, Nxit)],
                        beta, rho, Gamma,
                        ra, rd, lambda, eta, varphi_now,
                        Nu, Nx, Ndb, Nnb_max, Npsi, Nxit, Nd, d_stepsize,
                        err_flag, t_vec);

            // if same choice as at db_min
            if(i_d_max[0] == i_d_db_min){

                // then copy to all smaller db-nodes
                while(i_db >= i_db_lower){

                    i_db--;

                    index = index_func(i_u, i_x, i_db, i_nb, Nu, Nx, Ndb, Nnb_max);

                    v[index]     = v[index_db_min];
                    d[index]     = d[index_db_min];
                    Delta[index] = Delta[index_db_min];
                    c[index]     = c[index_db_min];
                    n[index]     = n[index_db_min];
                    a[index]     = a[index_db_min];
                    vt[index]    = vt[index_db_min];

                }
                break;
            }

                #if TIMING == 1
                end_val_inner = clock();
                t_vec[10] += (end_val_inner-start_val_inner)*1000/CLOCKS_PER_SEC;
                #endif // TIMING

        } // db


    } // nb
    } // x
    } // u (implicit barrier)

        #if TIMING == 1
        end_val = clock();
        t_vec[4] += (end_val-start_val)*1000/CLOCKS_PER_SEC;
        #endif // TIMING

        #if CHECKSUM == 1
        #pragma omp master
        {
        chekksum = 0;
        max_diff = 0;

            for(i_u = 0; i_u < Nu; i_u++){
            for(i_x = 0; i_x < Nx; i_x++){
            for(i_db = 0; i_db < Ndb; i_db++){
            for(i_nb = 0; i_nb < Nnb_db[i_db]; i_nb++){

                index = i_u*Nx*Ndb*Nnb_max
                        + i_x*Ndb*Nnb_max
                        + i_db*Nnb_max
                        + i_nb_f[i_u*Nx*Ndb + i_x*Ndb + i_db] + i_nb;
                chekksum += vt[index];

                if(t <= tcon && vt[index] > 0){
                    max_diff = MAX(max_diff, fabs(v[index]-v_plus[index]));
                }

            } } } }

        log_txt = fopen("log.txt", "a");
        if(t > tcon){
            fprintf(log_txt, "t = %d (checksum = %5.8f)\n", t, chekksum);
        } else {
            fprintf(log_txt, "t = %d (checksum = %5.8f) (max_diff = %5.8f)\n", t, chekksum, max_diff);
        }
        fclose(log_txt);
        }
        #endif // CHECKSUM

        #pragma omp barrier

    } // t

        #if VTICHECK == 1
        #pragma omp for schedule(dynamic) collapse(3) // dynamic is chosen due to unbalanced work load
        for(i_u  = 0; i_u  < Nu;      i_u++){
        for(i_x  = 0; i_x  < Nx;      i_x++){
        for(i_nb = 0; i_nb < Nnb_max; i_nb++){

            // continue if never in the state space (bad for work balancing)
            if(i_nb >= Nnb_db[i_u*Nx*Ndb + i_x*Ndb + (Ndb-1)]){ continue;}

                // shape(Nnb_db) = (Nu,Nx,Ndb)

            /////////////////////////////////
            // set up borrowing constraint //
            /////////////////////////////////

            if(i_x == 1){
                varphi_now = 0.0;
            } else if(i_u == 1){
                varphi_now = varphi*chi;
            } else {
                varphi_now = varphi;
            }

            ////////////////
            // minimum db //
            ////////////////

            i_db_lower = 0;

            // check if there is any border nodes
            if(i_nb <= i_nb_f[i_u*Nx*Ndb + i_x*Ndb + i_db_lower]){

                    // shape(i_nb_f) = (Nu,Nx,Ndb)

                while(i_nb <= i_nb_f[i_u*Nx*Ndb + i_x*Ndb + i_db_lower] && i_db_lower < Ndb){
                    i_db_lower++;
                }

                // continue if all db-nodes have been visited
                if(i_db_lower == Ndb){ continue; }

            }

        for(i_db = i_db_lower; i_db < Ndb; i_db++){

            vfi_check_now = grid_search_full(i_u, i_x, i_db, i_nb, &i_d_max[0],
                              &v[0], &db[0], &nb[0], &db_plus[0], &nb_plus[0], &vt_plus[0],
                              &i_nb_plus_f[0], &Nnb_plus_db[0],
                              &i_db_plus_exp[0], &nb_plus_exp[0], &vt_plus_exp[0],
                              &psi[0], &xit[0], &w_scale[i_w_scale_func(i_u, i_x, 0, 0, 0, Nu, Nx, Npsi, Nxit)],
                              beta, rho, Gamma,
                              ra, rd, lambda, eta, varphi_now,
                              Nu, Nx, Ndb, Nnb_max, Npsi, Nxit, Nd, d_stepsize,
                              err_flag, t_vec);

            #pragma omp critical
            vfi_check[0] = MAX(vfi_check_now, vfi_check[0]);

        } // db

        } // nb
        } // x
        } // u (implicit barrier)


        log_txt = fopen("log.txt", "a");
        fprintf(log_txt, "\n value function check: %5.8f \n", vfi_check[0]);
        fclose(log_txt);
        #endif // VTICHECK

        // free memory
        free(n_vec);
        free(c_vec);
        free(nb_vec);
        free(i_d_max);
        free(i_db_plus_exp);
        free(nb_plus_exp);
        free(c_plus_exp);
        free(vt_plus_exp);

	}

    //////////////////
    // end parallel //
    //////////////////

        #if TIMING_TOT == 1 && TIMING == 0
        end = clock();
        t_tot = (end-start)*1000/CLOCKS_PER_SEC;

        log_txt = fopen("log.txt", "a");
        fprintf(log_txt, "\n total: %d mins \n",
                         t_tot/1000/60);
        fclose(log_txt);
        #endif

        #if TIMING == 1
        end = clock();
        t_tot = (end-start)*1000/CLOCKS_PER_SEC;
        log_txt = fopen("log.txt", "a");

        fprintf(log_txt, "\n total: %d mins \n\n"
                "  part 1: ucon_c:      %d mins \n"
                "     ucon_c_func:       %d mins \n"
                "     interp:            %d mins \n"
                "     sum (pow):         %d mins \n"
                "  part 2: grid search: %d mins \n"
                "     eval:              %d mins \n"
                "     interp:            %d mins \n",
                t_tot/1000/60,
                t_vec[0]/1000/60,
                t_vec[1]/1000/60,
                t_vec[2]/1000/60,
                t_vec[3]/1000/60,
                t_vec[4]/1000/60,
                t_vec[5]/1000/60,
                t_vec[6]/1000/60);

            fprintf(log_txt, "\n choice-bounding: \n"
                    "  A:  %d mins \n"
                    "  B1: %d mins \n"
                    "  B2: %d mins \n"
                    "  B3: %d mins \n",
                    t_vec[7]/1000/60,
                    t_vec[8]/1000/60,
                    t_vec[9]/1000/60,
                    t_vec[10]/1000/60);

        free(t_vec);
        #endif

    /////////////////
    // free memory //
    /////////////////

    free(err_flag);

    return;

} // v_func
