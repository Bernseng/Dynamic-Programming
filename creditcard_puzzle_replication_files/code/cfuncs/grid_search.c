// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

// perform a grid-search of a range of d-values and saves the optimal choice and implications
void grid_search(int i_u, int i_x, int i_db, int i_nb,
                 long int i_d_lower, long int i_d_upper,  long int *i_d_max,
                 double *v, double *d, double *Delta, double *c, double *n, double *a, double *vt, double *ucon_c_d,
                 double *db, double *nb,  double *db_plus, double *nb_plus, double *vt_plus,
                 int *i_nb_plus_f, int *Nnb_plus_db,
                 int *i_db_plus_exp, double *nb_plus_exp, double *vt_plus_exp,
                 double *psi, double *xit, double *w_scale,
                 double beta, double rho, double Gamma,
                 double ra, double rd, double lambda, double eta, double varphi,
                 int Nu, int Nx, int Ndb, int Nnb_max, int Npsi, int  Nxit, long int Nd, double d_stepsize,
                 int *err_flag, long int *t_vec)
{

    long int i_d, index, indexc;
    double db_now, nb_now, d_now, c_d, n_d, value;
    double v_max, c_argmax, c_target;
    clock_t start_val, end_val;

    index = index_func(i_u, i_x, i_db, i_nb, Nu, Nx, Ndb, Nnb_max);

    v_max = -INFINITY;
    nb_now = nb[index];
    db_now = db[i_db];

    for(i_d = i_d_lower; i_d < i_d_upper+1; i_d++){

        // a. d-choice
        d_now = (double)i_d*d_stepsize;

            // simple checks
            if(d_now > MAX(db_now, eta*nb_now+varphi)){ break;};  // illegal, even worse for higher d
            if(d_now < -nb_now){ continue;};                      // forces c = 0, better for higher d

        // b. c-choice (implies a n-choice)
        c_d = con_c_func(i_u, i_x, i_nb, i_d, &ucon_c_d[0],
                         d_now, db_now, nb_now,
                         eta, varphi,
                         Nu, Nx, Nd, Nnb_max, err_flag);

            if(c_d <= 0){continue;}

        n_d = nb_now - c_d;

        // c. value of choice (given db, nb, d and c_d)

            #if TIMING == 1
            start_val = clock();
            #endif // TIMING

        value = value_func(i_x, d_now, n_d, c_d,
                           &db_plus[0], &nb_plus[0], &i_nb_plus_f[0], &Nnb_plus_db[0], &vt_plus[0],
                           &i_db_plus_exp[0], &nb_plus_exp[0], &vt_plus_exp[0],
                           &psi[0], &xit[0], &w_scale[0],
                           beta, rho, Gamma, lambda, ra, rd,
                           Nx, Ndb, Nnb_max, Npsi, Nxit, err_flag, &t_vec[0]);

                           // note: w is state dependent (both on s and u)
                           // shape(w_scale) = (Nx,Npsi,Nxit)

            #if TIMING == 1
            end_val = clock();
            t_vec[5] += (end_val - start_val)*1000/CLOCKS_PER_SEC;
            #endif // TIMING

        // d. max-iteration
        if(value > v_max){
            v_max = value;
            i_d_max[0] = i_d;
            c_argmax = c_d;
        }

    } // d

        v[index] = v_max;

        d[index]       = i_d_max[0]*d_stepsize;
        Delta[index]   = d[index]-db_now;
        c[index]       = c_argmax;
        n[index]       = nb_now - c_argmax;
        a[index]       = nb_now + d[index] - c_argmax;
        vt[index]      = -1.0/v_max;

        indexc = i_u*Nx*Nd*Nnb_max
                + i_x*Nd*Nnb_max
                + (i_d-1)*Nnb_max
                + i_nb;

        c_target = ucon_c_d[indexc];

        if(vt[index] <= 0){
           printf("\nEXIT_FAILURE, (v_max not positive (%g %g) (%d %d %d %d) (%g %g %g %g) (%g %g %g)))\n\n",
                  v_max, vt[index], i_u, i_x, i_db, i_nb, db_now, nb_now, nb[index-1], nb[index+1], d_now, c_d, c_target);
           exit(EXIT_FAILURE);

        }

    return;
}
