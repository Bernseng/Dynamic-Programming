// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

// perform a grid-search of a range of d-values and saves the optimal choice and implications
double grid_search_full(int i_u, int i_x, int i_db, int i_nb, long int *i_d_max,
                        double *v, double *db, double *nb,  double *db_plus, double *nb_plus, double *vt_plus,
                        int *i_nb_plus_f, int *Nnb_plus_db,
                        int *i_db_plus_exp, double *nb_plus_exp, double *vt_plus_exp,
                        double *psi, double *xit, double *w_scale,
                        double beta, double rho, double Gamma,
                        double ra, double rd, double lambda, double eta, double varphi,
                        int Nu, int Nx, int Ndb, int Nnb_max, int Npsi, int  Nxit, long int Nd, double d_stepsize,
                        int *err_flag, long int *t_vec)
{

    long int i_d, i_c,index;
    double db_now, nb_now, d_now, c_max_d, c_d, n_d, value;
    double v_max, c_argmax;

    index = index_func(i_u, i_x, i_db, i_nb, Nu, Nx, Ndb, Nnb_max);

    v_max = -INFINITY;
    nb_now = nb[index];
    db_now = db[i_db];

    for(i_d = 0; i_d < Nd; i_d++){

        // a. d-choice
        d_now = (double)i_d*d_stepsize;

        // b. find c_max
        if(d_now <= db_now){
            c_max_d = nb_now + d_now;

        } else {

            if(eta > 0){
                c_max_d = nb_now + MIN(d_now, 1.0/eta*(varphi-d_now));
            } else {
                c_max_d = nb_now + d_now;
            }

        }

            // simple checks
            if(d_now > MAX(db_now, eta*nb_now+varphi)){ break;};  // illegal, even worse for higher d
            if(c_max_d <= 0){ continue;};                         // forces c = 0, might be better for higher d

    for(i_c = 0; i_c < Nd; i_c++){

        // c. c and n choice
        c_d = ((double)i_c / (double)(Nd-1)) * c_max_d;
        n_d = nb_now - c_d;

        // d. value of choice (given db, nb, d and c_d)
        value = value_func(i_x, d_now, n_d, c_d,
                           &db_plus[0], &nb_plus[0], &i_nb_plus_f[0], &Nnb_plus_db[0], &vt_plus[0],
                           &i_db_plus_exp[0], &nb_plus_exp[0], &vt_plus_exp[0],
                           &psi[0], &xit[0], &w_scale[0],
                           beta, rho, Gamma, lambda, ra, rd,
                           Nx, Ndb, Nnb_max, Npsi, Nxit, err_flag, &t_vec[0]);

                           // note: w is state dependent (both on s and u)
                           // shape(w_scale) = (Nx,Npsi,Nxit)

        // e. max-iteration
        if(value > v_max){
            v_max      = value;
            i_d_max[0] = i_d;
            c_argmax   = c_d;
        }

    } // c
    } // d

    // f. return absolute difference
    return fabs(v[index] - v_max);


}
