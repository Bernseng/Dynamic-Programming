// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

// CONTAINS:
// 1. nb_plus_func
// 2. i_db_plus_exp_func
// 3. n_vec_func
// 4. con_c_func
// 5. index_func


////////////////
// 1. NB PLUS //
////////////////

// fills nb_plus_exp given realizations of psi and xit.
void nb_plus_func(double d, double n, double *nb_plus_exp,
                  double *psi, double *xit, int Npsi, int Nxit,
                  double Gamma, double ra, double rd, int *err_flag){

    int i_xit, i_psi, counter;

    counter = 0;
    for(i_psi = 0; i_psi < Npsi; i_psi++){ // permanent shock
    for(i_xit = 0; i_xit < Nxit; i_xit++){ // transitory shock

        nb_plus_exp[counter] = 1.0/(Gamma*psi[i_psi])*( (1.0+ra)*n - (rd-ra)*d ) + xit[i_xit];
        counter++;

    }
    }

        // OUTPUT CHECKS
        #if(CHECKS == 1)
        if(*err_flag == 1) return;

        // monotonicity: nb_plus_exp (asc)
        counter = 0;
        for(i_psi = 0; i_psi < Npsi; i_psi++){
        for(i_xit = 0; i_xit < Nxit; i_xit++){

            if(i_xit > 0 && nb_plus_exp[counter] < nb_plus_exp[counter-1]){
                *err_flag = 1;
                printf("nb_plus_func: nb_plus_exp[%d] = %g < %g = nb_plus_exp[%d]\n",
                        counter, nb_plus_exp[counter], nb_plus_exp[counter-1], counter-1);
            }
            if(*err_flag == 1) return;
            counter++;

        }
        }
        #endif // CHECKS

    return;

} // nb_plus_func

double nb_plus_func_min(double d, double n, double *nb_plus_exp,
                  double *psi, double *xit, int Npsi, int Nxit,
                  double Gamma, double ra, double rd, int *err_flag){

    int i_xit, i_psi, counter;
    double nb_min = HUGE_VAL;

    counter = 0;
    for(i_psi = 0; i_psi < Npsi; i_psi++){ // permanent shock
    for(i_xit = 0; i_xit < Nxit; i_xit++){ // transitory shock

        nb_plus_exp[counter] = 1.0/(Gamma*psi[i_psi])*( (1.0+ra)*n - (rd-ra)*d ) + xit[i_xit];
        nb_min = MIN(nb_min, nb_plus_exp[counter]);

    }
    }

    return nb_min;

}

//////////////////////////////
// 2. DB PLUS EXP (INDEXES) //
//////////////////////////////

// fills i_db_plus_exp.
void i_db_plus_exp_func(double d, double *db_plus, int *i_db_plus_exp,
                        double *psi, double Gamma, double lambda,
                        int Ndb, int Npsi, int *err_flag)
{

    int i_psi, imin=0, i_db_plus;
    double db_plus_now, diff_left, diff_right;

    for(i_psi = 0; i_psi < Npsi; i_psi++){

        db_plus_now = (1.0 - lambda) * d / (Gamma*psi[i_psi]);

        if(db_plus_now >= db_plus[Ndb-1]){

            i_db_plus_exp[i_psi] = Ndb-1;

        } else {

            // i. index below
            imin = binary_search(imin, Ndb, &db_plus[0], &db_plus_now, err_flag);

                #if CHECKS == 1
                if(*err_flag == 1) return;
                #endif // CHECKS

            // ii. nearest neighbor
            diff_left  = db_plus_now - db_plus[imin];
            diff_right = db_plus[imin+1] - db_plus_now;

            if(diff_left < diff_right){
                i_db_plus_exp[i_psi] = imin;
            } else {
                i_db_plus_exp[i_psi] = imin+1;
            }


        }

    } // psi

        // OUTPUT CHECKS
        #if CHECKS == 1
        // monotonicity: i_db_plus_exp (asc)
        for(i_psi = 1; i_psi < Npsi; i_psi++){
            if(i_db_plus_exp[i_psi] < i_db_plus_exp[i_psi-1]){
                *err_flag = 1;
                printf("i_db_plus_exp_func: i_db_plus[%d] < i_db_plus[%d]\n",i_psi, i_psi-1);
            }
            if(*err_flag == 1) return;
        }
        #endif // CHECKS

    return;
}


//////////////
// 3. N VEC //
//////////////

// fills n_vec.
void n_vec_func(double d, double *n_vec, double *kappa_plus, int *i_db_plus_exp,
                double *psi, double *xit, double Gamma, double ra, double rd,
                double n_max, double phi_n, double eps,
                int Nu, int Nx, int Ndb, int Npsi, int Nn, int *err_flag)

{

    int i_u, i_x, i_psi;
    double n_min=-INFINITY, n_min_now;

    // a. find n_min (unemployed)
    for(i_u = 0; i_u < Nu; i_u++){
    for(i_x = 0; i_x < Nx; i_x++){
    for(i_psi = 0; i_psi < Npsi; i_psi++){

        if(i_u == 0){
            n_min_now = (Gamma*psi[i_psi]*(kappa_plus[i_u*Nx*Ndb + i_x*Ndb + i_db_plus_exp[i_psi]]-xit[1]) + (rd-ra)*d) / (1+ra);
        } else {
            n_min_now = (Gamma*psi[i_psi]*(kappa_plus[i_u*Nx*Ndb + i_x*Ndb + i_db_plus_exp[i_psi]]-xit[0]) + (rd-ra)*d) / (1+ra);
        }
        n_min = MAX(n_min, n_min_now);

            // shape(kappa_plus) = (Nu,Nx,Ndb)

    } // psi
    } // x
    } // u

    // a. fill n_vec
    n_vec[0] = n_min;
    nonlinespace(&n_vec[1], n_min+eps, n_max, Nn-1, phi_n, err_flag);

        // OUTPUT CHECKS
        #if(CHECKS == 1)
        if(*err_flag == 1) return;
        #endif // CHECKS

    return;
}


//////////////////////
// 4. CONSTRAINED C //
//////////////////////

double con_c_func(int i_u, int i_x, int i_nb, long int i_d, double *ucon_c_d,
                  double d, double db, double nb,
                  double eta, double varphi,
                  int Nu, int Nx, long int Nd, int Nnb_max, int *err_flag)
{
    long int index;
    double c_max_d, c_target;

    // a. assets cannot be non-negative
    if(d <= db){

        c_max_d = nb + d;

    // b. satisfying the borrowing contract
    } else {

		if(eta > 0){
			c_max_d = nb + MIN(d, 1.0/eta*(varphi-d));
		} else {
			c_max_d = nb + d;
		}

    }

    // c. constrained consumption
    index = i_u*Nx*Nd*Nnb_max
            + i_x*Nd*Nnb_max
            + i_d*Nnb_max
            + i_nb;

    c_target = ucon_c_d[index];

        // shape(c_uncon_d) = (Nu,Nx,Nd,Nnb_max)

    return MIN(c_target, c_max_d);

}


//////////////
// 5. INDEX //
//////////////

long int index_func(int i_u, int i_x, int i_db, int i_nb,
                    int Nu,  int Nx,  int Ndb,  int Nnb_max)
{
    long int index;

    index = i_u*Nx*Ndb*Nnb_max +
            + i_x*Ndb*Nnb_max
            + i_db*Nnb_max
            + i_nb;

        // shape() = (Nu,Nx,Ndb,Nnb_max)

    return index;
}

long int i_w_scale_func(int i_u, int i_x, int i_x_plus, int i_psi, int i_xit,
                        int Nu, int Nx, int Npsi, int Nxit)
{
    long int index;

    index = i_u         *Nxit*Npsi*Nx*Nx
            + i_x       *Nxit*Npsi*Nx
            + i_x_plus  *Nxit*Npsi
            + i_psi     *Nxit
            + i_xit;

    return index;

}
