// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// number of threads in parallel region
#define MAXTHREADS MIN(24, omp_get_max_threads()-2)

#include "generic_funcs.c"

#define CHECKS 0
#define PRINTPROGRESS 0

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

void main(){}

void sim_func(int N, int T,
			  int *pop_age, int *pop_u, int *pop_x, double *pop_db, double *pop_nb,
			  double *pop_d, double *pop_c,
			  double *pop_n, double *pop_a,
			  double *pop_death, double *pop_u_val, double *pop_x_val,
			  double *pop_psi, double *pop_xi,
			  double *pop_P, double *pop_Y, double *pop_util,
			  int Nu, int Nx, int Ndb, double *db,
			  int Nnb_max, double *nb_T, int *i_nb_f_T, int *Nnb_db_T,
			  long int Nd, double d_stepsize,
			  double *d_T, double *ucon_c_d_T,
			  double rho, double Gamma, double mu,
			  double ra, double rd, double lambda, double eta, double varphi, double chi,
			  double pi_uw, double pi_uu, double *pi_x_lose, double *pi_x_gain, double pi_death)
{

    int *err_flag;

    // 1. memory allocation
    err_flag = (int*) calloc(1, sizeof(int));
    err_flag[0] = 0;

    // 2. loop: time
    //#pragma omp parallel num_threads(MAXTHREADS)
    //{

    FILE *log_txt;
    int t, i_u, i_x, i, i_plus, n, i_db_left, i_db_right, i_nb;
    int *i_nb_f, *Nnb_db;
    long i_d, i_d_below, i_d_above;
    double w_db, w_nb, below, above, left, right, diff_below, diff_above;
    double d_max, c_target, c_max, bc, varphi_now, checksum=0;
    double *nb, *d, *ucon_c_d;

    for(t = 0; t < T; t++){

	    #if PRINTPROGRESS == 1
        log_txt = fopen("log_sim.txt", "w");
        if(log_txt == NULL){exit(-1); }
        fprintf(log_txt, "log created\n");
        #endif // PRINT

    // 3. loop: households
    //#pragma omp for schedule(static)
    for(n = 0; n < N; n++){

        #if PRINTPROGRESS == 1
        fprintf(log_txt,"\n t = %d, n = %d\n",t,n);
        #endif // PRINT

        // a. current and future indexes
        i       = t    *N + n;
        i_plus  = (t+1)*N + n;

        // b. initialization
        if(t == 0 || pop_death[i] <= pi_death || pop_age[i] > 400 ){

            if(t == 0){
                pop_u[i] = 0;
                pop_x[i] = 0;
            }

            pop_P[i]  = pop_psi[i];
            pop_nb[i] = 1.0/12.0;

            if(pop_u[i] == 0){
                pop_nb[i] += pop_xi[i_plus];
                pop_Y[i]   = pop_xi[i_plus]*pop_P[i];
            } else {
                pop_nb[i] +=  mu;
                pop_Y[i]   = mu*pop_P[i];
            }

            pop_db[i]  = 0.0;
            pop_age[i] = 0;

        }

        i_u = pop_u[i];
        i_x = pop_x[i];

        if(i_x == 1){
            varphi_now = 0;
        } else if(i_u == 1){
            varphi_now = varphi*chi;
        } else {
            varphi_now = varphi;
        }

        #if PRINTPROGRESS == 1
        fprintf(log_txt,"  (u, x, db, nb) = (%d, %d, %g, %g)\n", i_u, i_x, pop_db[i], pop_nb[i]);
        #endif // PRINT

        // c. remove time and discrete dimensions (t = 0 due to inf horizon)
        nb = &nb_T[0*Nu*Nx*Ndb*Nnb_max + i_u*Nx*Ndb*Nnb_max + i_x*Ndb*Nnb_max];
        d  = &d_T[ 0*Nu*Nx*Ndb*Nnb_max + i_u*Nx*Ndb*Nnb_max + i_x*Ndb*Nnb_max];

            // shape(nb_T) = shape(d) = (T,Nu,Nx,Ndb,Nnb_max)

        i_nb_f = &i_nb_f_T[0*Nu*Nx*Ndb + i_u*Nx*Ndb + i_x*Ndb];
        Nnb_db = &Nnb_db_T[0*Nu*Nx*Ndb + i_u*Nx*Ndb + i_x*Ndb];

            // shape(i_nb_f_T) = shape(Nnb_db_T) = (T,Nu,Nx,Ndb)

        ucon_c_d = &ucon_c_d_T[0*Nu*Nx*Nd*Nnb_max + i_u*Nx*Nd*Nnb_max + i_x*Nd*Nnb_max];

            // shape(ucon_c_d_T) = (T,Nu,Nx,Nd,Nnb_max)

        // d. db-dimension (no extrapolation)

            #if CHECKS == 1
            if(pop_db[i] >= db[Ndb-1]){
                fprintf(log_txt,"pop_db[i] >= db[Ndb-1] \n");
                fclose(log_txt); exit(EXIT_FAILURE);
            }
            if(pop_db[i] < 0){
                fprintf(log_txt,"pop_db[i] >= db[Ndb-1] \n");
                fclose(log_txt); exit(EXIT_FAILURE);
            }
            #endif // PRINT

        i_db_left  = binary_search(0, Ndb, &db[0], &pop_db[i], err_flag);
        i_db_right = i_db_left + 1;

            // weights
            w_db = (pop_db[i]-db[i_db_left])/(db[i_db_right]-db[i_db_left]);

                #if CHECKS == 1
                if(w_db < 0){  *err_flag = 1; fprintf(log_txt,"w_db < 0\n");}
                if(w_db >= 1){ *err_flag = 1; fprintf(log_txt,"w_db >= 1\n");}
                if(*err_flag == 1) fclose(log_txt); exit(EXIT_FAILURE);
                #endif // CHECKS

        // e. nb-dimension (checks for extrapolation)

            #if CHECKS == 1
            if(pop_nb[i] <= nb[i_db_left*Nnb_max + i_nb_f[i_db_left]]){
                fprintf(log_txt,"pop_nb[i] too low at left \n");
                fclose(log_txt); exit(EXIT_FAILURE);
            }
            // no problem at left -> no problem at right
            #endif // PRINT

        if(pop_nb[i] >= nb[(Ndb-1)*Nnb_max + Nnb_db[Ndb-1]-2]){
            i_nb = Nnb_db[Ndb-1]-2;
        } else {
            i_nb = binary_search(0, Nnb_db[Ndb-1], &nb[(Ndb-1)*Nnb_max],
                                 &pop_nb[i], err_flag);
        }

            // weights
            w_nb = (pop_nb[i]-nb[(Ndb-1)*Nnb_max + i_nb])/(nb[(Ndb-1)*Nnb_max + i_nb+1]-nb[(Ndb-1)*Nnb_max + i_nb]);

                #if CHECKS == 1
                if(w_nb < 0){   *err_flag = 1; fprintf(log_txt,"w_nb < 0\n");}
                if(w_nb >= 1 && i_nb < Nnb_db[Ndb-1]-2){  *err_flag = 1; fprintf(log_txt,"w_nb = %g >= 1 (%g %g)\n", w_nb, nb[(Ndb-1)*Nnb_max + Nnb_db[Ndb-1]-2], pop_nb[i]);}
                if(*err_flag == 1) fclose(log_txt); exit(EXIT_FAILURE);
                #endif // CHECKS

            #if PRINTPROGRESS == 1
            fprintf(log_txt,"  nb in (%g %g)\n",nb[i_nb], nb[i_nb+1]);
            #endif // PRINT

        // f. choice: d

            below   = d[i_db_left*Nnb_max + i_nb    ];
            above   = d[i_db_left*Nnb_max + i_nb + 1];
            left    = below + w_nb * (above - below);

            below   = d[i_db_right*Nnb_max + i_nb    ];
            above   = d[i_db_right*Nnb_max + i_nb + 1];
            right   = below + w_nb * (above - below);

        pop_d[i] = left + w_db * (right - left);

		d_max = MAX(pop_db[i], eta*pop_nb[i] + varphi_now);
		pop_d[i] = MIN(d_max, pop_d[i]);

		#if PRINTPROGRESS == 1
        fprintf(log_txt,"  d = %g\n",pop_d[i]);
        #endif // PRINT

        // g. choice: c
        i_d_below = (int)floor(pop_d[i] / d_stepsize);
        i_d_above = i_d_below + 1;

        diff_below = pop_d[i] - i_d_below * d_stepsize;
        diff_above = i_d_above * d_stepsize - pop_d[i];

        if(diff_below < diff_above){
            i_d = i_d_below;
        } else {
            i_d = i_d_above;
        }

            below     = ucon_c_d[i_d*Nnb_max + i_nb    ];
            above     = ucon_c_d[i_d*Nnb_max + i_nb + 1];
            c_target  = below + w_nb * (above - below);

            #if PRINTPROGRESS == 1
            fprintf(log_txt,"  c_target = %g (%g %g %g)\n", c_target, below, above, w_nb);
            #endif // PRINT

            #if CHECKS == 1
            if(c_target <= 0){ fprintf(log_txt,"c_target =< 0, (%g, %g, %g, %d)\n",
                                       w_nb, below, above, i_d*Nnb_max + i_nb);
                               fclose(log_txt); exit(EXIT_FAILURE);}
            #endif // CHECKS

        if(pop_d[i] <= pop_db[i]){
            c_max = pop_nb[i] + pop_d[i];
        } else {
			if(eta > 0){
				c_max = pop_nb[i] + MIN(pop_d[i], 1.0/eta*(varphi_now-pop_d[i]));
			} else {
				c_max = pop_nb[i] + pop_d[i];
			}
		}

            #if PRINTPROGRESS == 1
            fprintf(log_txt,"  c_max = %g\n",c_max);
            #endif // PRINT

            #if CHECKS == 1
            if(c_max <= 0){ fprintf(log_txt,"c_max =< 0\n"); fclose(log_txt); exit(EXIT_FAILURE);}
            #endif // CHECKS

        pop_c[i] = MIN(c_target, c_max);

            #if PRINTPROGRESS == 1
            fprintf(log_txt,"  c = %g\n",pop_c[i]);
            #endif // PRINT

            #if CHECKS == 1
            if(pop_c[i] <= 0){ fprintf(log_txt,"pop_c[i] =< 0\n"); fclose(log_txt); exit(EXIT_FAILURE);}
            #endif // CHECKS

        // h. accumulation
        pop_util[i] = pow(pop_c[i], 1.0-rho) / (1.0-rho);
        pop_n[i]    = pop_nb[i] - pop_c[i];
        pop_a[i]    = pop_n[i]  + pop_d[i];

            // check borrowing constraint
            #if CHECKS == 1
            if(pop_d[i] > pop_db[i]){

                bc = pop_d[i] - (eta*pop_n[i] + varphi_now);
                if(bc > pow(10, -8)){
                    fprintf(log_txt,"d = %g, n = %g, nb = %g\n", pop_d[i], pop_n[i], pop_nb[i]);
                    fprintf(log_txt,"bc = %g > 0\n", bc);
                    exit(EXIT_FAILURE);
                }

            }
            #endif // CHECKS

        // g. next period
        pop_age[i_plus] = pop_age[i] + 1;
        pop_nb[i_plus]  =  1.0 / (Gamma*pop_psi[i_plus]) * ((1.0+ra)*pop_n[i]-(rd-ra)*pop_d[i]);
        pop_db[i_plus]  = (1.0-lambda) * pop_d[i] / (Gamma*pop_psi[i_plus]);
        pop_P[i_plus]   = pop_P[i]*pop_psi[i_plus];

            // u
            if(i_u == 0){

                if(pop_u_val[i_plus] <= pi_uw){
                    pop_u[i_plus] = 1;
                } else {
                    pop_u[i_plus] = 0;
                }

            } else {

                if(pop_u_val[i_plus] <= pi_uu){
                    pop_u[i_plus] = 1;
                } else {
                    pop_u[i_plus] = 0;
                }

            }

            // x
            if(pop_x[i] == 0){

                if(pop_x_val[i_plus] <= pi_x_lose[pop_u[i_plus]]){
                    pop_x[i_plus] = 1;
                } else {
                    pop_x[i_plus] = 0;
                }

            } else {

                if(pop_x_val[i_plus] <= pi_x_gain[pop_u[i_plus]]){
                    pop_x[i_plus] = 0;
                } else {
                    pop_x[i_plus] = 1;
                }

            }

            if(pop_d[i] <= 1e-6 && pop_x[i] == 0){
                pop_x[i_plus] = 0;
            }

            // nb, Y
            if(pop_u[i_plus] == 0){
                pop_nb[i_plus] += pop_xi[i_plus];
                pop_Y[i_plus]   = pop_xi[i_plus]*pop_P[i_plus];
            } else {
                pop_Y[i_plus]   = mu*pop_P[i_plus];
                pop_nb[i_plus] += mu;
            }

        #if PRINTPROGRESS == 1
            fprintf(log_txt,"  n    = %5.2f\n",pop_n[i]);
            fprintf(log_txt,"  a    = %5.2f\n",pop_a[i]);
            fprintf(log_txt,"  psi+ = %5.2f\n",pop_psi[i_plus]);
            fprintf(log_txt,"  xi+  = %5.2f\n",pop_xi[i_plus]);
            fprintf(log_txt,"  u+   = %d [%g %g %g]\n",pop_u[i_plus],pop_u_val[i_plus],pi_uw,pi_uu);
            fprintf(log_txt,"  x+   = %d [%g %g %g]\n",pop_x[i_plus],pop_x_val[i_plus],pi_xw,pi_xu);
            fprintf(log_txt,"  db+  = %5.2f\n",pop_db[i_plus]);
            fprintf(log_txt,"  nb+  = %5.2f\n",pop_nb[i_plus]);
        #endif // PRINT

		checksum += pop_n[i];

        } // N

        #if PRINTPROGRESS == 1
        fprintf(log_txt,"checksum = %g\n", checksum);
        fclose(log_txt);
        #endif // PRINT

    } // T
    //}

    free(err_flag);

} // sim_func
