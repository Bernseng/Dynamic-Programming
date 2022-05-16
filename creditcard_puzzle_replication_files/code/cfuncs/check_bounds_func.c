// version: 1.0.
// @author: Jeppe Druedahl, 2017.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "generic_funcs.c"
#include "base_funcs.c"

#define PRINTPROGRESS 0

#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#define MIN(X,Y) ((X)<(Y)?(X):(Y))

void main(){}

void check_bounds_func(int *check, int t, int i_u, int i_x, double db, double nb, double *db_plus,
                       double *kappa_plus,
                       double *psi, double *xit,
                       double Gamma, double ra, double rd, double lambda,
                       double eta, double varphi,
                       int Nu, int Nx, int Ndb, long int Nd, int Npsi, int Nxit, double d_stepsize)
{

    long int i_d;
    int i_psi, *i_db_plus_exp, index, i_u_plus, i_x_plus;
    double c, n, d_now, *nb_plus_exp, nb_plus_min;
    int *err_flag;
    FILE *log_txt;

    // a. allocate memory
    i_db_plus_exp   = (int*) calloc(Npsi, sizeof(int));
    nb_plus_exp     = (double*) calloc(Npsi, sizeof(double));
    err_flag        = (int*) calloc(1, sizeof(int));
    err_flag[0]     = 0;

    // b. calculate n
    c = 0;
    n = nb-c;

    for(i_d = 0; i_d < Nd; i_d++){

        d_now = (double)i_d*d_stepsize;
		nb_plus_min = HUGE_VAL;

        // c. check borrowing constraint
        if(d_now < -nb){
            continue;
        }

        if(d_now > db && d_now > eta*n  + varphi){ continue; }

        // d. find db plus indexes
        i_db_plus_exp_func(d_now, &db_plus[0], &i_db_plus_exp[0],
                           &psi[0], Gamma, lambda,
                           Ndb, Npsi, err_flag);

        // e. find nb plus (xit fixed, Nxit=1)
        //nb_plus_func(d_now, n, &nb_plus_exp[0], &psi[0], &xit[0], Npsi, 1, Gamma, ra, rd, err_flag);

        // f. run test
        check[0] = 1;

        for(i_psi    = 0; i_psi    < Npsi;  i_psi++){
        for(i_u_plus = 0; i_u_plus < Nu;    i_u_plus++){
        for(i_x_plus = 1; i_x_plus < Nx;    i_x_plus++){

            index = i_u_plus*Nx*Ndb + i_x_plus*Ndb + i_db_plus_exp[i_psi];

            if(i_u_plus == 0){
                nb_plus_min = nb_plus_func_min(d_now, n, &nb_plus_exp[0], &psi[i_psi], &xit[1], 1, Nxit-1, Gamma, ra, rd, err_flag);
            } else {
                nb_plus_min = nb_plus_func_min(d_now, n, &nb_plus_exp[0], &psi[i_psi], &xit[0], 1, 1, Gamma, ra, rd, err_flag);
            }

            //nb_plus_min = nb_plus_exp[i_psi];

			if(nb_plus_min < kappa_plus[index]){
                check[0] = 0;
                break;
            } // shape(kappa_plus) = (Nu,Nx,Ndb)

        } } }

        // g. if all survive the test: positive exit
        if(check[0] == 1){

            free(i_db_plus_exp);
            free(nb_plus_exp);
            free(err_flag);

			#if PRINTPROGRESS == 1
			log_txt = fopen("log_check_bounds.txt", "a");
			fprintf(log_txt, "positive (%d,%d,%d,%g,%g): d = %g, nb_plus_min = %g\n", t, i_u, i_x, db, nb, d_now, nb_plus_min);
			fclose(log_txt);
			#endif // PRINT

            return;
		}

    } // d

    // h. free memory
    free(i_db_plus_exp);
    free(nb_plus_exp);
    free(err_flag);

    // i. no d where all survive the test: negative exit
    check[0] = 0;

			#if PRINTPROGRESS == 1
			log_txt = fopen("log_check_bounds.txt", "a");
			fprintf(log_txt, "negative (%d,%d,%d,%g,%g): d = %g, nb_plus_min = %g\n", t, i_u, i_x, db, nb, d_now, nb_plus_min);
			fclose(log_txt);
			#endif // PRINT

	return;

}
