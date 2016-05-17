#include <math.h>
#include <stdio.h>

void erste_Abl_func( double *erste_Abl, int *zufallszahlen, double theta, double mean1, double mean2, double mean3, int n1, int n2, int n3, int n){

  double erste_ableitung = 0;
  int i,j;

  for( i = 0; i < n; i++){
  	for( j = 0; j < zufallszahlen[i]; j++){
  		erste_ableitung = erste_ableitung + 1/(j + theta);
		}
	}

	erste_ableitung = erste_ableitung + n1 * log( theta / (mean1 + theta)) + n2 * log( theta / (mean2 + theta)) + n3 * log( theta / (mean3 + theta));

  erste_Abl[0] = erste_ableitung;

}




void newton_Shape( int *zufallszahlen, double *mean1_in, double *mean2_in, double *mean3_in, int *n1_in, int *n2_in, int *n3_in, double *theta_out){

  double mean1 = mean1_in[0];
  double mean2 = mean2_in[0];
  double mean3 = mean3_in[0];
  int n1 = n1_in[0];
  int n2 = n2_in[0];
  int n3 = n3_in[0];
  int n = n1 + n2 + n3;
  double erste_ableitung = 0;
  double erste_ableitung_help = 0;
  double erste_abl_func_in[1];
  double zweite_ableitung = 0;
  double theta = 0;
  double theta_min = pow(10, -2);
  double theta_max = pow(10, 6);
  //double theta_alt;
  double eps = pow(10,-6);
  int max_iter = 500;
  int counter_iter = 0;


  int i, j;
	// Startwertsuche: Multipliziere theta solange mit 10, bis die erste Ableitung groesser Null ist
  // Ueberpruefe, ob die erste Ableitung fuer theta_min bereits kleiner als Null ist, wenn ja, ist MLE fuer theta
  //  kleiner als theta_min und theta wird im Bereich kleiner als theta_min gesucht
  erste_abl_func_in[0] = 0;
  erste_Abl_func( erste_abl_func_in, zufallszahlen, theta_min, mean1, mean2, mean3, n1, n2, n3, n);
  erste_ableitung = erste_abl_func_in[0];
  theta = theta_min;
  if( erste_ableitung < 0 ){
    while( erste_ableitung < 0 ){
      theta = theta / 10;
      erste_abl_func_in[0] = 0;
      erste_Abl_func( erste_abl_func_in, zufallszahlen, theta, mean1, mean2, mean3, n1, n2, n3, n);
      erste_ableitung = erste_abl_func_in[0];
    }
	}
  else{
	  while( theta <= theta_max ){
      erste_ableitung_help = erste_ableitung;
      erste_abl_func_in[0] = 0;
      erste_Abl_func( erste_abl_func_in, zufallszahlen, theta, mean1, mean2, mean3, n1, n2, n3, n);
      erste_ableitung = erste_abl_func_in[0];

		  if(erste_ableitung < 0 && erste_ableitung_help > 0){
			  break;
		  }
		  else{
		    theta = theta * 5;
      }
    }
    theta = theta / 5;
	}

	if( theta*5 >  theta_max ){
		theta_out[0] = INFINITY;
	}
	else{
	  while( counter_iter < max_iter ){
	    //theta_alt = theta;
	    // Berechnen der ersten Ableitung
      erste_abl_func_in[0] = 0;
      erste_Abl_func( erste_abl_func_in, zufallszahlen, theta, mean1, mean2, mean3, n1, n2, n3, n);
      erste_ableitung = erste_abl_func_in[0];

      if( fabs( erste_ableitung ) < eps ){
        break;
      }
	    // Berechnen der zweiten Ableitung
	    zweite_ableitung = 0;
	    for( i = 0; i < n; i++){
	      for( j = 0; j < zufallszahlen[i]; j++){
	        zweite_ableitung = zweite_ableitung - 1/pow(j + theta,2);
	      }
	    }
	    zweite_ableitung = zweite_ableitung + n1 * mean1 / ( theta * (mean1 + theta)) + n2 * mean2 / ( theta * (mean2 + theta)) + n3 * mean3 / ( theta * (mean3 + theta));
      // Newton-Iteration
	    theta = theta - erste_ableitung / zweite_ableitung;
      counter_iter++;
	  }

    	theta_out[0] = (double) theta;
}


}
