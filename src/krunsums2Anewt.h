struct distribution_parameters{
  double exp;
  double var;
};

int *K_run_sums;


double a1, a2, a3, a31, a4, a41, a5, a6, a7, a8, a4_d1, p_i,
       a2_d, a3_d, a4_d, a5_d, a6_d;
 
void calculation_of_series( double i, double p) 

{
   
   p_i=pow(p, i);
  a2=(1-p*p_i)/(1-p)-1;
  a31=1-pow(p, 2)*p_i-(i+2)*(1-p)*p*p_i;
  a3=a31/pow(1-p, 2)-1;
  a41=2*(1-pow(p, 3)*p_i-(i+3)*pow(p, 2)*p_i*(1-p))-
      pow(1-p, 2)*(i+3)*(i+2)*p*p_i;
  a4=a41/pow(1-p, 3)-2;
  a5=a3-a2;
  a6=a4-3*a3+a2;
  a2_d=(1-p*p_i-(i+1)*p_i*(1-p))/pow(1-p, 2);
  a3_d=(2*a31-(i+2)*(i+1)*p_i*pow(1-p, 2))/pow(1-p, 3);
  a4_d1=(i+3)*(i+2)*(i+1)*p_i*pow(1-p, 3);
  a4_d=(3*a41-a4_d1)/pow(1-p, 4);
  a5_d=a3_d-a2_d;
  a6_d=a4_d-3*a3_d+a2_d;
  a8=(p*p-pow(p, 2)*p_i-i*p*p_i*(1-p))/pow(1-p, 2);
  a7=p*a6_d-a6+p*a5_d-a5+a8+a2;  
}

struct distribution_parameters distribution_k_run_sums_version4
					 (double p, int n, int k)
{

 struct distribution_parameters result;
  int j, b_site, e_site, site, loop;
 double i, m, r, b1, b2, q, a2_i, a5_i, a6_i, a7_i, a2_m, a5_m, a6_m,
        a7_m, a2_r, a5_r, a6_r, a7_r, TJ2, TJ1, TH2, TH1, TH2_false,
        TH1_false, tem, tem1, Mean, Var, Covar, Covar2, Covar3,
        parameter, number, judge, real_j, real_n;

 q=1-p;
 loop=k;
 judge=0;
 Mean=0;
 Var=0;
 Covar=0;
 Covar2=0;
 Covar3=0;

/* Calculate the Mean and the sum of variances of the  individual variables */


 for (j=k; j<=n-1; j++)
   {
     tem=0;
     real_j=j;
     real_n=n;
     tem=q*q*pow(p, real_j);
     if (tem<0)  tem=0;
     Mean=Mean+real_j*(real_n-real_j-1)*tem;
     Var=Var+real_j*real_j*(real_n-real_j-1)*tem*(1-tem);
     tem=tem/q;
     if (tem<=0)  break; 
     Mean=Mean+2*real_j*tem;
     Var=Var+2*real_j*real_j*tem*(1-tem);
   }


/* Calculate the covariances between the variables  denoting  patterns with 
   same size of heads. The parameter Covar denotes the  covariances between
   the mutually exclusive variables, and the Covar2 denotes the covariances
   between the reinforcing variables  */
 real_j=0;
 real_n=0;

 for (j=k; j<=n-1; j++)
   {
     real_j=j;
     real_n=n;

	     tem=0; 
	     b2=real_j+1;                  /* 1T vs. 2T & 1T */
	     if (b2>real_n-real_j+1) b2=real_n-real_j+1;
             e_site=b2;
	    
	     if (e_site==real_n-real_j+1)
	        tem=2*pow(real_j, 2)*q*pow(p, 2*real_j)*
		  ((real_n-real_j-1)*q*q+q);
             else
         	tem=2*pow(real_j, 2)*q*pow(p, 2*real_j)*real_j*q*q;  
	     if (tem<0) tem=0;
	     Covar=Covar+tem;

	     TJ2=2*real_j*pow(q, 2)*pow(p, real_j);
	     if (TJ2<0) TJ2=0;
	     TH2=real_j*pow(q, 2)*pow(p, real_j);
	     if (TH2<0) TH2=0;
	     TH1=real_j*q*pow(p, real_j);
	     if (TH1<0) TH1=0;
	     judge=pow(q, 4)*pow(p, 2*real_j);
	     if (judge<=0)  break;  

	                                    /* 2T vs. 2T & 1T */
	     a1=(real_n-2*real_j)-1;
	     if (a1>=1)
	       {
		 tem=a1*real_j*TJ2*TH2; 
		 if (tem<0) tem=0;
		 Covar=Covar+tem; 
		 tem=real_j*(real_j*TJ2*TH2+TJ2*TH1)-
		     (real_j*(real_j+1)/2)*TJ2*TH2;
		 if (tem<0) tem=0;
		 Covar=Covar+tem;
	       }
	     else
	       {
		 a2=(real_n-real_j)-1;
		 if (a2>=1)
		   {
		     tem=a2*((real_n-real_j-1)*TJ2*TH2+TJ2*TH1)
		         -(a2*(a2+1)/2)*TJ2*TH2;
		     if (tem<0) tem=0;
		     Covar=Covar+tem;
		   }
	       }

	     e_site=real_j+2;   /* 1T vs 2T & 1T (reinforcement) */
       	     tem1=0;
	     if (e_site==n-j+1)
	        tem1=2*pow(real_j, 2)*q*(1-q)*pow(p, 2*real_j);
	     else
	        if (e_site<n-j+1)
	          tem1=2*pow(real_j, 2)*q*q*(1-q)*pow(p, 2*real_j);
	     if (tem1<0) tem1=0;
	     Covar2=Covar2+tem1;

	                  /* 2T vs 2T & 1T (reinforcement) */
	     a3=(real_n-2*real_j)-2;
	     if (a3>=0)
	       {
		 tem1=(a3-1)*(1/q-1)*TJ2*TH2;
		 if (tem1<0) tem1=0;
		 Covar2=Covar2+tem1;
		 if (a3==0)
		   tem1=(1/q-1)*TJ2*TH1;
		 else
		   tem1=(1/q-1)*(TJ2*TH2+TJ2*TH1);
		 if (tem1<0) tem1=0;
		 Covar2=Covar2+tem1;
	       }
   }   


/* Calculate the covariances between the variables denoting patterns with
   unequal size of heads  */


 real_n=0;
 real_j=0;
 for (j=k; j<=n-1; j++)  
   {
     real_n=n;
     real_j=j;
     loop=loop+1;  
     tem=0; tem1=0;                      
     a2_i=0; a5_i=0; a6_i=0; a7_i=0;
     a2_m=0; a5_m=0; a6_m=0; a7_m=0;

     TJ2=2*real_j*pow(q, 2)*pow(p, real_j);
     if (TJ2<0) TJ2=0;
     TJ1=TJ2/q;
     if (TJ1<0) TJ1=0;
     TH2=TJ2/2;
     if (TH2<0) TH2=0;
     TH2_false=TH2/real_j;
     if (TH2_false<0) TH2_false=0;
     TH1=TH2/q;
     if (TH1<0) TH1=0;
     TH1_false=TH1/real_j;
     if (TH1_false<0) TH1_false=0;
     judge=pow(q, 4)*pow(p, 2*real_j);
     if (judge<=0) break; 
                                          /* 1T vs 2T & 1T */ 
     i=real_n-2*real_j-1;
     if (i>=1)
       {
	 calculation_of_series(i,p);
         a2_i=a2; a5_i=a5; a6_i=a6; a7_i=a7;
	 tem=(real_j*a2*TH2+real_j*a5*TH2_false+a2*TH1+a5*TH1_false)*TJ1;
	 if (tem<0) tem=0;
	 Covar=Covar+2*tem;
	                                 /* reinforcement */
	 tem1=((a2-p_i)*TH2+(a5-i*p_i)*TH2_false)*(1/q-1)*TJ1+
	      (1/q-1)*p_i*TJ1*TH1+(1/q-1)*i*p_i*TJ1*TH1_false;
	 if (tem1<0) tem1=0;
	 Covar2=Covar2+2*tem1;
       }
     m=n-j-1;
     calculation_of_series(m,p);
     a2_m=a2; a5_m=a5; a6_m=a6; a7_m=a7;
     a2=a2_m-a2_i;  a5=a5_m-a5_i; a6=a6_m-a6_i; a7=a7_m-a7_i;
     tem=(((real_n-real_j-1)*a2-a5)*TH2+((real_n-real_j-1)*a5-a6)*TH2_false+
          2*a2*TH1+2*a5*TH1_false)*TJ1;
     if (tem<0) tem=0;
     Covar=Covar+2*tem;


                                      /* 2T vs 2T & 1T */     
     tem=0; tem1=0;
     a2_i=0; a5_i=0; a6_i=0; a7_i=0;
     a2_m=0; a5_m=0; a6_m=0; a7_m=0;
     parameter=(real_n-real_j-2)/2.0;
     a1=modf(parameter, &number);
     if (number==(real_n-real_j-2)/2.0)  
       i=(real_n-3*real_j-2)/2.0;
     else
       i=(real_n-3*real_j-3)/2.0;

     if (i>=1)
       {
	 calculation_of_series(i,p);
         a2_i=a2; a5_i=a5; a6_i=a6; a7_i=a7;
	 tem=(((real_n-3*real_j-1)*(2*real_j+1))*a2+
	      (real_n-7*real_j-3)*a5-2*a6)*TJ2*TH2+
             (((real_n-3*real_j-1)*(2*real_j+1))*a5+
	      (real_n-7*real_j-3)*a6-2*a7)*TJ2*TH2_false;
	 if (tem<0) tem=0;
	 Covar=Covar+tem;
	 tem=(((3*real_j+1)*real_j*a2+(4*real_j+1)*a5+a6)/2)*TJ2*TH2+
	     (((3*real_j+1)*real_j*a5+(4*real_j+1)*a6+a7)/2)*TJ2*TH2_false+
             (real_j*a2+a5)*TJ2*TH1+(real_j*a5+a6)*TJ2*TH1_false;
	 if (tem<0) tem=0;
         Covar=Covar+2*tem;
	 
         /* calculate the covariances between the reinforced variables */
	 
	 tem1=((real_n-3*real_j-3)*(a2-p_i)-2*(a5-i*p_i))*2*(1/q-1)*TJ2*TH2+
	      ((real_n-3*real_j-3)*(a5-i*p_i)-
	       2*(a6-i*i*p_i))*2*(1/q-1)*TJ2*TH2_false;
	 if (tem1<0) tem1=0;
	 Covar2=Covar2+tem1;
	 if (number==(real_n-real_j-2)/2.0)
	   tem1=((a2-p_i)*TH2+(a5-i*p_i)*TH2_false)*2*(1/q-1)*TJ2+
	        ((a2-p_i)*TH1+(a5-i*p_i)*TH1_false)*2*(1/q-1)*TJ2+
		2*(1/q-1)*TJ2*p_i*(TH1+i*TH1_false);
	 else
	   tem1=(a2*TH2+a5*TH2_false)*2*(1/q-1)*TJ2+
	        (a2*TH1+a5*TH1_false)*2*(1/q-1)*TJ2;
	 if (tem1<0) tem1=0;
	 Covar2=Covar2+tem1;
	 tem1=(j*a2+a5)*2*(1/q-1)*TJ2*TH2+
	      (j*a5+a6)*2*(1/q-1)*TJ2*TH2_false;
	 if (tem1<0) tem1=0;
	 Covar2=Covar2+tem1;
       }

     m=real_n-2*real_j-2;
     if (m>=1)
        {
	  calculation_of_series(m,p);
          a2_m=a2; a5_m=a5; a6_m=a6; a7_m=a7;
	  a2=a2_m-a2_i;  a5=a5_m-a5_i; a6=a6_m-a6_i; a7=a7_m-a7_i;
	  tem=2*((0.5*real_n*(real_n-2*real_j-1)*a2-(real_n-real_j-0.5)*a5+
		  0.5*a6)*TJ2*TH2+
		 (0.5*real_n*(real_n-2*real_j-1)*a5-(real_n-real_j-0.5)*a6+
		  0.5*a7)*TJ2*TH2_false+
		((real_n-2*real_j-1)*a2-a5)*TJ2*TH1+
		 ((real_n-2*real_j-1)*a5-a6)*TJ2*TH1_false); 
	  if (tem<0) tem=0;
          Covar=Covar+tem;
          tem=((3*real_j+1-real_n)*(real_n-real_j-1)*a2+
	       (3*real_n-5*real_j-3)*a5-2*a6)*TJ2*TH2+
              ((3*real_j+1-real_n)*(real_n-real_j-1)*a5+
	       (3*real_n-5*real_j-3)*a6-2*a7)*TJ2*TH2_false+
	      (2*(3*real_j+1-real_n)*a2+4*a5)*TJ2*TH1+
	      (2*(3*real_j+1-real_n)*a5+4*a6)*TJ2*TH1_false;
	  if (tem<0) tem=0;
           Covar=Covar+tem;


          /* calculate the covariances between the reinforced variables */

	  tem1=((real_n-2*real_j-2)*a2-a5)*(1/q-1)*TJ2*TH2+
	       ((real_n-2*real_j-2)*a5-a6)*(1/q-1)*TJ2*TH2_false+
	       (a2*TH1+a5*TH1_false)*(1/q-1)*TJ2;
	  if (tem1<0) tem1=0;
	  Covar2=Covar2+2*tem1;
	}

     r=real_n-real_j-1;
     calculation_of_series(r,p);
     a2_r=a2; a5_r=a5; a6_r=a6; a7_r=a7;
     a2=a2_r-a2_m; a5=a5_r-a5_m; a6=a6_r-a6_m; a7=a7_r-a7_m;
     tem=((real_n-real_j-1)*a2-a5)*(real_n-real_j-1)*TJ2*TH2+
         ((real_n-real_j-1)*a5-a6)*(real_n-real_j-1)*TJ2*TH2_false+
	 2*(real_n-real_j-1)*a2*TJ2*TH1+2*(real_n-real_j-1)*a5*TJ2*TH1_false;
     if (tem<0) tem=0;
     Covar=Covar+tem;
   }

real_n=0;
real_j=0;
if (loop==n)
  {
    Covar3=0;                       /* 0T vs 2T & 1T  */
    for (j=k; j<=n-1; j++)
      {
	real_n=n;
	real_j=j;
	tem=0;
	tem=2*real_n*real_j*pow(q, 2)*pow(p, real_n+real_j)*(real_n-real_j-1);
	if (tem<0) tem=0;
	Covar3=Covar3+tem;
	tem=2*2*real_n*real_j*q*pow(p, real_n+real_j);
	if (tem<0) tem=0;
	Covar3=Covar3+tem;
      }
    if (Covar3<0) Covar3=0;
    tem=pow(p, real_n);
    if (tem<0) tem=0;
    Mean=Mean+real_n*tem;
    Var=Var+real_n*real_n*tem*(1-tem);
  }
 result.exp=Mean;
 result.var=Var+Covar2-Covar3-Covar;
 return(result);

}















