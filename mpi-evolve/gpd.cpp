//**********************************************************//
//                                                          //
//    Diquark Calculation of GPDs                           //
// ---------------------------------------------------------//
//                                                          //
//    Author: Brandon Kriesten                              //
//                                                          //
//    Calculates initial scale GPDs in the Diquark model    //
//    using helicity amplitude formalism. Parameters are    //
//    fitted using form factor elastic scattering data.     //
//                                                          //
//    16 total parameters per GPD, each parameter is        //
//    fitted using Cates data, following reference          //
//    Gonzalez,Goldstein,Kathuria,Liuti.                    //
//                                                          //
//    Reference: 10.1103/PhysRevC.88.065206                 //
//                                                          //
//    Input Variables:                                      //
//                                                          //
//    X = longitudinal variable                             //
//    zeta = x_Bj                                           //
//    t = 4-momentum transfer squared                       //
//        always negative                                   //
//                                                          //
//    Output Variables:                                     //
//                                                          //
//    hu = H_u in DGLAP region                              //
//    hd = H_d in DGLAP region                              //
//    eu = E_u in DGLAP region                              //
//    ed = E_d in DGLAP region                              //
//    hutil = H_u tilde in DGLAP region                     //
//    hdtil = H_d tilde in DGLAP region                     //
//    eutil = E_u tilde in DGLAP region                     //
//    edtil = E_d tilde in DGLAP region                     //
//                                                          //
//    hu_plus = H_u upper curve for error band in DGLAP     //
//              region                                      //
//    hd_plus = H_d upper curve for error band in DGLAP     //
//              region                                      //
//    eu_plus = E_u upper curve for error band in DGLAP     //
//              region                                      //
//    ed_plus = E_d upper curve for error band in DGLAP     //
//              region                                      //
//    hutil_plus = H_u tilde upper curve for error band in  //
//              DGLAP region                                //
//    hdtil_plus = H_d tilde upper curve for error band in  //
//              DGLAP region                                //
//    eutil_plus = E_u tilde upper curve for error band in  //
//              DGLAP region                                //
//    edtil_plus = E_d tilde upper curve for error band in  //
//              DGLAP region                                //
//                                                          //
//    contact: btk8bh@virginia.edu                          //
//                                                          //
//**********************************************************//

//Necessary includes
#include <iostream>
#include <cmath>
#include <fstream>

//Include Header File
#include "gpd.h"

//Necessary for output
using namespace std;

//**********************************************************//
//    DGLAP REGION PARAMETERS                               //
//----------------------------------------------------------//
//                                                          //
//    M = Proton Mass                                       //
//    m = Quark Mass                                        //
//    Mx = Diquark Mass                                     //
//    Ml = Diquark Form Factor Mass Cutoff                  //
//    alpha = Regge t-independent exponent                  //
//    alphap = Regge t-depenedent exponent                  //
//    p = Regge t-dependent exponent                        //
//    beta = 0 (this should be fixed if higher moments      //
//           need better fit                                //
//    a = 1 can be fitted later, modifies t-dependent       //
//           regge parameters                               //
//    N = Normalization constant                            //
//    MMp = Mass factor with zeta dependence                //
//    MM = Mass factor zeta-independent                     //
//    Xp = X' which is zeta dependent                       //
//    t0 = minimum t based on zeta                          //
//    dT = Transverse momentum transfer between initial and //
//         final proton states, fourier transform of        //
//         transverse spatial position of quarks.           //
//                                                          //
//**********************************************************//

//Calculation of GPD H for u-quarks fitted to Dirac form factor data
double gpdHu(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  //10.1103/PhysRevC.88.065206 Parameters
  double m = .420;
  double Mx = .604;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.448;
  double p = 0.620;
  //  double alphap = 1.814;
  //double p = 0.449;
  double N = 2.043;
  double beta = 0.;
  double a = 1.;
  double MMp = (X-zeta)/(1.-zeta)*M*M - ((X-zeta)/(1.-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1.-X))*Mx*Mx-Ml*Ml;
  double Hu = 0;
  double Xp = (X-zeta)/(1.-zeta);
  double t0 = -zeta*zeta*M*M/(1.-zeta);
  double dT = sqrt((t0-t)*(1.-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {

    double D0 = (1.-X)*MM-kT*kT;
    double D1 = (1.-Xp)*MMp-kT*kT-dT*dT*(1.-Xp)*(1.-Xp);
    double D2 = (1.-Xp)*kT*dT;
    double omx3 = (1.-X)*(1.-X)*(1.-X);
    double zeta2 = zeta*zeta;
    double omz2 = (1.-zeta)*(1.-zeta);
    double mas = m+M*X;
    double masp = m+M*Xp;
    double kt2 = kT*kT;
    double reg1 = pow(X,-alpha);
    double regp = -alphap*pow((1.-X),p);
    double reg2 = pow(X,(regp*t));
    double reg3 = pow(X,-beta*zeta*zeta/(1.-zeta)*t);
    Hu += -N*2*pi*(1.-zeta/2)*reg1*reg2*reg3*(omx3/omz2)*kT		\
      *(((mas*masp+kt2)*D1-2*D2*D2)/(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Hu*.001 + ((zeta*zeta)/(4*(1.-zeta)))*gpdEu(X,zeta,t);
}

double gpdHuplus(double X,double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .603;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.448;
  double p = 0.620;
  double beta = 0;
  double a = 1;
  double N = 2.043;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg = reg1*reg2;

  double alpha_err = 0.07*alpha;
  double N_err = 0.1*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = .0220;
  double dp = .0170;
  
  double hu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
		       +pow(reg_err*alpha_err,2)+pow(N_err,2));
  return gpdHu(X,zeta,t) + gpdHu(X,zeta,t)*hu_err/(reg);
  
}

double gpdHuminus(double X,double zeta,double t){

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .603;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 1.814;
  double p = .449;
  //double alphap = 2.448;
  //double p = 0.620;
  double beta = 0;
  double a = 1;
  double N = 2.043;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg = reg1*reg2;

  double alpha_err = 0.07*alpha;
  double N_err = 0*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = .0220;
  double dp = .0170;

  double hu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
		       +pow(reg_err*alpha_err,2)+pow(N_err,2));
  
  return gpdHu(X,zeta,t) - hu_err*gpdHu(X,zeta,t)/(reg);
  
}

//Calculation of GPD H for d-quarks fitted to Dirac form factor data
double gpdHd(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = 0.860;
  double alpha = .0317;
  double alphap = 2.209;
  double p = 0.658;
  /*double alphap = 1.139;
    double p = -0.113;*/
  double N = 1.570;
  double beta = 0;
  double a = 1;
  double MMp = (X-zeta)/(1.-zeta)*M*M - ((X-zeta)/(1.-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1.-X))*Mx*Mx-Ml*Ml;
  double Hd = 0;
  double Xp = (X-zeta)/(1.-zeta);
  double t0 = -zeta*zeta*M*M/(1.-zeta);
  double dT = sqrt((t0-t)*(1.-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double D0 = (1.-X)*MM-kT*kT;
    double D1 = (1.-Xp)*MMp-kT*kT-dT*dT*(1.-Xp)*(1.-Xp);
    double D2 = (1.-Xp)*kT*dT;
    double omx3 = (1.-X)*(1.-X)*(1.-X);
    double zeta2 = zeta*zeta;
    double omz2 = (1.-zeta)*(1.-zeta);
    double mas = m+M*X;
    double masp = m+M*Xp;
    double kt2 = kT*kT;
    double reg1 = pow(X,-alpha);
    double regp = -alphap*pow(1.-X,p);
    double reg2 = pow(X,(regp*t));
    double reg3 = pow(X,-beta*zeta*zeta/(1.-zeta)*t);

    Hd += -N*2*pi*(1.-zeta/2)*reg1*reg2*reg3*(omx3/omz2)*kT\
      *(((mas*masp+kt2)*D1-2*D2*D2)/(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Hd*.001 + ((zeta*zeta)/(4*(1-zeta)))*gpdEd(X,zeta,t);
}

double gpdHdplus(double X, double zeta, double t) {

  
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = 0.860;
  double alpha = .0317;
  //double alphap = 1.139;
  //double p = -0.113;
  double alphap = 2.209;
  double p = 0.658;
  double beta = 0;
  double a = 1;
  double N = 1.570;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hd = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  
  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);
  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0.1*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = 0.0564;
  double dp = 0.1040;

  double hd_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHd(X,zeta,t) + gpdHd(X,zeta,t)*hd_err/(reg);
    
}

double gpdHdminus(double X, double zeta, double t) {


  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = 0.860;
  double alpha = .0317;
  double alphap = 1.139;
  double p = -0.113;
  //double alphap = 2.209;
  //double p = 0.658;
  double beta = 0;
  double a = 1;
  double N = 1.570;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hd = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);
  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = 0.0564 ;
  double dp = 0.1040 ;

  double hd_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHd(X,zeta,t) - gpdHd(X,zeta,t)*hd_err/(reg);

}

//Calculation of GPD E for u-quarks fitted to Pauli form factor data
double gpdEu(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .604;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.835;
  double p = .969;
  double beta = 0.;
  double a = 1;
  double N = 1.803;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Eu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double omx3 = (1-X)*(1-X)*(1-X);
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
    double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    Eu += N*2*pi*(1-zeta/2)*reg1*reg2*reg3*((omx3)/(1-zeta))	\
      *kT*((2*M*((-2*M)*(X-Xp)*kT*kT-(m+M*X)*(1-Xp)*D1))\
	   /(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Eu*.001;
}

double gpdEuplus(double X, double zeta, double t) {

  
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .604;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.835;
  double p = .969;
  double beta = 0;
  double a = 1;
  double N = 1.803;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Eu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0.11*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = 0.0509;
  double dp = 0.0307;
  
  double eu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEu(X,zeta,t) + gpdEu(X,zeta,t)*eu_err/reg;
  
}

double gpdEuminus(double X, double zeta, double t) {


  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .604;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.835;
  double p = .969;
  double beta = 0;
  double a = 1;
  double N = 1.803;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Eu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0.11*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double dap = 0.0509;
  double dp = 0.0307;
  
  double eu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEu(X,zeta,t) - gpdEu(X,zeta,t)*eu_err/reg;

}


//Calculation of GPD E for d-quarks fitted to Pauli form factor data
double gpdEd(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = .860;
  double alpha = .0317;
  double alphap = 1.281;
  double p = 0.726;
  double beta = 0;
  double a = 1;
  double N = -2.780;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Ed = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double omx3 = (1-X)*(1-X)*(1-X);
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
    double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    Ed += N*2*pi*(1-zeta/2)*reg1*reg2*reg3*((omx3)/(1-zeta))	\
      *kT*((2*M*((-2*M)*(X-Xp)*kT*kT-(m+M*X)*(1-Xp)*D1))\
	   /(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Ed*.001;
}

double gpdEdplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = .860;
  double alpha = .0317;
  double alphap = 1.281;
  double p = 0.726;
  double beta = 0;
  double a = 1;
  double N = -2.780;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Ed = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.1*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double dap = 0.0310;
  double dp = 0.0631;

  double ed_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEd(X,zeta,t) + gpdEd(X,zeta,t)*ed_err/reg;
  
}


double gpdEdminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = .860;
  double alpha = .0317;
  double alphap = 1.281;
  double p = 0.726;
  double beta = 0;
  double a = 1;
  double N = -2.780;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Ed = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.1*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double dap = 0.0310;
  double dp = 0.0631;

  double ed_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEd(X,zeta,t) - gpdEd(X,zeta,t)*ed_err/reg;

}

//Calculation of GPD H_t for u-quarks fitted to axial form factor data
double gpdHutil(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 1.543;
  double p = .346;
  double N = .0504;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
    Hutil += -N*2*pi*(1-zeta/2)*reg1*reg2*((pow((1-X),3))\
      /(pow((1-zeta),2)))*kT*((((m+M*X)*(m+M*Xp)-kT*kT)*D1+2*D2*D2)\
      /(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  if(zeta == 0) {
    return Hutil*0.001;
  }
  else {
    return Hutil*.001 + ((zeta*zeta)/(4*(1-zeta)))*gpdEutil(X,zeta,t);
  }
}

double gpdHutilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 1.543;
  double p = .346;
  double N = .0504;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHutil(X,zeta,t) + gpdHutil(X,zeta,t)*hutil_err/reg;

}

double gpdHutilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 1.543;
  double p = .346;
  double N = .0504;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHutil(X,zeta,t) - gpdHutil(X,zeta,t)*hutil_err/reg;

}


//Calculation of GPD H for d-quarks fitted to axial form factor data
double gpdHdtil(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 1.298;
  double p = .974;
  double N = -.0262;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hdtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
    Hdtil += -N*2*pi*(1-zeta/2)*reg1*reg2*((pow((1-X),3))\
      /(pow((1-zeta),2)))*kT*((((m+M*X)*(m+M*Xp)-kT*kT)*D1+2*D2*D2)\
      /(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  if(zeta == 0) {
    return Hdtil*.001;
  }
  else{
    return Hdtil*.001 + ((zeta*zeta)/(4*(1-zeta)))*gpdEutil(X,zeta,t);
  }
}

double gpdHdtilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 1.298;
  double p = .974;
  double N = -.0262;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hdtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hdtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHdtil(X,zeta,t) + gpdHdtil(X,zeta,t)*hdtil_err/reg;
  
}

double gpdHdtilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 1.298;
  double p = .974;
  double N = -.0262;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hdtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hdtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHdtil(X,zeta,t) - gpdHdtil(X,zeta,t)*hdtil_err/reg;

}

//Calculation of GPD E_t for u-quarks fitted to pseudo-scalar form factor data
double gpdEutil(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 5.130;
  double p = 3.507;
  double N = 1.074;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Eutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0.001; kT < 5; kT+=.001) {
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1.-X),p)*t));
    double xtil = (X-zeta) / (1.-zeta);
    /*
    double xtil1 = (X-zeta) / (1.-X);
    double emtil = xtil*M*M - Ml*Ml - xtil1*Mx*Mx;
    double emtilt = xtil*M*M-Ml*Ml-xtil1*Mx*Mx;
    double emtil0t = X*M*M-Ml*Ml-Mx*Mx*X/(1.-X);
    double antile1 = -4*M*((m+M*X)+(m+M*xtil))*kT;
    double antile2 = -4*M*(m+M*X)*(1.-X)/(1.-zeta);
    double aaat = emtilt - kT*kT*((1.-zeta)/(1.-X))-dT*dT*(1.-X)/(1.-zeta);
    double bbb = dT*kT*2;
    double akntile1 = antile1*2*pi/pow(aaat*aaat-bbb*bbb,1.5)*kT*2;
    double akntile2 = antile2*2*pi/pow(aaat*aaat - bbb*bbb,1.5)*aaat;
    double denomt = pow(emtil0t-kT*kT/(1.-X),2);
    Eutil += N*reg1*reg2*(1./(1.-X))*kT*(1.-zeta)*(akntile1+akntile2)/denomt*(1/zeta);
    */    
    double cM = X*(1-X)*M*M - (1-X)*Ml*Ml - Mx*Mx*X;
    double cMp = Xp*(1-Xp)*M*M-(1-Xp)*Ml*Ml-Mx*Mx*Xp;
    double D = cM-kT*kT;
    double aaa = cMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double bbb = 2*kT*dT*(1-Xp);
    double num = 4*M*2*kT*kT*((m+M*X)+(m+M*Xp))-aaa*(4*M*(m+M*X))*(1-Xp);
    double denom = D*D*pow(aaa*aaa-bbb*bbb,1.5);

    Eutil += 2*pi*N*((1-zeta)/zeta)*kT*(1/(1-X))*reg1*reg2*(num/denom)*(1-X)*(1-X)*(1-Xp)*(1-Xp)*(1-Xp);
  }
  
  return Eutil*.001;
  
}

double gpdEutilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 5.130;
  double p = 3.507;
  double N = 1.074;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Eutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double eutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEutil(X,zeta,t) + gpdEutil(X,zeta,t)*eutil_err/reg;
  
}

double gpdEutilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 5.130;
  double p = 3.507;
  double N = 1.074;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Eutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double eutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEutil(X,zeta,t) - gpdEutil(X,zeta,t)*eutil_err/reg;

}



//Calculation of GPD E_t for d-quarks fitted to pseudo-scalar form factor data
double gpdEdtil(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 3.385;
  double p = 2.326;
  double N = -0.966;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Edtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1.-X),p)*t));
    double xtil = (X-zeta) / (1.-zeta);

    /*
    double xtil1 = (X-zeta) / (1.-X);
    double emtil = xtil*M*M - Ml*Ml - xtil1*Mx*Mx;
    double emtilt = xtil*M*M-Ml*Ml-xtil1*Mx*Mx;
    double emtil0t = X*M*M-Ml*Ml-Mx*Mx*X/(1.-X);
    double antile1 = -4*M*((m+M*X)+(m+M*xtil))*kT;
    double antile2 = -4*M*(m+M*X)*(1.-X)/(1.-zeta);
    double aaat = emtilt - kT*kT*((1.-zeta)/(1.-X))-dT*dT*(1.-X)/(1.-zeta);
    double bbb = dT*kT*2;
    double akntile1 = antile1*2*pi/pow(aaat*aaat-bbb*bbb,1.5)*kT*2;
    double akntile2 = antile2*2*pi/pow(aaat*aaat - bbb*bbb,1.5)*aaat;
    double denomt = pow(emtil0t-kT*kT/(1.-X),2);
    Edtil += N*reg1*reg2*(1./(1.-X))*kT*(1.-zeta)*(akntile1+akntile2)/denomt * (1/zeta);
    */


    double cM = X*(1-X)*M*M - (1-X)*Ml*Ml - Mx*Mx*X;
    double cMp = Xp*(1-Xp)*M*M-(1-Xp)*Ml*Ml-Mx*Mx*Xp;
    double D = cM-kT*kT;
    double aaa = cMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double bbb = 2*kT*dT*(1-Xp);
    double num = 4*M*2*kT*kT*((m+M*X)+(m+M*Xp))-aaa*(4*M*(m+M*X))*(1-Xp);
    double denom = D*D*pow(aaa*aaa-bbb*bbb,1.5);

    Edtil += 2*pi*N*((1-zeta)/zeta)*kT*(1/(1-X))*reg1*reg2*(num/denom)*(1-X)*(1-X)*(1-Xp)*(1-Xp)*(1-Xp);
  }

  return Edtil*.001;
  
}

double gpdEdtilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 3.385;
  double p = 2.326;
  double N = -0.966;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Edtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double edtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEdtil(X,zeta,t) + gpdEdtil(X,zeta,t)*edtil_err/reg;
  
}

double gpdEdtilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 3.385;
  double p = 2.326;
  double N = -0.966;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Edtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double edtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
		       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEdtil(X,zeta,t) - gpdEdtil(X,zeta,t)*edtil_err/reg;

}
