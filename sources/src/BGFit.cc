// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *               BGFit.cc                 *
// *             Jonas Kasper              *
// *     kasper@physik.rwth-aachen.de      *
// *           Created in 2018             *
// *                                       *
// *****************************************

#include "BGFit.hh"
#include "TF1.h"

//------------------------------------------------------------------

BGFit::BGFit():low_border(215),up_border(325){}

BGFit::BGFit(double low_b, double up_b):low_border(low_b),up_border(up_b){	
}

double BGFit::EvaluatePol3(double *x, double *par){
	if(x[0]>low_border && x[0]<up_border) TF1::RejectPoint();
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}

double BGFit::EvaluatePol2(double *x, double *par){
	if(x[0]>low_border && x[0]<up_border) TF1::RejectPoint();
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

double BGFit::EvaluateExpo(double *x, double *par){
	if(x[0]>low_border && x[0]<up_border) TF1::RejectPoint();
	return TMath::Exp(par[0]+par[1]*x[0]);
}