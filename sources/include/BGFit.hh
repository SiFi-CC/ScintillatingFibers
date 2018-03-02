// *****************************************
// *                                       *
// *          ScintillatingFibers          *
// *              BGFit.hh                  *
// *            Jonas Kasper               *
// *    kasper@physik.rwth-aachen.de       *        
// *          Created in 2018              *
// *                                       *
// *****************************************


#ifndef __BGFit_H_
#define __BGFit_H_ 1


class BGFit {
 private: 
	double low_border;
	double up_border;

 public:
	BGFit();
	BGFit(double low_b, double up_b);
	double EvaluatePol3(double *x, double *par);
	double EvaluatePol2(double *x, double *par);
	double EvaluateExpo(double *x, double *par);
  
};

#endif
