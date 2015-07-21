// Evaluate a function analytic on a rectangle
// by using iterated analytic continuation
#include "iRRAM.h"
#include "Functions/Analytic/BA_ANA.h"
REAL factorial(const int n){
	if(n==0) return REAL(1);
  static vector<REAL> ans;
  //if(ans.size() > n)
  //  return ans[n];
  ans.resize(n+1);
	return REAL(n)*factorial(n-1);
  return ans[n];
}


REAL inv_factorial(const int n){
	if ((n!=0)&&((n+1)*n+2*n>= -ACTUAL_STACK.actual_prec)){
		REAL return_value(0);
		sizetype error;
		sizetype_set(error,1,ACTUAL_STACK.actual_prec);
		return_value.seterror(error);
		return return_value;
	}
	if (n==0)
		return REAL(1);
	REAL inv_fact=inv_factorial(n-1)/REAL(n);
	return inv_fact;
}


// evaluate the d-th derivative of f at some point z.
// Requirement for the error estimate: |z|<=1/(2l)
template<class ARG>
ARG eval_derivative(shared_ptr<const POWERSERIES<ARG>> series, const int l,
    const REAL& B,const REAL& center,  const ARG& z, const int d) {
  // first compute the parameters B, l 
	ARG sum(series->get_coeff(d)); // partial sum
	sizetype sum_error; // and its error
	sum.geterror(sum_error); // get error of first coefficient
	ARG factor(1); // (i+n ncr n)*x^i
  // Result with the smallest error so far will be returned in the end
  ARG best(sum);
  sizetype best_error;
  // error bound approximating the error made by considering
  // only J many coefficients (truncation error)
  // it is given by B*l^d*2^(-J+1)*d*choose(J+d,d)
  REAL error_factor = B*power(l,d)*2*d;
  // initial estimate for the number of terms to be considered
  int j=0,J=d+l-ACTUAL_STACK.actual_prec; 
  sizetype trunc_error; // error bound for the rest sum
  REAL error = error_factor;
  bool first_iterate=true;
  do { // Calculate the new partial sum
       // and error=error_factor*choose(J+d,d)/2^J:
       // iteratively to avoid excessively long intermediate results
     while (j++<J){
        error*=REAL(j+d)/REAL(2*j);
        factor = factor*z*REAL(j+d)/REAL(j);
        sum = sum+series->get_coeff(j+d)*factor; 
        }
     sizetype_add(trunc_error,error.vsize,error.error);
     sum.geterror(sum_error);
     sizetype local_error; // this holds the error for the current sum
     sizetype_add(local_error, sum_error, trunc_error);
     // check whether the local error is better than the best one.
     if (first_iterate || sizetype_less(local_error, best_error)) { 
       first_iterate = false;
       best = sum; // if so, adjust the return value.
       best_error = local_error;
       best.seterror(best_error); // and its error. 
       }
    J+=d;  // Next time try more terms.
    } // Continue until achieving sufficent precision
      // or until more terms actually INcrease the error
    while (sizetype_less(sum_error, trunc_error) &&
      (best_error.exponent >= ACTUAL_STACK.actual_prec));
  return best;
}

template<class ARG>
BA_ANA<ARG> analytic_continuation(const BA_ANA<ARG>& f, const ARG& z){
   // define poweseries by closure
   POWERSERIES<ARG> ps([f,z] 
                          (unsigned int i) -> ARG {
                            return eval(f->get_coeff(), f->get_k(), f->get_A(),z,i);
                          });
   return BA_ANA<ARG>(ps, f.get_k(), f.get_A(),1);
}


REAL sinseries(const int n){
	if (0 == n%2)
		return 0;
	else {
		if (0 == (n-1)%4)
			return inv_factorial(n);
		else
			return -inv_factorial(n);
	 }
}

void compute(){
  BA_ANA<REAL> real_test(sinseries, 2, 2, 1);
  REAL x = "0.05";
  int n;
  iRRAM::cin >> n;
  REAL z = "0.1";
  iRRAM::cout <<setRwidth(50)   << eval(real_test, z, 0) << endl;
  for(int i=0; i<n; i++)
  {
    real_test = analytic_continuation(real_test, z);
  }
  REAL real_x = x+n*REAL(0.1);
  iRRAM::cout << "After " << n << " continuations; evaluating at " 
              << real_x  << endl;
  rwrite(eval(real_test, x, 0),20);
  iRRAM::cout << endl;
  iRRAM::cout << "should be:" << endl;
  rwrite(sin(real_x), 20);
  iRRAM::cout << endl;
}

