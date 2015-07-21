#include "iRRAM.h"
using namespace iRRAM;
void compute(){
  REAL a0=REAL(11)/2;
  REAL a1=REAL(61)/11;
  REAL ans;
  for(int i=2; i<1000; i++){
    ans = REAL(111)-(REAL(1130)-REAL(3000)/a0)/a1;
    a0=a1;
    a1=ans;
  }
  cout<<ans<<std::endl;
}
