#include <iostream>
using namespace std;
int main(){
  double a0=double(11)/2;
  double a1=double(61)/11;
  double ans;
  for(int i=2; i<1000; i++){
    ans = double(111)-(double(1130)-double(3000)/a0)/a1;
    a0=a1;
    a1=ans;
  }
  cout<<ans<<std::endl;
}
