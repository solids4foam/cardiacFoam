#include <cmath>
#include <cstdio>
typedef double scalar;
namespace constant { namespace mathematical { const double pi = M_PI; } }
int main() {
  struct P { double X,Y,Z,t; };
  P pts[3] = {{0.2,0.3,0.4,0.1},{0.55,0.65,0.45,0.05},{0.8,0.15,0.7,0.08}};
  const double Ax=0.02,Ay=0.02,Az=0.02,V0=1.0,gamma=1.0,mu=3846.15,K=8333.33,rho0=1060.0;
  const char* names[2] = {"PASSIVE_Tmax0","FULL_Tmax1000"};
  double Tmaxs[2] = {0.0, 1000.0};
  for (int blk=0; blk<2; ++blk) {
    double Tmax = Tmaxs[blk];
    for (int i=0;i<3;++i) {
      double X=pts[i].X,Y=pts[i].Y,Z=pts[i].Z,t=pts[i].t;
      #include "B_expr.H"
      printf("CPP %s X=%.2f Y=%.2f Z=%.2f t=%.2f  B=(% .10e % .10e % .10e)\n",
             names[blk],X,Y,Z,t,Bx,By,Bz);
    }
  }
  return 0;
}
