#include "../attachment_law.h"
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
void require(bool v){if(!v) throw std::runtime_error("attachment law test failed");}
int main(){
 std::mt19937 g(17);std::uniform_real_distribution<double>d(-1,1);
 double worst=0;
 for(int i=0;i<10000;i++){
  double x=d(g),y=d(g),tx=d(g),ty=d(g),ox=d(g),oy=d(g),vx=d(g),vy=d(g),ux=d(g),uy=d(g),w=d(g),k=3,c=.7;
  auto a=attachmentLaw(x,y,vx,vy,tx,ty,ox,oy,ux,uy,w,k,c);
  require(a.dissipation>=0);worst=std::max(worst,std::abs(a.power_residual));require(std::abs(a.power_residual)<1e-12);
  require(std::abs(x*a.fy-y*a.fx + a.torque-ox*a.fy+oy*a.fx)<1e-12);
  // Adding a common rigid velocity must not change forces or dissipation.
  double W=.8,U=.3,V=-.2;
  auto b=attachmentLaw(x,y,vx+U-W*y,vy+V+W*x,tx,ty,ox,oy,ux+U-W*oy,uy+V+W*ox,w+W,k,c);
  require(std::abs(a.fx-b.fx)<1e-12 && std::abs(a.fy-b.fy)<1e-12);
  // Independent finite-difference derivative of the spring energy under bone rotation.
  auto energy=[&](double angle){double rx=tx-ox,ry=ty-oy;double X=ox+cos(angle)*rx-sin(angle)*ry,Y=oy+sin(angle)*rx+cos(angle)*ry;return .5*k*((x-X)*(x-X)+(y-Y)*(y-Y));};
  auto elastic=attachmentLaw(x,y,0,0,tx,ty,ox,oy,0,0,0,k,0);
  double deriv=(energy(1e-6)-energy(-1e-6))/2e-6;
  require(std::abs(elastic.torque+deriv)<1e-7);
 }
 auto rest=attachmentLaw(.2,.3,0,0,.2,.3,0,0,0,0,0,3,.7);
 require(rest.fx==0 && rest.fy==0 && rest.torque==0);
 std::cout<<"PASS 10000 randomized states: angular balance, power identity, rigid-motion objectivity, finite-difference energy torque; max power residual="<<worst<<"\n";
}

