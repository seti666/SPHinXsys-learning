#pragma once
// Corotating vector penalty: e = particle - bone-carried target.
// Rayleigh dissipation uses e_dot - omega cross e (objective rate).
#include <cmath>
struct AttachmentLawResult { double fx, fy, torque, energy, dissipation, power_residual; };
inline AttachmentLawResult attachmentLaw(double x,double y,double vx,double vy,
 double tx,double ty,double ox,double oy,double ovx,double ovy,double omega,double k,double c) {
 const double ex=x-tx, ey=y-ty;
 const double rx=x-ox, ry=y-oy;
 const double wx=vx-ovx+omega*ry, wy=vy-ovy-omega*rx;
 const double fx=-k*ex-c*wx, fy=-k*ey-c*wy;
 const double torque=-(rx*fy-ry*fx);
 const double edx=vx-ovx+omega*(ty-oy), edy=vy-ovy-omega*(tx-ox);
 const double diss=c*(wx*wx+wy*wy);
 const double power=fx*vx+fy*vy-fx*ovx-fy*ovy+torque*omega;
 return {fx,fy,torque,0.5*k*(ex*ex+ey*ey),diss,power+k*(ex*edx+ey*edy)+diss};
}
