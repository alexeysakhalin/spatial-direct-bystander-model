#ifndef IKM_TARGET_RESPONSE_MATH_H
#define IKM_TARGET_RESPONSE_MATH_H
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace PhysiCell {
struct IKM_ResponseStep { double first, second, mean_second; };
inline void ikm_unit_interval(double x) {
    if(!std::isfinite(x) || x<0 || x>1) throw std::runtime_error("Invalid IKM response fraction");
}
inline double ikm_hill_occupancy(double concentration,double half,double power) {
    if(!std::isfinite(concentration) || !std::isfinite(half) || !std::isfinite(power)
       || concentration < -1e-12 || half<=0 || power<=0)
        throw std::runtime_error("Invalid IKM cytokine input");
    if(concentration<=0) return 0;
    // A logistic form avoids overflowing pow(c/K,n).
    double x=power*(std::log(concentration)-std::log(half));
    if(x>=0) return 1/(1+std::exp(-x));
    double y=std::exp(x); return y/(1+y);
}
inline IKM_ResponseStep ikm_two_stage_response(double first,double second,double input,
                                              double dt,double tau) {
    ikm_unit_interval(first); ikm_unit_interval(second); ikm_unit_interval(input);
    if(!std::isfinite(dt) || !std::isfinite(tau) || dt<=0 || tau<=0)
        throw std::runtime_error("Invalid IKM response timescale");
    double x=dt/tau;
    if(!std::isfinite(x)) throw std::runtime_error("IKM response timestep overflow");
    if(x==0) return {first,second,second}; // positive ratio below machine range
    double e=std::exp(-x);
    double a=-std::expm1(-x)/x;
    double b, w;
    if(x<1e-4) {
        b=x*(.5-x/3+x*x/8-x*x*x/30+x*x*x*x/144);
        w=x*x*(1./6-x/12+x*x/40-x*x*x/180+x*x*x*x/1008);
    } else { b=(-std::expm1(-x)-x*e)/x; w=1-a-b; }
    IKM_ResponseStep r;
    r.first=input+(first-input)*e;
    r.second=input+(second-input+(first-input)*x)*e;
    r.mean_second=input*w+second*a+first*b;
    // Restrict floating-point roundoff only; invalid input was rejected above.
    r.first=std::max(0.,std::min(1.,r.first));
    r.second=std::max(0.,std::min(1.,r.second));
    r.mean_second=std::max(0.,std::min(1.,r.mean_second));
    return r;
}
}
#endif
