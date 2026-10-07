#ifndef IKM_REFLECTIVE_DOMAIN_H
#define IKM_REFLECTIVE_DOMAIN_H
#include <cmath>
#include <stdexcept>
namespace PhysiCell {
// Mirror the cell centre into the numerical box, without a new wall-force parameter.
inline double ikm_reflect_coordinate(double x, double lo, double hi, double& sign)
{
    if(!std::isfinite(x) || !std::isfinite(lo) || !std::isfinite(hi) || hi<=lo)
        throw std::runtime_error("Invalid reflecting-domain coordinate");
    sign=1.0;
    if(x>=lo && x<hi) return x;
    const double length=hi-lo;
    double phase=std::fmod(x-lo,2.0*length);
    if(phase<0) phase+=2.0*length;
    double answer;
    if(phase<length) answer=lo+phase;
    else {answer=hi-(phase-length);sign=-1.0;}
    // Keep voxel lookup away from the upper endpoint.
    return std::min(std::nextafter(hi,lo),std::max(lo,answer));
}
inline bool ikm_reflect_position(Cell* c)
{
    int enabled=c->custom_data.find_variable_index("ikm_reflect_domain");
    if(enabled<0 || c->custom_data[enabled]<0.5) return false;
    const auto& box=c->get_container()->underlying_mesh.bounding_box;
    for(int axis=0;axis<3;++axis)
    {
        if(axis==2 && BioFVM::default_microenvironment_options.simulate_2D)
        {
            if(box[2]>0 || box[5]<=0) throw std::runtime_error("2D reflecting domain must contain z=0");
            c->position[2]=0; c->velocity[2]=0;
            continue;
        }
        double sign;
        c->position[axis]=ikm_reflect_coordinate(c->position[axis],box[axis],box[axis+3],sign);
        c->velocity[axis]*=sign;
    }
    return true;
}
}
#endif
