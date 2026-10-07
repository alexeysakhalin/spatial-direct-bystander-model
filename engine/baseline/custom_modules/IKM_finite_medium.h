#ifndef IKM_FINITE_MEDIUM_H
#define IKM_FINITE_MEDIUM_H
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include "../BioFVM/BioFVM.h"
#include "../modules/PhysiCell_settings.h"

namespace IKM {
// Finite-volume reservoir exchange. No biological response parameters or random draws.
class FiniteMedium {
public:
    struct Species {
        std::string key;
        int substrate;
        double reservoir, external, decay, initial_amount, last_domain;
        double native_delta=0, exported=0, reservoir_lost=0, max_residual=0;
        std::vector<int> nodes;
        std::vector<double> volume, conductance;
    };
    bool enabled=false;
    double reservoir_volume=0, external_flow=0;
    BioFVM::Microenvironment* environment=nullptr;
    std::vector<Species> species;
    std::ofstream ledger;

    static double required(const std::string& key,bool positive=false) {
        auto& d=PhysiCell::parameters.doubles;
        if(d.find_index(key)<0) throw std::runtime_error("Missing medium parameter: "+key);
        double v=d(key);
        if(!std::isfinite(v) || v<0 || (positive && v<=0))
            throw std::runtime_error("Invalid medium parameter: "+key);
        return v;
    }
    double domain_amount(int substrate) const {
        long double total=0;
        for(unsigned i=0;i<environment->number_of_voxels();++i) {
            double c=environment->density_vector(i)[substrate];
            if(!std::isfinite(c) || c<0) throw std::runtime_error("Invalid medium concentration");
            total+=(long double)environment->mesh.voxels[i].volume*c;
        }
        return (double)total;
    }
    void configure(BioFVM::Microenvironment& m,bool write_ledger=true) {
        auto& b=PhysiCell::parameters.bools;
        if(b.find_index("ikm_medium_enabled")<0 || !b("ikm_medium_enabled")) return;
        enabled=true;environment=&m;
        reservoir_volume=required("ikm_medium_reservoir_volume_um3",true);
        external_flow=required("ikm_medium_external_flow_um3_per_min");
        const double ell=required("ikm_medium_external_distance_um");
        const bool two=BioFVM::default_microenvironment_options.simulate_2D;
        auto& mesh=m.mesh;
        if(two && mesh.z_coordinates.size()!=1) throw std::runtime_error("Medium2D requires one z layer");
        for(int n=0;n<2;++n) {
            Species s;s.key=n==0?"ifng":"tnf";
            s.substrate=m.find_density_index(n==0?"IFN-gamma":"TNF");
            if(s.substrate<0) throw std::runtime_error("Missing medium substrate");
            s.reservoir=required("ikm_medium_"+s.key+"_initial");
            s.external=required("ikm_medium_"+s.key+"_external");
            s.decay=required("ikm_medium_"+s.key+"_decay_per_min");
            const double D=m.diffusion_coefficients[s.substrate];
            if(!std::isfinite(D)||D<0) throw std::runtime_error("Invalid diffusion coefficient");
            for(unsigned i=0;i<m.number_of_voxels();++i) {
                if(m.is_dirichlet_node(i) && m.get_substrate_dirichlet_activation(s.substrate,i))
                    throw std::runtime_error("Medium substrate also has an active Dirichlet boundary");
                auto idx=mesh.cartesian_indices(i);
                double g=0;
                if(two) g=D*mesh.dx*mesh.dy/(ell+0.5*mesh.dz);
                else {
                    // Count both faces when an axis has only one voxel.
                    if(idx[0]==0) g+=D*mesh.dy*mesh.dz/(ell+0.5*mesh.dx);
                    if(idx[0]+1==mesh.x_coordinates.size()) g+=D*mesh.dy*mesh.dz/(ell+0.5*mesh.dx);
                    if(idx[1]==0) g+=D*mesh.dx*mesh.dz/(ell+0.5*mesh.dy);
                    if(idx[1]+1==mesh.y_coordinates.size()) g+=D*mesh.dx*mesh.dz/(ell+0.5*mesh.dy);
                    if(idx[2]==0) g+=D*mesh.dx*mesh.dy/(ell+0.5*mesh.dz);
                    if(idx[2]+1==mesh.z_coordinates.size()) g+=D*mesh.dx*mesh.dy/(ell+0.5*mesh.dz);
                }
                if(g>0) {s.nodes.push_back(i);s.volume.push_back(mesh.voxels[i].volume);s.conductance.push_back(g);}
            }
            s.last_domain=domain_amount(s.substrate);
            s.initial_amount=s.last_domain+reservoir_volume*s.reservoir;
            species.push_back(s);
        }
        if(write_ledger) {
            std::string name=PhysiCell::PhysiCell_settings.folder+"/medium_balance.csv";
            ledger.open(name);
            if(!ledger) throw std::runtime_error("Cannot open medium ledger");
            ledger<<"time_min,substrate,reservoir_concentration,domain_amount,reservoir_amount,native_net_change,external_net_export,reservoir_decay_loss,balance_residual,max_abs_residual\n";
            ledger<<std::setprecision(17);
        }
    }
    void advance(double dt) {
        if(!enabled) return;
        if(!std::isfinite(dt)||dt<0) throw std::runtime_error("Invalid medium timestep");
        for(auto& s:species) {
            const double before=domain_amount(s.substrate);
            s.native_delta+=before-s.last_domain;
            const double old_r=s.reservoir;
            s.reservoir*=std::exp(-s.decay*dt);
            s.reservoir_lost+=reservoir_volume*(old_r-s.reservoir);
            long double numerator=reservoir_volume*s.reservoir+dt*external_flow*s.external;
            long double denominator=reservoir_volume+dt*external_flow;
            for(unsigned k=0;k<s.nodes.size();++k) {
                double a=1/(1+dt*s.conductance[k]/s.volume[k]);
                double w=dt*s.conductance[k]*a;
                numerator+=(long double)w*environment->density_vector(s.nodes[k])[s.substrate];
                denominator+=w;
            }
            const double rnew=(double)(numerator/denominator);
            long double domain_change=0;
            for(unsigned k=0;k<s.nodes.size();++k) {
                double a=1/(1+dt*s.conductance[k]/s.volume[k]);
                double& c=environment->density_vector(s.nodes[k])[s.substrate];
                double next=a*c+(1-a)*rnew;
                domain_change+=(long double)s.volume[k]*(next-c);c=next;
            }
            s.exported+=dt*external_flow*(rnew-s.external);s.reservoir=rnew;
            s.last_domain=before+(double)domain_change;
            const double residual=s.last_domain+reservoir_volume*s.reservoir-s.initial_amount-s.native_delta+s.exported+s.reservoir_lost;
            s.max_residual=std::max(s.max_residual,std::abs(residual));
            const double scale=1+s.initial_amount+std::abs(s.native_delta)+std::abs(s.exported)+s.reservoir_lost;
            if(!std::isfinite(residual)||std::abs(residual)>1e-8*scale)
                throw std::runtime_error("Finite-medium amount balance failed");
        }
    }
    void write(double time) {
        if(!enabled||!ledger.is_open()) return;
        for(const auto& s:species) {
            const double current=domain_amount(s.substrate);
            const double residual=current+reservoir_volume*s.reservoir-s.initial_amount-s.native_delta+s.exported+s.reservoir_lost;
            ledger<<time<<','<<s.key<<','<<s.reservoir<<','<<current<<','<<reservoir_volume*s.reservoir<<','<<s.native_delta<<','<<s.exported<<','<<s.reservoir_lost<<','<<residual<<','<<s.max_residual<<'\n';
        }
        ledger.flush();
        if(!ledger) throw std::runtime_error("Medium ledger write failed");
    }
};
}
#endif
