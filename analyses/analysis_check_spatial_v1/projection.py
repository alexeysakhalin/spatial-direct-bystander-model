"""Select portable numerical results while preserving all hourly observations."""
def project(case, access, gradients):
    akeys=['rows','first_saved_positive_hours']
    envelope={k:access['initial_target_envelope'][k]for k in ['dimension','halfspace_equations','numerical_boundary_tolerance_um']}
    gkeys=['rows','grid_xyz','spacing_xyz_um','D_um2_min','loss_per_min','decay_only_half_life_min','homogeneous_diffusion_loss_length_um']
    return {'case':case['id'],'ratio':case['ratio'],'arm':case['arm'],'seed':case['seed'],'access':dict({k:access[k]for k in akeys},initial_target_envelope=envelope),'gradients':{k:gradients[k]for k in gkeys}}
