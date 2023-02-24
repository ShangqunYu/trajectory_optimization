function [osim_toe,osim_heel] = get_operational_sp_inertia(model, q,is_opt)

[H,~] = get_mass_matrix( model, q, is_opt);
[~,~,~,~,Jc_toe,Jc_heel] = get_toe_heel_jacobians( model, q, is_opt);

if is_opt==1
    H_inv = inv(H,'symbolicqr');
else
    H_inv = inv(H);

end
for ii=1:model.NLEGS
    if is_opt==1
        osim_toe{ii} = inv(Jc_toe{ii}*H_inv*Jc_toe{ii}','symbolicqr');%operational space inertia matrix
        osim_heel{ii} = inv(Jc_heel{ii}*H_inv*Jc_heel{ii}','symbolicqr');%operational space inertia matrix
    else
        osim_toe{ii} = inv(Jc_toe{ii}*H_inv*Jc_toe{ii}');%operational space inertia matrix
        osim_heel{ii} = inv(Jc_heel{ii}*H_inv*Jc_heel{ii}');%operational space inertia matrix
    end
end