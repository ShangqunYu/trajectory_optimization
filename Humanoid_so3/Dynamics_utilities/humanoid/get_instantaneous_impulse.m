function [imp_toe,imp_heel,...
    v_prior_toe,v_prior_heel,...
    v_post_toe,v_post_heel] = get_instantaneous_impulse(q,qdot,model,is_opt)


[osim_toe,osim_heel] = get_operational_sp_inertia(model, q,is_opt);

[~,~,~,~,Jc_toe,Jc_heel] = get_toe_heel_jacobians( model, q, is_opt);

for ii = 1:model.NLEGS
    v_prior_toe{ii} = Jc_toe{ii}*qdot;
    v_prior_heel{ii} = Jc_heel{ii}*qdot;
    
    imp_toe{ii} = -osim_toe{ii}*v_prior_toe{ii};
    imp_heel{ii} = -osim_heel{ii}*v_prior_heel{ii};
    
    v_post_toe{ii} = v_prior_toe{ii} + osim_toe{ii}\imp_toe{ii};
    v_post_heel{ii} = v_prior_heel{ii} + osim_heel{ii}\imp_heel{ii};
end


    