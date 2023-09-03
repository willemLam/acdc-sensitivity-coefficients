function [E_star_augmented, YY_augmented, Grid_para_augmented] = Augment_system_with_filter_n_losses(E_star, Grid_para, Filter_para, idx1)


    %% Augment YY to include the filter and IGBT losses into the admitance matrix
    type = 0; %equivalent PI with filter and losses
    [Yac_augmented, A_augmented, Zloss,Zfilter,linedata1] = include_losses_filter_in_Y('linedata_AC.txt',Grid_para,Filter_para,idx1,E_star,type);
    Yac_augmented = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac_augmented,'UniformOutput',false));
    YY_augmented = blkdiag(Yac_augmented,Grid_para.Ydc);

    Grid_para_augmented = Grid_para;

    Grid_para_augmented.YY = YY_augmented;
    Grid_para_augmented.G = real(YY_augmented);
    Grid_para_augmented.B = imag(YY_augmented);
    Grid_para_augmented.Yac = Yac_augmented;
    Grid_para_augmented.Ydc = Grid_para.Ydc;
    Grid_para_augmented.n_ac = Grid_para.n_ac+Grid_para_augmented.n_AFE; %Every AFE gets one extra virtual node
    Grid_para_augmented.n_nodes = Grid_para_augmented.n_ac*Grid_para_augmented.n_ph + Grid_para_augmented.n_dc;


    Zloss = cell2mat(arrayfun(@(x) x*ones(Grid_para.n_ph,1),Zloss,'UniformOutput',false));
    Zfilter = cell2mat(arrayfun(@(x) x*ones(Grid_para.n_ph,1),Zfilter,'UniformOutput',false));

    idx3 = Get_multiphase_Node_indices(idx1,Grid_para);
    E_augment_loss = E_star(sort([idx3.vscac_pq;idx3.vscac_vq])) + (Zloss+Zfilter).*(Grid_para.YY(sort([idx3.vscac_pq;idx3.vscac_vq]),:)*E_star);
    E_star_augmented = [E_star(1:Grid_para.n_ac*Grid_para.n_ph); E_augment_loss; E_star(Grid_para.n_ac*Grid_para.n_ph+1 : end)  ];
    
    Grid_para_augmented.pos_ac3 = Grid_para.pos_ac3 + Grid_para.n_AFE*Grid_para.n_ph;
    Grid_para_augmented.pos_dc3 = Grid_para.pos_dc3 + Grid_para.n_AFE*Grid_para.n_ph;
end