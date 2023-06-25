function [Yac, A, Zloss,Zfilter,linedata1] = include_losses_filter_in_Y(Text,Grid_para,Filter_para,idx1,E_star,type)
%     Text = 'linedata_AC.txt';

    %% get filter and loss parameters
   
    [Zloss,Zfilter] = get_filter_loss_impedance(Grid_para,Filter_para,E_star);
%     Zloss = 0*Zloss + 1e-6;

    %% make Y augmented
    [Yac, YYL, YL, YT, YYT, ~, ~, ~, ~, A, linedata1] = make_Y_augmented(Text,Grid_para,idx1,Zloss,Zfilter,type);

    
%     %% make Y prime
%     [Yac_p, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = make_Y_prime(Text,Grid_para,idx1,Zloss);
%         %% make Y prime prime
%     [Yac_pp,~ , ~, ~, ~, ~, ~, ~, ~, ~, ~] = make_Y_prime(Text,Grid_para,idx1,Zloss+Zfilter);

end
