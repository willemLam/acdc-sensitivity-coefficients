% REDO EASIER 3 for loops!

function [IKp_abs, IKp, IKq_abs,IKq, IKv_abs, IKv, Currents]=Coeffs_Currents(E,K,Grid_para, idx)

YYL = Grid_para.YYL;
YYT = Grid_para.YYT;
nph = Grid_para.n_ph;

[COEFFpcomplex_alph, COEFFqcomplex_alph, COEFFvcomplex_alph] = transform_K_complex(K, Grid_para, idx);


    Currents = cell(size(YYL)/nph);
    IKp_abs = cell(1,size(COEFFpcomplex_alph,2));
    IKq_abs = cell(1,size(COEFFpcomplex_alph,2));
    IKp = cell(1,size(COEFFpcomplex_alph,2));
    IKq = cell(1,size(COEFFpcomplex_alph,2));
    
    % For all active/reactive power injections
    for l = 1:size(COEFFpcomplex_alph,2)
        % Re-initialize Temp Matrices
        temp_p = cell(size(YYL)/nph);
        temp_q = cell(size(YYL)/nph); 
        temp_v = cell(size(YYL)/nph); 
        temp_p_complex = cell(size(YYL)/nph);
        temp_q_complex = cell(size(YYL)/nph);
        temp_v_complex = cell(size(YYL)/nph);
        
        % For all rows of Admittance Matrix
        for i = 1:nph:size(YYL,1)
            for j = 1:nph:size(YYL,2)
                
                % Compute Currents only once
                if (l == 1)
                    Currents{ceil(i/nph),ceil(j/nph)} = YYL(i:i+(nph-1),j:j+(nph-1))*(E(i:i+(nph-1)) - E(j:j+(nph-1))) + ...
                                    YYT(i:i+(nph-1),j:j+(nph-1))*E(i:i+(nph-1));
                end
                
                if( max(abs(Currents{ceil(i/nph),ceil(j/nph)})) ~= 0)
                    % Compute Complex Coeffs
                    dIij_dPl = YYL(i:i+(nph-1),j:j+(nph-1))*(COEFFpcomplex_alph(i:i+(nph-1),l) - COEFFpcomplex_alph(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))* COEFFpcomplex_alph(i:i+(nph-1),l);
                    dIij_dQl = YYL(i:i+(nph-1),j:j+(nph-1))*(COEFFqcomplex_alph(i:i+(nph-1),l) - COEFFqcomplex_alph(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))* COEFFqcomplex_alph(i:i+(nph-1),l);
                    dIij_dVl = YYL(i:i+(nph-1),j:j+(nph-1))*(COEFFvcomplex_alph(i:i+(nph-1),l) - COEFFvcomplex_alph(j:j+(nph-1),l)) + ...
                               YYT(i:i+(nph-1),j:j+(nph-1))* COEFFvcomplex_alph(i:i+(nph-1),l);
                           
                    % Compute Coeffs
                    temp_p{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dPl));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dQl));
                    temp_v{ceil(i/nph),ceil(j/nph)} =  ((1./abs(Currents{ceil(i/nph),ceil(j/nph)})).*real(conj(Currents{ceil(i/nph),ceil(j/nph)}).*dIij_dVl));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dPl;
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dQl;
                    temp_v_complex{ceil(i/nph),ceil(j/nph)} =  dIij_dVl;
                else
                    temp_p{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_v{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_p_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_q_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                    temp_v_complex{ceil(i/nph),ceil(j/nph)} =  (zeros(nph,1));
                end
                
            end
        end
        
        IKp_abs{l} = cell2mat(temp_p);
        IKq_abs{l} = cell2mat(temp_q);
        IKv_abs{l} = cell2mat(temp_v);
        IKp{l} = cell2mat(temp_p_complex);
        IKq{l} = cell2mat(temp_p_complex);
        IKv{l} = cell2mat(temp_v_complex); 
    end
        
end
