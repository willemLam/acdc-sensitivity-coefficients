function [Zloss,Zfilter] = get_filter_loss_impedance(Grid_para,Filter_para,E_star)

    n_dc = Grid_para.n_dc;
    n_ac = Grid_para.n_ac;
    V_b = Grid_para.V_b;
    A_b = Grid_para.A_b;
    Y_b = Grid_para.Y_b;
    
    %% RL FILTER
    R = Filter_para.R;
    X = Filter_para.X;
    
    %% IGBT losses
    n_AFE = Grid_para.n_AFE;
    pos_ac3 = Grid_para.pos_ac3;
    G = Grid_para.G;
    B = Grid_para.B;
    Y = complex(G,B);
    I_b=A_b/(V_b);
    alp = exp(2*pi/3*1i);
    A = 1/3*[1 1     1; 
             1 alp   alp^2; 
             1 alp^2 alp];

    ACell =  repmat({A}, 1, n_ac);
    ICell =  repmat({1}, 1, n_dc); 
    Atot = blkdiag(ACell{:},ICell{:});

    Imag = abs(Y(pos_ac3(:,1),:)*Atot*E_star);
    Imag = Imag(2:3:end);
%     R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b)./Grid_para.V_b./Imag;
    R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*sqrt(2)*I_b)./Grid_para.V_b./(Imag/sqrt(3)*sqrt(2));
  
    R_eq = (R_eq_ctu*4/pi); 
    R_eq(isnan(R_eq))=0;
    R_eq(isinf(R_eq))=0;


    Zloss = R_eq.*Filter_para.Include_losses;
    Zfilter = repmat(complex(R,X),n_AFE,1);
    
%     I_star = Y(pos_ac3(:,1),:)*Atot*E_star;
%     I_star = I_star(2:3:end);
%     E_star(pos_ac3(1:3:end,1),:) .* conj(Zloss.*I_star)
    
end
