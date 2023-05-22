function [K, Time] = SC_Voltage_V6_rectangular_unbalanced(S_star,E_star,idx1,idx3,Grid_para,Filter_para,idxCtrl,unblanced_3ph,filter)

Yac = Grid_para.Yac;
Ydc = Grid_para.Ydc;
Sac = S_star(1:Grid_para.n_ph*Grid_para.n_ac);
Eac = E_star(1:Grid_para.n_ph*Grid_para.n_ac);
Sdc = S_star(Grid_para.n_ph*Grid_para.n_ac+1:end);
Edc = E_star(Grid_para.n_ph*Grid_para.n_ac+1:end);
n_ph = Grid_para.n_ph;
Fl = Grid_para.pos_ac3;

% the function computes analytically the SC of the nodes when one or more
% VSC's are present that operate on QV control. The function takes the filter into account

% This function computes the nodal voltage sensitivity coefficients for
% each control variable of the nodes listed in idxCtrl
% INPUT
% - Y           nodal admittance matrix
% - S0          apparent power injections at nominal voltage for all buses
% - idx1ph      structure containing network indices of all nodes  
%    .slack    indices of the slack buses
%    .pq       indices of the PQ buses
%    .qv       indices of the QV buses
% - idx3ph      structure containing 3ph expanded indices of all nodes
%    .slack    indices of the slack buses
%    .pq       indices of the PQ buses
%    .qv       indices of the QV buses
% - idxCtrl     vector containing netowrk indices of the nodes for which
%               the user wants to compute the voltage sensitivity
%               coefficients (SCs)
% - nph         number of phases in considered grid
% - vdep        structure containing the voltage-dependent weights and
%               exponents
%    .alpha    weights active power voltage-dependency
%    .lambda    exponents active power voltage-dependency
%    .beta     weights reactive power voltage-dependency
%    .omega    exponents reactive power voltage-dependency
% - tol         tolerance for Newton-Raphson convergence criterion
% - n_max       maximum number of iterations
%
% OUTPUT
% - K           A cell structure with complex, magnitude and angle voltage SCs
%               the size is length(idxCtrl) x 4, where
%   - column 1, has the index of the node
%   - column 2 a 2x1 cell containing each the complex voltage SCs
%   pertaining to the control variables of the node. 
%   . If a node is PQ then the 1st cell contains SCs w.r.t. P_{0,l}^{\phi} 
%   (active power) and the 2nd cell contains SCs w.r.t. Q_{0,l}^{\phi} (reactive power). 
%   . If a node is QV then the 1st cell contains SCs w.r.t. P_{0,m}^{\phi} (active power) and the
%   2nd cell contains SCs w.r.t. |\bar{E}^{\phi}_m| (voltage magnitude). 
%   . If a node is slack then the 1st cell contains SCs w.r.t. |\bar{E}^{\phi}_k| (voltage magnitude) and the
%   2nd cell contains SCs w.r.t. \angle(\bar{E}^{\phi}_k) (phase-angle).
%   (!!!) Each cell contains a nph|N| x nph matrix where |N| is the number of nodes (1ph) and
%   each column corresponds to the SCs pertaining to a phase and a control variable as explained above.
%   - column 3&4 have the same structure as column 2 but contain,
%   respectively, the magnitude and phase-angle nodal voltage SCs.


% Y = YY_augmented3;
% S0 = S_star;
% E = E_star;




% Bcell = mat2cell(B, 3*ones(1,15), 3*ones(1,15));
% A = arrayfun(@(x) A_inv*cell2mat(x), Bcell, 'UniformOutput', false)
% A = cell2mat(A)


%% Construct A
% Timing
T = tic;

if unblanced_3ph
    alp = exp(2*pi/3*1i);
    A_inv = (1/3*[1 1     1; 
         1 alp   alp^2; 
         1 alp^2 alp]);
    b_inv = [0;1;0];
else
    A_inv = eye(n_ph);
    b_inv = ones(n_ph,1);
end

if n_ph == 1
    idx = idx1;
else
    idx = idx3;
end



% X_mat = repmat({[ones(3,1),zeros(3,1),zeros(3,1)]}, 1, 2*(size(Eac,1)/nph)-2);
% X_mat = blkdiag(X_mat{:},eye(4));


A_invl = repmat({A_inv}, 1, size(Eac,1)/n_ph);
A_invlm = blkdiag(A_invl{:}); %cell2mat(A_invl); %

b_inv4 = repmat(b_inv,size(Fl,1)/n_ph,1);

for i = 1:n_ph:size(Yac,1)
    for j = 1:n_ph:size(Yac,2)
        Yac_s(i:i+n_ph-1,j:j+n_ph-1) = A_inv.*(transpose(sum(Yac(i:i+n_ph-1,j:j+n_ph-1))) );
    end
end
    
% Create F_{in}^{\phi\phi'}, P0 and Q0
Fac = diag(Eac)*conj(Yac)*diag(conj(Eac));
Fac_s = (A_invlm.*Eac).*conj(Yac)*(conj(A_invlm.*Eac));% diag(A_invlm*Eac)*conj(Yac)*diag(conj(A_invlm*Eac)); %diag(A_invlm*Eac)*conj(A_invlm*(Yac)*diag(Eac));%
Fac_s = diag(Eac)*conj(Yac_s)*diag(conj(Eac));

for i = 1:n_ph:size(Fac,1)
    for j = 1:n_ph:size(Fac,2)
        Fac_s(i:i+n_ph-1,j:j+n_ph-1) = A_inv.*(transpose(sum(Fac(i:i+n_ph-1,j:j+n_ph-1))) );
    end
end

P0ac = real(Sac);
Q0ac = imag(Sac);

% Fdc = diag(Edc)*conj(Ydc(19:26,19:26))*diag(conj(Edc));
Fdc = real(diag(Edc)*conj(Ydc)*diag(conj(Edc)));
P0dc = real(Sdc);
Q0dc = imag(Sdc);
      
F = blkdiag(Fac,Fdc);
F_s = blkdiag(Fac_s,Fdc);
E = [Eac;Edc];
E_s = [Eac;Edc];


%% Filter
I_b =Grid_para.Y_b*Grid_para.V_b;
Imag = abs(Grid_para.YY*E_star);
R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b)./Grid_para.V_b./Imag;
R_eq = (R_eq_ctu*4/pi); % ???
R_eq(isnan(R_eq))=0;
R_eq(isinf(R_eq))=0;
Zf = (Filter_para.R + R_eq(sort([idx3.vscac_pq;idx3.vscac_vq]))) + 1i*Filter_para.X;

    
% influence of the filter

% alpha = ((conj(diag(Yac(Fl(:,1),Fl(:,2)))).*conj(Eac(Fl(:,2))) + conj(diag(Yac(Fl(:,1),Fl(:,1) ))).*conj(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*diag(Yac(Fl(:,1),Fl(:,1))).*Eac(Fl(:,1)) +  ...
%         ((diag(Yac(Fl(:,1),Fl(:,2)))).*(Eac(Fl(:,2))) + (diag(Yac(Fl(:,1),Fl(:,1) ))).*(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*conj(diag(Yac(Fl(:,1),Fl(:,1)))).*conj(Eac(Fl(:,1))) )./abs(Eac(Fl(:,1)));
% alpha = alpha(1:size(Fl,1));

% beta = ((conj(diag(Yac(Fl(:,1),Fl(:,2)))).*conj(Eac(Fl(:,2))) + conj(diag(Yac(Fl(:,1),Fl(:,1) ))).*conj(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*diag(Yac(Fl(:,1),Fl(:,1))).*Eac(Fl(:,1)) -  ...
%         ((diag(Yac(Fl(:,1),Fl(:,2)))).*(Eac(Fl(:,2))) + (diag(Yac(Fl(:,1),Fl(:,1) ))).*(Eac(Fl(:,1)))).*conj(Zf(1:3:end)).*conj(diag(Yac(Fl(:,1),Fl(:,1) ))).*conj(Eac(Fl(:,1))) );
% beta = beta(1:size(Fl,1));

if filter
    alpha = (conj(Yac(Fl(:,1),Fl(:,2)))*conj(Eac(Fl(:,2))) + conj(Yac(Fl(:,1),Fl(:,1)))*conj(Eac(Fl(:,1)))).*conj(Zf);
    beta = (Yac(Fl(:,1),Fl(:,2))*Eac(Fl(:,2)) + Yac(Fl(:,1),Fl(:,1))*Eac(Fl(:,1))).*conj(Zf);
    a = real(alpha.*(Yac(Fl(:,1),Fl(:,2))*Eac(Fl(:,2))) + beta.*(conj(Yac(Fl(:,1),Fl(:,2)))*conj(Eac(Fl(:,2)))));
    b = real(alpha.*(Yac(Fl(:,1),Fl(:,1))*Eac(Fl(:,1))) + beta.*(conj(Yac(Fl(:,1),Fl(:,1)))*conj(Eac(Fl(:,1)))));
    c = -imag(alpha.*(Yac(Fl(:,1),Fl(:,2))*Eac(Fl(:,2))) + beta.*(conj(Yac(Fl(:,1),Fl(:,2)))*conj(Eac(Fl(:,2)))));
    d = -imag(alpha.*(Yac(Fl(:,1),Fl(:,1))*Eac(Fl(:,1))) + beta.*(conj(Yac(Fl(:,1),Fl(:,1)))*conj(Eac(Fl(:,1)))));
    aa = zeros(length(idx.vscac_vq),length(idx.pqac)); aa(sub2ind(size(aa), (1:size(Fl,1))',Fl(:,2)-n_ph)) = a;
    bb = zeros(length(idx.vscac_vq),length(idx.vscac_vq)); bb(sub2ind(size(bb), (1:size(Fl,1))',(1:size(Fl,1))')) = b;
    cc = zeros(length(idx.vscac_vq),length(idx.pqac)); cc(sub2ind(size(cc), (1:size(Fl,1))',Fl(:,2)-n_ph)) = c;
    dd = zeros(length(idx.vscac_vq),length(idx.vscac_vq)); dd(sub2ind(size(dd), (1:size(Fl,1))',(1:size(Fl,1))')) = d;
else
    aa = 0;
    bb = 0;
    cc = 0;
    dd = 0;
end


% idx3ph.vscdc_vq = idx3ph.vscdc_vq-n_nodesac;
% idx3ph.pdc  =idx3ph.pdc-n_nodesac;
% Create each sub-matrix

dP_real = real(F)./repmat(abs(E).',[size(F,1) 1]) + diag(sum(real(F),2)./abs(E));
dP_imag = imag(F) - diag(sum(imag(F),2));

dQ_real = imag(F)./repmat(abs(E).',[size(F,1) 1]) + diag(sum(imag(F),2)./abs(E));
dQ_imag = -real(F) + diag(sum(real(F),2));

dV_real = -eye(size(F));% +  (dQ_real - diag(diag(dQ_real)));
dV_imag = dQ_imag;

dVdc = -eye(size(F));% +  (dP_real - diag(diag(dP_real)));


%% Balanced
F_sym = conj(Grid_para.YY).*(E_star);
H_sym = diag(Grid_para.YY*(E_star));

dP_real = real(F_sym) + real(H_sym);
dP_imag = imag(F_sym) + imag(H_sym);

dQ_real = -imag(F_sym) + imag(H_sym);
dQ_imag = real(F_sym) - real(H_sym);


dV_real = -eye(size(F));
dV_imag = eye(size(F));

%% Unbalanced

alp = exp(2*pi/3*1i);
T_inv = 1/3*   [1 1     1; 
                1 alp   alp^2; 
                1 alp^2 alp];
t_inv = [0;1;0];
TCell =  repmat({T_inv}, 1, Grid_para.n_ac);
tCell =  repmat({t_inv}, 1, Grid_para.n_ac);
ICell =  repmat({1}, 1, Grid_para.n_dc); 
Ttot = blkdiag(TCell{:},ICell{:});
ttot = vertcat(tCell{:},ICell{:});

% sequence voltages

dP_symE_real = real(Ttot);
dP_symE_imag = -imag(Ttot);

dQ_symE_real = imag(Ttot);
dQ_symE_imag = real(Ttot);



% sequence powers
F_sym = conj(Grid_para.YY).*(Ttot*E_star);
H_sym = diag(Grid_para.YY*(Ttot*E_star));

for i = 1:n_ph:size(Fac,1)
    for j = 1:n_ph:size(Fac,2)
        F_sym(i:i+n_ph-1,j:j+n_ph-1) = inv(T_inv).*(transpose(sum(F_sym(i:i+n_ph-1,j:j+n_ph-1))) );
        H_sym(i:i+n_ph-1,j:j+n_ph-1) = inv(T_inv).*(transpose(sum(H_sym(i:i+n_ph-1,j:j+n_ph-1))) );
    end  
end

dP_sym_real = real(F_sym) + real(H_sym);
dP_sym_imag = imag(F_sym) + imag(H_sym);

dQ_sym_real = -imag(F_sym) + imag(H_sym);
dQ_sym_imag = real(F_sym) - real(H_sym);

%Combine

indexes_0 = 1:3:Grid_para.n_ac*Grid_para.n_ph;
indexes_n = 3:3:Grid_para.n_ac*Grid_para.n_ph;

dP_sym_real(indexes_0,:) = dP_symE_real(indexes_0,:);
dP_sym_real(indexes_n,:) = dP_symE_real(indexes_n,:);

dP_sym_imag(indexes_0,:) = dP_symE_imag(indexes_0,:);
dP_sym_imag(indexes_n,:) = dP_symE_imag(indexes_n,:);

dQ_sym_real(indexes_0,:) = dQ_symE_real(indexes_0,:);
dQ_sym_real(indexes_n,:) = dQ_symE_real(indexes_n,:);

dQ_sym_imag(indexes_0,:) = dQ_symE_imag(indexes_0,:);
dQ_sym_imag(indexes_n,:) = dQ_symE_imag(indexes_n,:);



dV_sym_real = -diag(ttot);% +  (dQ_real - diag(diag(dQ_real)));
dV_sym_imag = diag(ttot);

dVdc_sym = -diag(ttot);% +  (dP_real - diag(diag(dP_real)));


dVdc_sym = dVdc;% +  (dP_real - diag(diag(dP_real)));

%% P - Q 
% A11   magE, angE (idx.pqac)
A11 = [dP_real(idx.pqac,idx.pqac) dP_imag(idx.pqac,idx.pqac);
       dQ_real(idx.pqac,idx.pqac) dQ_imag(idx.pqac,idx.pqac)];
% A12   Q   ,angE (idx.pvac)
A12 = [zeros(length(idx.pqac),length(idx.pvac)) dP_imag(idx.pqac,idx.pvac);
       zeros(length(idx.pqac),length(idx.pvac)) dQ_imag(idx.pqac,idx.pvac)];
% A13   Edc (idx.pdc)
A13 = [dP_real(idx.pqac,idx.pdc);
       dQ_real(idx.pqac,idx.pdc)];
% A14   Pdc (idx.vdc)
A14 = [zeros(length(idx.pqac),length(idx.vdc));
       zeros(length(idx.pqac),length(idx.vdc))];
% A15   magE, angE (vscac_pq)
A15 = [dP_real(idx.pqac,idx.vscac_pq) dP_imag(idx.pqac,idx.vscac_pq);
       dQ_real(idx.pqac,idx.vscac_pq) dQ_imag(idx.pqac,idx.vscac_pq)];
% A16   magE, angE (vscac_vq)
A16 = [dP_real(idx.pqac,idx.vscac_vq) dP_imag(idx.pqac,idx.vscac_vq);
       dQ_real(idx.pqac,idx.vscac_vq) dQ_imag(idx.pqac,idx.vscac_vq)];
% A17   Edc (vscdc_p)
A17 = [dP_real(idx.pqac,idx.vscdc_pq);
       dQ_real(idx.pqac,idx.vscdc_pq)];
   
%% P - V
% A21   magE, angE (idx.pqac)
A21 = [dP_real(idx.pvac,idx.pqac) dP_imag(idx.pvac,idx.pqac);
       dV_real(idx.pvac,idx.pqac) dV_imag(idx.pvac,idx.pqac)];
% A22   Q   ,angE (idx.pvac)
A22 = [zeros(length(idx.pvac),length(idx.pvac)) dP_imag(idx.pvac,idx.pvac);
       zeros(length(idx.pvac),length(idx.pvac)) dV_imag(idx.pvac,idx.pvac)];
% A23   Edc (idx.pdc)
A23 = [dP_real(idx.pvac,idx.pdc);
       dV_real(idx.pvac,idx.pdc)];
% A24   Pdc (idx.vdc)
A24 = [zeros(length(idx.pvac),length(idx.vdc));
       zeros(length(idx.pvac),length(idx.vdc))];
% A25   magE, angE (vscac_pq)
A25 = [dP_real(idx.pvac,idx.vscac_pq) dP_imag(idx.pvac,idx.vscac_pq);
       dV_real(idx.pvac,idx.vscac_pq) dV_imag(idx.pvac,idx.vscac_pq)];
% A26   magE, angE (vscac_vq)
A26 = [dP_real(idx.pvac,idx.vscac_vq) dP_imag(idx.pvac,idx.vscac_vq);
       dV_real(idx.pvac,idx.vscac_vq) dV_imag(idx.pvac,idx.vscac_vq)];
% A27   Edc (vscdc_p)
A27 = [dP_real(idx.pvac,idx.vscdc_pq);
       dV_real(idx.pvac,idx.vscdc_pq)];

%% Pdc 
% A11   magE, angE (idx.pqac)
A31 = [dP_real(idx.pdc,idx.pqac) dP_imag(idx.pdc,idx.pqac)];
% A12   Q   ,angE (idx.pvac)
A32 = [zeros(length(idx.pdc),length(idx.pvac)) dP_imag(idx.pdc,idx.pvac)];
% A13   Edc (idx.pdc)
A33 = [dP_real(idx.pdc,idx.pdc)];
% A14   Pdc (idx.vdc)
A34 = [zeros(length(idx.pdc),length(idx.vdc))];
% A15   magE, angE (vscac_pq)
A35 = [dP_real(idx.pdc,idx.vscac_pq) dP_imag(idx.pdc,idx.vscac_pq)];
% A16   magE, angE (vscac_vq)
A36 = [dP_real(idx.pdc,idx.vscac_vq) dP_imag(idx.pdc,idx.vscac_vq)];
% A17   Edc (vscdc_p)
A37 = [dP_real(idx.pdc,idx.vscdc_pq)];
   
%% Vdc 
% A11   magE, angE (idx.pqac)
A41 = [dP_real(idx.vdc,idx.pqac) dP_imag(idx.vdc,idx.pqac)];
% A12   Q   ,angE (idx.pvac)
A42 = [zeros(length(idx.vdc),length(idx.pvac)) dP_imag(idx.vdc,idx.pvac)];
% A13   Edc (idx.pdc)
A43 = [dP_real(idx.vdc,idx.pdc)];
% A14   Pdc (idx.vdc)
A44 = [dVdc(idx.vdc,idx.pdc)];
% A15   magE, angE (vscac_pq)
A45 = [dP_real(idx.vdc,idx.vscac_pq) dP_imag(idx.vdc,idx.vscac_pq)];
% A16   magE, angE (vscac_vq)
A46 = [dP_real(idx.vdc,idx.vscac_vq) dP_imag(idx.vdc,idx.vscac_vq)];
% A17   Edc (vscdc_p)
A47 = [dP_real(idx.vdc,idx.vscdc_pq)];

%% VSC PQ
% A51   magE, angE (idx.pqac)
A51 = [dP_sym_real(idx.vscac_pq,idx.pqac) dP_sym_imag(idx.vscac_pq,idx.pqac);
       dQ_sym_real(idx.vscac_pq,idx.pqac) dQ_sym_imag(idx.vscac_pq,idx.pqac)];
% A52   Q   ,angE (idx.pvac)
A52 = [zeros(length(idx.vscac_pq),length(idx.pvac)) dP_sym_imag(idx.vscac_pq,idx.pvac);
       zeros(length(idx.vscac_pq),length(idx.pvac)) dQ_sym_imag(idx.vscac_pq,idx.pvac)];
% A53   Edc (idx.pdc)
A53 = [dP_sym_real(idx.vscac_pq,idx.pdc);
       dQ_sym_real(idx.vscac_pq,idx.pdc)];
% A54   Pdc (idx.vdc)
A54 = [zeros(length(idx.vscac_pq),length(idx.vdc));
       zeros(length(idx.vscac_pq),length(idx.vdc))];
% A55   magE, angE (vscac_pq)
A55 = [dP_sym_real(idx.vscac_pq,idx.vscac_pq) dP_sym_imag(idx.vscac_pq,idx.vscac_pq);
       dQ_sym_real(idx.vscac_pq,idx.vscac_pq) dQ_sym_imag(idx.vscac_pq,idx.vscac_pq)];
% A56   magE, angE (vscac_vq)
A56 = [dP_sym_real(idx.vscac_pq,idx.vscac_vq) dP_sym_imag(idx.vscac_pq,idx.vscac_vq);
       dQ_sym_real(idx.vscac_pq,idx.vscac_vq) dQ_sym_imag(idx.vscac_pq,idx.vscac_vq)];
% A57   Edc (vscdc_p)
A57 = [dP_sym_real(idx.vscac_pq,idx.vscdc_pq);
       dQ_sym_real(idx.vscac_pq,idx.vscdc_pq)];
   
%% VSC VdcQ
% A61   magE, angE (idx.pqac)
A61 = [dP_sym_real(idx.vscac_vq,idx.pqac) dP_sym_imag(idx.vscac_vq,idx.pqac);
       dQ_sym_real(idx.vscac_vq,idx.pqac) dQ_sym_imag(idx.vscac_vq,idx.pqac)];
% A62   Q   ,angE (idx.pvac)
A62 = [zeros(length(idx.vscac_vq),length(idx.pvac)) dP_sym_imag(idx.vscac_vq,idx.pvac);
       zeros(length(idx.vscac_vq),length(idx.pvac)) dQ_sym_imag(idx.vscac_vq,idx.pvac)];
% A63   Edc (idx.pdc)

index_vsc = repmat(idx.vscdc_vq,1,n_ph)';
A63 = [dP_sym_real(index_vsc,idx.pdc); %!!!! maybe x and y have to be swapped
       dQ_sym_real(index_vsc,idx.pdc)];
% A64   Pdc (idx.vdc)
A64 = [zeros(length(idx.vscac_vq),length(idx.vdc));
       zeros(length(idx.vscac_vq),length(idx.vdc))];
% A65   magE, angE (vscac_pq)
A65 = [dP_sym_real(idx.vscac_vq,idx.vscac_pq) dP_sym_imag(idx.vscac_vq,idx.vscac_pq);
       dQ_sym_real(idx.vscac_vq,idx.vscac_pq) dQ_sym_imag(idx.vscac_vq,idx.vscac_pq)];
% A66   magE, angE (vscac_vq)
A66 = [dP_sym_real(idx.vscac_vq,idx.vscac_vq) dP_sym_imag(idx.vscac_vq,idx.vscac_vq);
       dQ_sym_real(idx.vscac_vq,idx.vscac_vq) dQ_sym_imag(idx.vscac_vq,idx.vscac_vq)];
% A67   Edc (vscdc_p)
A67 = [dP_sym_real(idx.vscac_vq,idx.vscdc_pq);
       dQ_sym_real(idx.vscac_vq,idx.vscdc_pq)];
   
%% Pdc vsc
% A71   magE, angE (idx.pqac)
A71 = [dP_real(idx.vscdc_pq,idx.pqac) dP_imag(idx.vscdc_pq,idx.pqac)];
% A72   Q   ,angE (idx.pvac)
A72 = [zeros(length(idx.vscdc_pq),length(idx.pvac)) dP_imag(idx.vscdc_pq,idx.pvac)];
% A73   Edc (idx.pdc)
A73 = [dP_real(idx.vscdc_pq,idx.pdc)];
% A74   Pdc (idx.vdc)
A74 = [zeros(length(idx.vscdc_pq),length(idx.vdc))];
% A75   magE, angE (vscac_pq)
A75 = [dP_real(idx.vscdc_pq,idx.vscac_pq) dP_imag(idx.vscdc_pq,idx.vscac_pq)];
% A76   magE, angE (vscac_vq)
A76 = [dP_real(idx.vscdc_pq,idx.vscac_vq) dP_imag(idx.vscdc_pq,idx.vscac_vq)];
% A77   Edc (vscdc_p)
A77 = [dP_real(idx.vscdc_pq,idx.vscdc_pq)];




%% Assemble A
A = [ A11 A12 A13 A14 A15 A16 A17; ...
      A21 A22 A23 A24 A25 A26 A27; ...
      A31 A32 A33 A34 A35 A36 A37; ...
      A41 A42 A43 A44 A45 A46 A47; ...
      A51 A52 A53 A54 A55 A56 A57; ...
      A61 A62 A63 A64 A65 A66 A67; ...
      A71 A72 A73 A74 A75 A76 A77];

% Anew=[A11 A16 A13;...
%       A61 A66 A63;...
%       A31 A36 A33];
  
Time.A = toc(T);

%%
% 
% %% eq 3
% 
% sc = complex(J_VR(:,1),J_VX(:,1))
% sc_real = 1./abs(E(2:end)) .* real(conj(E(2:end)) .* sc)
% sc_imag = 1./abs(E(2:end)).^2 .* imag(conj(E(2:end)) .* sc)
% 
% sum(sum(real(F(15,:)))/abs(E(15)) + real(F(15,15))/abs(E(15)))*sc_real(15-1)+...  %15
% sum(real(F(15,[1:14,16]))./transpose(abs(E([1:14,16]))))*sc_real(7-1)+...  %1:14,16
% sum(-sum(imag(F(15,:))) - imag(F(15,15)))*sc_imag(15-1)+...  %15
% sum(imag(F(15,[1:14,16])))*sc_imag(7-1)+... %1:14,16
% sum(-real(F(17,[18:20]))./transpose(abs(E([18:20]))))*sc_real(19-1) %18:20
% 
% sum(sum(real(F(17,:)))/abs(E(17)) + real(F(17,17))/abs(E(17)))*sc_real(17-1) %17
% 
% 
% %validated -> eq3 is correct
% 
% sum(sum(real(F(15,:)))/abs(E(15)) + real(F(15,15))/abs(E(15)))*J_VR(15-1,1) +... %15
% sum(real(F(15,[1:14,16]))./transpose(abs(E([1:14,16]))))*J_VR(7-1,1) +...%1:14,16
% sum(-sum(imag(F(15,:))) - imag(F(15,15)))*J_VX(15-1,1) +...%15
% sum(imag(F(15,[1:14,16])))*J_VX(7-1,1) %1:14,16
% 
% sum(sum(real(F(17,:)))/abs(E(17)) + real(F(17,17))/abs(E(17)))*J_VR(17-1,1) +...%17
% sum(-real(F(17,[18:20]))./transpose(abs(E([18:20]))))*J_VR(19-1,1) %18:20
% 
% 
% sum(real(F(16,:)))/abs(E(16)) + real(F(16,16))/abs(E(16)) %16
% real(F(16,[1:14,15]))./transpose(abs(E([1:14,15]))) %1:14,15
% 
% sum(real(F(18,:)))/abs(E(18)) + real(F(18,18))/abs(E(18)) %18
% -real(F(18,[17,19:20]))./transpose(abs(E([17,19:20]))) %17,19:20
% 
% -sum(imag(F(16,:))) - imag(F(16,16)) %16
% imag(F(16,[1:14,15])) %1:14,15
% 
% %% eq 4
% 
% sc = complex(J_QR(:,14),J_QX(:,14))
% sc_real = 1./abs(E(2:end)) .* real(conj(E(2:end)) .* sc)
% sc_imag = 1./abs(E(2:end)).^2 .* imag(conj(E(2:end)) .* sc)
% 
% (sum(imag(F(15,:)))/abs(E(15)) + imag(F(15,15))/abs(E(15)))*sc_real(15-1)+... %15
% sum(imag(F(15,[1:14,16]))./transpose(abs(E([1:14,16]))))*sc_real(7-1)+... %1:14,16
% (sum(real(F(15,:))) - real(F(15,15)))*sc_imag(15-1)+... %15
% sum(-real(F(15,[1:14,16])))*sc_imag(7-1) %1:14,16
% 
% %validated -> sum = -1
% A4 = [ A41 A42 A43 A44 A45]
% %%
% 
% %% eq 5
% 
% sc = complex(J_PR(:,14),J_PX(:,14))
% sc_real = 1./abs(E(2:end)) .* real(conj(E(2:end)) .* sc)
% sc_imag = 1./abs(E(2:end)).^2 .* imag(conj(E(2:end)) .* sc)
% 
% (sum(real(F(19,:)))/abs(E(19)) + real(F(19,19))/abs(E(19)))*sc_real(19-1)+... %15
% sum(real(F(19,[17,18,20]))./transpose(abs(E([17,18,20])))).*sc_real([17,18,20]-1) %1:14,16
% 
% %validated -> sum = 1
% A5 = [ A51 A52 A53 A54 A55]



 %% Compute nodal voltage SCs for each control variable

% Initialize Outputs
K = cell(length(idxCtrl),4);
Time.K = [];

for id_x = 1:length(idxCtrl)
    % Initialize cell entry
    K{id_x,1} = idxCtrl(id_x); % 1ph index
    K{id_x,2} = cell(2,1);     % Nodal Voltage SCs (complex)
    K{id_x,3} = cell(2,1);     % Magnitude Voltage SCs
    K{id_x,4} = cell(2,1);     % Angle Voltage SCs
    
    % Identify node-type of idxCtrl [1: slack, 2: pq, 3: qv]
     node_type = 1*sum( idxCtrl(id_x) == idx.slack ) + ...
                 2*sum( idxCtrl(id_x) == idx.pqac ) + ...
                 3*sum( idxCtrl(id_x) == idx.pvac ) + ...
                 4*sum( idxCtrl(id_x) == idx.pdc ) + ...
                 5*sum( idxCtrl(id_x) == idx.vdc ) + ...
                 6*sum( idxCtrl(id_x) == idx.vscac_pq ) + ...
                 7*sum( idxCtrl(id_x) == idx.vscac_vq ) + ...
                 8*sum( idxCtrl(id_x) == idx.vscdc_pq ) + ...
                 9*sum( idxCtrl(id_x) == idx.vscdc_vq ); 
    for ctrl_var = 1:2
       
     T = tic;   
           
     % Initialize outputs
     K{id_x,2}{ctrl_var,1} = zeros(size(E,1),1); % complex
     K{id_x,3}{ctrl_var,1} = zeros(size(E,1),1); % magnitude
     K{id_x,4}{ctrl_var,1} = zeros(size(E,1),1); % angle

     % Get 3ph index
     tmp_idx = idxCtrl(id_x);
     
     %% V controlled nodes
     switch( node_type )
         case 1 % Slack node
             if (ctrl_var == 1) 
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % magnitude
             elseif (ctrl_var == 2) 
                K{id_x,4}{ctrl_var,1}(tmp_idx,1) = 1; % angle
             end
             
         case 3 % PV node
             if(ctrl_var == 1) % P
                % Do nothing as derivative is zero 
             elseif (ctrl_var == 2) % magnitude V
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % magnitude
             end
             
         case 5 % Vdc node
             if (ctrl_var == 1)
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % magnitude
             elseif (ctrl_var == 1)
                 % Do nothing as derivative is zero
             end
             
         case 9 % VSCdc_vq node
             if (ctrl_var == 1) 
                 K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % magnitude
             elseif (ctrl_var == 1) 
                 % Do nothing as derivative is zero
             end
             
         otherwise 
             % Do nothing as derivatives are all zeros.
     end
     
     %% Construct u(X)
      switch( node_type )
         case 1 % slack node
             if (ctrl_var == 1) 
                u1 = [- dP_real(idx.pqac,tmp_idx);
                      - dQ_real(idx.pqac,tmp_idx)];
                u2 = [- dP_real(idx.pvac,tmp_idx);
                      - dV_real(idx.pvac,tmp_idx)];
                u3 = [- dP_real(idx.pdc,tmp_idx)];
                u4 = [- dVdc(idx.vdc,tmp_idx)];
                u5 = [- dP_sym_real(idx.vscac_pq,tmp_idx);
                      - dQ_sym_real(idx.vscac_pq,tmp_idx)];
                u6 = [- dV_sym_real(idx.vscac_vq,tmp_idx);
                      - dQ_sym_real(idx.vscac_vq,tmp_idx)];
                u7 = [- dP_real(idx.vscdc_pq,tmp_idx)];
             elseif (ctrl_var == 2) 
                u1 = [- dP_imag(idx.pqac,tmp_idx);
                      - dQ_imag(idx.pqac,tmp_idx)];
                u2 = [- dP_imag(idx.pvac,tmp_idx);
                      - dV_imag(idx.pvac,tmp_idx)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [- dP_sym_imag(idx.vscac_pq,tmp_idx);
                      - dQ_sym_imag(idx.vscac_pq,tmp_idx)];
                u6 = [- dV_sym_imag(idx.vscac_vq,tmp_idx);
                      - dQ_sym_imag(idx.vscac_vq,tmp_idx)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             end
             
          case 2 % pqac node
             if (ctrl_var == 1) 
                u1 = [-dV_real(idx.pqac,tmp_idx);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             elseif (ctrl_var == 2) 
                u1 = [zeros(length(idx.pqac),1);                      
                      -dV_imag(idx.pqac,tmp_idx)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             end
             
          case 3 % pvac node
             if (ctrl_var == 1) % if X = P0^{\phi}_{0,l}
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [-dV_real(idx.pvac,tmp_idx);
                      zeros(length(idx.pqac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];      
             elseif (ctrl_var == 2)
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      -dP_imag(idx.pvac,tmp_idx)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];     
             end  
             
          case 4 % pdc mode
              if (ctrl_var == 1) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [-dV_real(idx.pdc,tmp_idx)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             elseif (ctrl_var == 2) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
              end
              
          case 5 % vdc mode
              if (ctrl_var == 1) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [-dP_real(idx.vdc,tmp_idx)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             elseif (ctrl_var == 2) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
              end
              
          case 6 % vscac_pq node
             if (ctrl_var == 1) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [-dV_sym_real(idx.vscac_pq,tmp_idx);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             elseif (ctrl_var == 2) 
                u1 = [zeros(length(idx.pqac),1);                      
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      -dV_sym_imag(idx.vscac_pq,tmp_idx)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             end
             
          case 7 % vscac_vq node
              if (ctrl_var == 1) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             elseif (ctrl_var == 2) 
                u1 = [zeros(length(idx.pqac),1);                      
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)]; 
                u6 = [zeros(length(idx.vscac_vq),1);
                      -dV_sym_imag(idx.vscac_vq,tmp_idx)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
              end
             
          case 8 % vscdc_pq mode
             if (ctrl_var == 1) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [-dV_real(idx.vscdc_pq,tmp_idx)];
             elseif (ctrl_var == 2) 
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             end
              
          case 9 % vscdc_vq node
              if (ctrl_var == 1) 
                index_vsc_v = repmat(idx.vscdc_vq,1,n_ph)';
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [-dP_real(idx.pdc,tmp_idx)]; %might be off
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [-dP_sym_real(index_vsc_v,tmp_idx);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
             elseif (ctrl_var == 2) 
                u1 = [zeros(length(idx.pqac),1);                      
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [zeros(length(idx.pdc),1)];
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [zeros(length(idx.vscac_vq),1);
                      zeros(length(idx.vscac_vq),1)];
                u7 = [zeros(length(idx.vscdc_pq),1)];
              end
              
          otherwise % "error"
              warning('something is wrong, mate!')
      end    
     

     
     %% Solve A*x(X)=u(X) 
     u = [u1;u2;u3;u4;u5;u6;u7];
     x = linsolve(A,u);
   
     
     %% Assemble magnitude nodal voltage SCs
     % K{id_x,3}{ctrl_var,1}(idx.slack,1)       %Already good (not taken into account)
     K{id_x,3}{ctrl_var,1}(idx.pqac,1) =     x( 1: length(idx.pqac) ); 
     %K{id_x,3}{ctrl_var,1}(idx.pvac,1) =       %Already good; 
     K{id_x,3}{ctrl_var,1}(idx.pdc,1) =      x( 2*length(idx.pqac) + 2*length(idx.pvac) + 1: 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) ); 
     %K{id_x,3}{ctrl_var,1}(idx.vdc,1) =        %Already good
     K{id_x,3}{ctrl_var,1}(idx.vscac_pq,1) = x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 1: 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + length(idx.vscac_pq) );
     K{id_x,3}{ctrl_var,1}(idx.vscac_vq,1) = x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 1 : 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + length(idx.vscac_vq) );
     K{id_x,3}{ctrl_var,1}(idx.vscdc_pq,1) = x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 2*length(idx.vscac_vq) + 1 : 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 2*length(idx.vscac_vq) + length(idx.vscdc_pq) );
     % K{id_x,3}{ctrl_var,1}(idx.vscdc_pq,1) =  %Already good (not taken into account)
     
     
     %% Assemble phase-angle nodal voltage SCs
     % K{id_x,4}{ctrl_var,1}(idx.slack,1)       %Already good
     K{id_x,4}{ctrl_var,1}(idx.pqac,1) =     x( length(idx.pqac) + 1: 2*length(idx.pqac) ); 
     K{id_x,4}{ctrl_var,1}(idx.pvac,1) =     x( 2*length(idx.pqac) + length(idx.pvac) +1 : 2*length(idx.pqac) + 2*length(idx.pvac) ); 
     %K{id_x,4}{ctrl_var,1}(idx.pdc,1) =        %DC quantity, angle does not exists 
     %K{id_x,4}{ctrl_var,1}(idx.vdc,1) =        %DC quantity, angle does not exists
     K{id_x,4}{ctrl_var,1}(idx.vscac_pq,1) = x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + length(idx.vscac_pq) + 1: 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) );
     K{id_x,4}{ctrl_var,1}(idx.vscac_vq,1) = x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + length(idx.vscac_vq) + 1 : 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 2*length(idx.vscac_vq) );
     % K{id_x,4}{ctrl_var,1}(idx.vscdc_pq,1) =  %DC quantity, angle does not exists 
     % K{id_x,4}{ctrl_var,1}(idx.vscdc_pq,1) =  %DC quantity, angle does not exists

     % Infer complex nodal voltage SCs 
     
     K{id_x,2}{ctrl_var,1}(:,1) = complex(K{id_x,3}{ctrl_var,1}(:,1) , K{id_x,4}{ctrl_var,1}(:,1)); %
%      K{id_x,2}{ctrl_var,1}(:,1) = E.*( (1./abs(E)).*K{id_x,3}{ctrl_var,1}(:,1) + 1i*K{id_x,4}{ctrl_var,1}(:,1) );
     
     % Time
     Time.K = [Time.K; toc(T)];
     
    end
 end

end