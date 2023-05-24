function [K, Time] = SC_voltage_rectangular(E,idx,Grid_para,idxCtrl)

%% For debugging purposes
% E = E_star_augmented;
% idx = idx3_augmented;
% Grid_para = Grid_para_augmented;


% the function computes analytically the SC of the nodes when one or more
% VSC's are present that operate on QV control. The function takes the filter into account

% INPUT
% - E           Nodal voltage phasors (including virtual node) 
% - idx         Structure containing network indices of all nodes
% - Grid_para   Structure of all networks parameters (e.g. nodal admittance matrix, number of nodes...)
% - S0          Epparent power injections at nominal voltage for all buses
% - idxCtrl     Index of the nodes (3-ph) of the nodes for which the SC have to be computed 
%
% OUTPUT
% - K           A cell structure with complex, real, imaginary, magnitude and angle voltage SCs
%               the size is length(idxCtrl) x 4, where
%   - column 1, has the index of the node
%   - column 2 a 2x1 cell containing each the complex voltage SCs
%   - column 3&4&5&6 have the same structure as column 2 but contain,
%     respectively, the real, imaginary, magnitude and phase-angle nodal voltage SCs.


%% Construct A
% Timing
T = tic;
n_ph = Grid_para.n_ph;
n_ac = Grid_para.n_ac;

%% Balanced
F = conj(Grid_para.YY).*E;
H = diag(Grid_para.YY*E);

dP_real = real(F) + real(H);
dP_imag = imag(F) + imag(H);

dQ_real = -imag(F) + imag(H);
dQ_imag = real(F) - real(H);

dV_real = -eye(size(F));% +  (dQ_real - diag(diag(dQ_real)));
dV_imag = eye(size(F));

dVdc = -eye(size(F));% +  (dP_real - diag(diag(dP_real)));

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
F_sym = conj(Grid_para.YY).*(Ttot*E);
H_sym = diag(Grid_para.YY*(Ttot*E));

for i = 1:n_ph:n_ac*n_ph
    for j = 1:n_ph:n_ac*n_ph
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
M_zero = zeros(size(dP_sym_real));

dP_symp_real = dP_sym_real;
dP_symp_imag = dP_sym_imag;
dQ_symp_real = dQ_sym_real;
dQ_symp_imag = dQ_sym_imag;

dP_symp_real(indexes_0,:) = M_zero(indexes_0,:);
dP_symp_real(indexes_n,:) = M_zero(indexes_n,:);
dP_symp_imag(indexes_0,:) = M_zero(indexes_0,:);
dP_symp_imag(indexes_n,:) = M_zero(indexes_n,:);
dQ_symp_real(indexes_0,:) = M_zero(indexes_0,:);
dQ_symp_real(indexes_n,:) = M_zero(indexes_n,:);
dQ_symp_imag(indexes_0,:) = M_zero(indexes_0,:);
dQ_symp_imag(indexes_n,:) = M_zero(indexes_n,:);

dP_sym_real(indexes_0,:) = dP_symE_real(indexes_0,:);
dP_sym_real(indexes_n,:) = dP_symE_real(indexes_n,:);
dP_sym_imag(indexes_0,:) = dP_symE_imag(indexes_0,:);
dP_sym_imag(indexes_n,:) = dP_symE_imag(indexes_n,:);
dQ_sym_real(indexes_0,:) = dQ_symE_real(indexes_0,:);
dQ_sym_real(indexes_n,:) = dQ_symE_real(indexes_n,:);
dQ_sym_imag(indexes_0,:) = dQ_symE_imag(indexes_0,:);
dQ_sym_imag(indexes_n,:) = dQ_symE_imag(indexes_n,:);

dV_sym_real = -diag(ttot);
dV_sym_imag = diag(ttot);

dVdc_sym = -diag(ttot);
dVdc_sym = dVdc;

%% P - Q nodes
% A11   Ereal, Eimag (idx.pqac)
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
% A15   Ereal, Eimag (vscac_pq)
A15 = [dP_real(idx.pqac,idx.vscac_pq) dP_imag(idx.pqac,idx.vscac_pq);
       dQ_real(idx.pqac,idx.vscac_pq) dQ_imag(idx.pqac,idx.vscac_pq)];
% A16   Ereal, Eimag (vscac_vq_v
A16 = [dP_real(idx.pqac,idx.vscac_vq) dP_imag(idx.pqac,idx.vscac_vq);
       dQ_real(idx.pqac,idx.vscac_vq) dQ_imag(idx.pqac,idx.vscac_vq)];
% A17   Edc (vscdc_p)
A17 = [dP_real(idx.pqac,idx.vscdc_pq);
       dQ_real(idx.pqac,idx.vscdc_pq)];
   
%% P - V
% A21   Ereal, Eimag (idx.pqac)
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
% A25   Ereal, Eimag (vscac_pq)
A25 = [dP_real(idx.pvac,idx.vscac_pq) dP_imag(idx.pvac,idx.vscac_pq);
       dV_real(idx.pvac,idx.vscac_pq) dV_imag(idx.pvac,idx.vscac_pq)];
% A26   Ereal, Eimag (vscac_vq)
A26 = [dP_real(idx.pvac,idx.vscac_vq) dP_imag(idx.pvac,idx.vscac_vq);
       dV_real(idx.pvac,idx.vscac_vq) dV_imag(idx.pvac,idx.vscac_vq)];
% A27   Edc (vscdc_p)
A27 = [dP_real(idx.pvac,idx.vscdc_pq);
       dV_real(idx.pvac,idx.vscdc_pq)];

%% Pdc 
% A11   Ereal, Eimag (idx.pqac)
A31 = [dP_real(idx.pdc,idx.pqac) dP_imag(idx.pdc,idx.pqac)];
% A12   Q   ,angE (idx.pvac)
A32 = [zeros(length(idx.pdc),length(idx.pvac)) dP_imag(idx.pdc,idx.pvac)];
% A13   Edc (idx.pdc)
A33 = [dP_real(idx.pdc,idx.pdc)];
% A14   Pdc (idx.vdc)
A34 = [zeros(length(idx.pdc),length(idx.vdc))];
% A15   Ereal, Eimag (vscac_pq)
A35 = [dP_real(idx.pdc,idx.vscac_pq) dP_imag(idx.pdc,idx.vscac_pq)];
% A16   Ereal, Eimag (vscac_vq)
A36 = [dP_real(idx.pdc,idx.vscac_vq) dP_imag(idx.pdc,idx.vscac_vq)];
% A17   Edc (vscdc_p)
A37 = [dP_real(idx.pdc,idx.vscdc_pq)];
   
%% Vdc 
% A11   Ereal, Eimag (idx.pqac)
A41 = [dP_real(idx.vdc,idx.pqac) dP_imag(idx.vdc,idx.pqac)];
% A12   Q   ,angE (idx.pvac)
A42 = [zeros(length(idx.vdc),length(idx.pvac)) dP_imag(idx.vdc,idx.pvac)];
% A13   Edc (idx.pdc)
A43 = [dP_real(idx.vdc,idx.pdc)];
% A14   Pdc (idx.vdc)
A44 = [dVdc(idx.vdc,idx.vdc)];
% A15   Ereal, Eimag (vscac_pq)
A45 = [dP_real(idx.vdc,idx.vscac_pq) dP_imag(idx.vdc,idx.vscac_pq)];
% A16   Ereal, Eimag (vscac_vq)
A46 = [dP_real(idx.vdc,idx.vscac_vq) dP_imag(idx.vdc,idx.vscac_vq)];
% A17   Edc (vscdc_p)
A47 = [dP_real(idx.vdc,idx.vscdc_pq)];

%% VSC PQ
% A51   Ereal, Eimag (idx.pqac)
A51 = [dP_symp_real(idx.vscac_pq,idx.pqac) dP_symp_imag(idx.vscac_pq,idx.pqac);
       dQ_symp_real(idx.vscac_pq,idx.pqac) dQ_symp_imag(idx.vscac_pq,idx.pqac)];
% A52   Q   ,angE (idx.pvac)
A52 = [zeros(length(idx.vscac_pq),length(idx.pvac)) dP_symp_imag(idx.vscac_pq,idx.pvac);
       zeros(length(idx.vscac_pq),length(idx.pvac)) dQ_symp_imag(idx.vscac_pq,idx.pvac)];
% A53   Edc (idx.pdc)
A53 = [dP_symp_real(idx.vscac_pq,idx.pdc);
       dQ_symp_real(idx.vscac_pq,idx.pdc)];
% A54   Pdc (idx.vdc)
A54 = [zeros(length(idx.vscac_pq),length(idx.vdc));
       zeros(length(idx.vscac_pq),length(idx.vdc))];
% A55   Ereal, Eimag (vscac_pq)
A55 = [dP_sym_real(idx.vscac_pq,idx.vscac_pq) dP_sym_imag(idx.vscac_pq,idx.vscac_pq);
       dQ_sym_real(idx.vscac_pq,idx.vscac_pq) dQ_sym_imag(idx.vscac_pq,idx.vscac_pq)];
% A56   Ereal, Eimag (vscac_vq)
A56 = [dP_symp_real(idx.vscac_pq,idx.vscac_vq) dP_symp_imag(idx.vscac_pq,idx.vscac_vq);
       dQ_symp_real(idx.vscac_pq,idx.vscac_vq) dQ_symp_imag(idx.vscac_pq,idx.vscac_vq)];
% A57   Edc (vscdc_p)
A57 = [dP_symp_real(idx.vscac_pq,idx.vscdc_pq);
       dQ_symp_real(idx.vscac_pq,idx.vscdc_pq)];
   
%% VSC VdcQ_V
% A61   Ereal, Eimag (idx.pqac)
A61 = [dP_symp_real(idx.vscac_vq,idx.pqac) dP_symp_imag(idx.vscac_vq,idx.pqac);
        dQ_symp_real(idx.vscac_vq,idx.pqac) dQ_symp_imag(idx.vscac_vq,idx.pqac)];
% A62   Q   ,angE (idx.pvac)
A62 = [zeros(length(idx.vscac_vq),length(idx.pvac)) dP_symp_imag(idx.vscac_vq,idx.pvac);
       zeros(length(idx.vscac_vq),length(idx.pvac)) dQ_symp_imag(idx.vscac_vq,idx.pvac)];
% A63   Edc (idx.pdc)
index_vsc = repmat(idx.vscdc_vq,1,n_ph)';
A63 = [dP_symp_real(index_vsc,idx.pdc).*repmat([0;1;0],length(idx.vscdc_vq),1); %!!!! maybe x and y have to be swapped
        dQ_symp_real(idx.vscac_vq,idx.pdc).*repmat([0;1;0],length(idx.vscdc_vq),1)]; % ??
% A64   Pdc (idx.vdc)
A64 = [zeros(length(idx.vscac_vq),length(idx.vdc));
       zeros(length(idx.vscac_vq),length(idx.vdc))];
% A65   Ereal, Eimag (vscac_pq)
A65 = [dP_symp_real(idx.vscac_vq,idx.vscac_pq) dP_symp_imag(idx.vscac_vq,idx.vscac_pq);
        dQ_symp_real(idx.vscac_vq,idx.vscac_pq) dQ_symp_imag(idx.vscac_vq,idx.vscac_pq)];
% A66   Ereal, Eimag (vscac_vq)
A66 = [dP_sym_real(idx.vscac_vq,idx.vscac_vq) dP_sym_imag(idx.vscac_vq,idx.vscac_vq);
         dQ_sym_real(idx.vscac_vq,idx.vscac_vq) dQ_sym_imag(idx.vscac_vq,idx.vscac_vq)];
% A67   Edc (vscdc_p)
A67 = [dP_symp_real(idx.vscac_vq,idx.vscdc_pq);
        dQ_symp_real(idx.vscac_vq,idx.vscdc_pq)];
   
%% Pdc vsc
% A71   Ereal, Eimag (idx.pqac)
A71 = [dP_real(idx.vscdc_pq,idx.pqac) dP_imag(idx.vscdc_pq,idx.pqac)];
% A72   Q   ,angE (idx.pvac)
A72 = [zeros(length(idx.vscdc_pq),length(idx.pvac)) dP_imag(idx.vscdc_pq,idx.pvac)];
% A73   Edc (idx.pdc)
A73 = [dP_real(idx.vscdc_pq,idx.pdc)];
% A74   Pdc (idx.vdc)
A74 = [zeros(length(idx.vscdc_pq),length(idx.vdc))];
% A75   Ereal, Eimag (vscac_pq)
A75 = [dP_real(idx.vscdc_pq,idx.vscac_pq) dP_imag(idx.vscdc_pq,idx.vscac_pq)];
% A76   Ereal, Eimag (vscac_vq)
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

Time.A = toc(T);


 %% Compute nodal voltage SCs for each control variable

% Initialize Outputs
K = cell(length(idxCtrl),4);
Time.K = [];

for id_x = 1:length(idxCtrl)
    % Initialize cell entry
    K{id_x,1} = idxCtrl(id_x); % 3ph index
    K{id_x,2} = cell(2,1);     % Nodal Voltage SCs (complex)
    K{id_x,3} = cell(2,1);     % Real Voltage SCs
    K{id_x,4} = cell(2,1);     % Imaginary Voltage SCs
    K{id_x,5} = cell(2,1);     % Magnitude Voltage SCs
    K{id_x,6} = cell(2,1);     % Angle Voltage SCs
    
    % Identify node-type of idxCtrl
     node_type = 1*sum( idxCtrl(id_x) == idx.slack ) + ...
                 2*sum( idxCtrl(id_x) == idx.pqac) + ...
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
     K{id_x,3}{ctrl_var,1} = zeros(size(E,1),1); % real
     K{id_x,4}{ctrl_var,1} = zeros(size(E,1),1); % imag
     K{id_x,5}{ctrl_var,1} = zeros(size(E,1),1); % magnitude
     K{id_x,6}{ctrl_var,1} = zeros(size(E,1),1); % angle

     % Get 3ph index
     tmp_idx = idxCtrl(id_x);
     
     %% V controlled nodes
     switch( node_type )
         case 1 % Slack node
             if (ctrl_var == 1) 
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % real
             elseif (ctrl_var == 2) 
                K{id_x,4}{ctrl_var,1}(tmp_idx,1) = 0; % imag
             end
             
         case 3 % PV node
             if(ctrl_var == 1) % P
                % Do nothing as derivative is zero 
             elseif (ctrl_var == 2) % magnitude V
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % magnitude
             end
             warning('nope')
             
         case 5 % Vdc node
             if (ctrl_var == 1)
                K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % real
             elseif (ctrl_var == 1)
                 % Do nothing as derivative is zero
             end
             
         case 9 % VSCdc_vq node
             if (ctrl_var == 1) 
                 K{id_x,3}{ctrl_var,1}(tmp_idx,1) = 1; % real
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
                      - dV_imag(idx.pvac,tmp_idx)];
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
                u6 = [- dV_sym_real(idx.vscac_vq,tmp_idx);
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
                      -dP_real(idx.pvac,tmp_idx)];
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
                u3 = [-dVdc(idx.pdc,tmp_idx)];
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

           case 7 % vscac_vq_v node
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
                index_vsc = repmat(idx.vscdc_vq,1,n_ph)';
                u1 = [zeros(length(idx.pqac),1);
                      zeros(length(idx.pqac),1)];
                u2 = [zeros(length(idx.pvac),1);
                      zeros(length(idx.pvac),1)];
                u3 = [-dP_real(idx.pdc,tmp_idx)]; %might be off
                u4 = [zeros(length(idx.vdc),1)];
                u5 = [zeros(length(idx.vscac_pq),1);
                      zeros(length(idx.vscac_pq),1)];
                u6 = [-dP_sym_real(index_vsc,tmp_idx).*repmat([0;1;0],length(idx.vscdc_vq),1);
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
     K{id_x,3}{ctrl_var,1}(idx.pqac,1) =       x( 1: length(idx.pqac) ); 
     %K{id_x,3}{ctrl_var,1}(idx.pvac,1) =       %Already good; 
     K{id_x,3}{ctrl_var,1}(idx.pdc,1) =        x( 2*length(idx.pqac) + 2*length(idx.pvac) + 1: 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) ); 
     %K{id_x,3}{ctrl_var,1}(idx.vdc,1) =        %Already good
     K{id_x,3}{ctrl_var,1}(idx.vscac_pq,1) =   x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 1: 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + length(idx.vscac_pq) );
     K{id_x,3}{ctrl_var,1}(idx.vscac_vq,1) =   x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 1 : 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + length(idx.vscac_vq) );
     K{id_x,3}{ctrl_var,1}(idx.vscdc_pq,1) =   x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 2*length(idx.vscac_vq) + 1 : 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 2*length(idx.vscac_vq) + length(idx.vscdc_pq) );
     % K{id_x,3}{ctrl_var,1}(idx.vscdc_pq,1) =  %Already good (not taken into account)
     
     %% Assemble phase-angle nodal voltage SCs
     % K{id_x,4}{ctrl_var,1}(idx.slack,1)       %Already good
     K{id_x,4}{ctrl_var,1}(idx.pqac,1) =     x( length(idx.pqac) + 1: 2*length(idx.pqac) ); 
     K{id_x,4}{ctrl_var,1}(idx.pvac,1) =     x( 2*length(idx.pqac) + length(idx.pvac) +1 : 2*length(idx.pqac) + 2*length(idx.pvac) ); 
     %K{id_x,4}{ctrl_var,1}(idx.pdc,1) =        %DC quantity, imag part is zero
     %K{id_x,4}{ctrl_var,1}(idx.vdc,1) =        %DC quantity, imag part is zero
     K{id_x,4}{ctrl_var,1}(idx.vscac_pq,1) = x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + length(idx.vscac_pq) + 1: 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) );
     K{id_x,4}{ctrl_var,1}(idx.vscac_vq,1) = x( 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + length(idx.vscac_vq) + 1 : 2*length(idx.pqac) + 2*length(idx.pvac) + length(idx.pdc) + length(idx.vdc) + 2*length(idx.vscac_pq) + 2*length(idx.vscac_vq) );
     % K{id_x,4}{ctrl_var,1}(idx.vscdc_pq,1) =  %DC quantity, imag part is zero
     % K{id_x,4}{ctrl_var,1}(idx.vscdc_pq,1) =  %DC quantity, imag part is zero

     % Infer complex nodal voltage SCs 
     K{id_x,2}{ctrl_var,1}(:,1) = complex(K{id_x,3}{ctrl_var,1}(:,1) , K{id_x,4}{ctrl_var,1}(:,1)); %complex
     
     K{id_x,5}{ctrl_var,1}(:,1) = (1./abs(E)).*     real(conj(E).*K{id_x,2}{ctrl_var,1}(:,1)) ; %magnitude
     K{id_x,6}{ctrl_var,1}(:,1) = (1./(abs(E).^2)).*imag(conj(E).*K{id_x,2}{ctrl_var,1}(:,1)) ; % angle
     
     
     % Time
     Time.K = [Time.K; toc(T)];

    end
 end

end