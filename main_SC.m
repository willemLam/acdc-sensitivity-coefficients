%main script for linear power system state estimation

clear all;
close all;
clc;

addpath(genpath(pwd))

%% Base values
A_b = 1e5;
V_b= 400;
Y_b = A_b/V_b^2; 

Vdc_b = 800;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 

%% Set the Grid parameters
Grid_para.n_dc = 8;
Grid_para.n_ac = 18;
Grid_para.n_ph = 3;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;
Grid_para.A_b = A_b;


[Yac, YYL, YL, YT, YYT, I_b, Ampacities, y_lx, y_tx, A, linedata_ac]  = Ymatrix('linedata_AC.txt',A_b,V_b,[]);
[Ydc, YYLdc, YLdc, YT_dc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC.txt',Adc_b,Vdc_b,[]);

Ydc = Ydc(19:26,19:26)/2;
Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
YY = blkdiag(Yac,Ydc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);
Grid_para.YY = YY;
Grid_para.Yac = Yac;
Grid_para.Ydc = Ydc;


%% Set the filter parameters
Filter_para.R = 0.08*Y_b; %checked
Filter_para.X = 0.04*Y_b;  %checked
Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Exclude_losses = 1;

   
%% Get the EMTP measurements
repeat=1;
ZIN_polyphase = [];
[Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);

V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
Vdc_LF = transpose(Vdc_LF);


%% Set the nodes types
idx1.slack = 1;
idx1.pqac = [2:14]';
idx1.pvac = []';

idx1.pdc = [23:26]';
idx1.vdc = []';

idx1.vscac_pq = []';
idx1.vscac_vq = [15:18]';

idx1.vscdc_pq = []';
idx1.vscdc_vq = [19:22]';

idx3 = Get_multiphase_Node_indices(idx1,Grid_para);
linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para(idx1,linedata,Grid_para);


% balanced
modes = {'P23';'E22';'Q18';'Q9';'P9'};
% unbalaced
% modes = {'P9';'Q9'};


for m=1:length(modes)
    
    mode = char(modes(m));
    
    z1 = 150:295;
    %% unbalanced
    % switch mode
    %     case 'P9'
    %         z2 = 400:445;
    %     case 'Q9'
    %         z2 = 550:595;
    %     otherwise
    %         disp('NOPE') 
    % end

    %% balanced
    switch mode
        case 'P23'
            z2 = 400:445;
        case 'E22'
            z2 = 550:595;
        case 'Q18'
            z2 = 700:745;
        case 'Q9'
            z2 = 850:895;
        case 'P9'
            z2 = 1000:1045;
        otherwise
            disp('NOPE') 
    end

E_star = mean([V_complex_LF(:,z1); Vdc_LF(:,z1)],2);
S_star = mean([transpose(complex(Nodal_P(z1,:), Nodal_Q(z1,:))); transpose(complex(Pdc_inj(z1,:),0))],2);
S_star(sort([idx3.vscac_pq;idx3.vscac_vq])) = -1i*imag(S_star(sort([idx3.vscac_pq;idx3.vscac_vq])));
S_star(sort([idx3.vscdc_pq;idx3.vscdc_vq])) = 0;

E_star2 = mean([V_complex_LF(:,z2); Vdc_LF(:,z2)],2);
S_star2 = mean([transpose(complex(Nodal_P(z2,:), Nodal_Q(z2,:))); transpose(complex(Pdc_inj(z2,:),0))],2);
S_star2(sort([idx3.vscac_pq;idx3.vscac_vq])) = -1i*imag(S_star2(sort([idx3.vscac_pq;idx3.vscac_vq])));
S_star2(sort([idx3.vscdc_pq;idx3.vscdc_vq])) = 0;


alp = exp(2*pi/3*1i);
A = 1/3*[1 1     1; 
         1 alp   alp^2; 
         1 alp^2 alp];
ACell =  repmat({A}, 1, Grid_para.n_ac);
ICell =  repmat({1}, 1, Grid_para.n_dc); 
Atot = blkdiag(ACell{:},ICell{:});
Atot*E_star



filter = 0;
unblanced_3ph  = 0;
idxCtrl = 1:Grid_para.n_nodes;

%% Compute SC
% [K, Time] = SC_Voltage_V5(Yac,S0ac,Eac,Ydc,S0dc,Edc,idx1ph,idx3ph,idxCtrl,nph,vdep,Zf,Fl,unblanced_3ph,filter);
[K, Time] = SC_Voltage_V5_3(S_star,E_star,idx1,idx3,Grid_para,Filter_para,idxCtrl,unblanced_3ph,filter);

J_PR = zeros(Grid_para.n_nodes);
J_PX = zeros(Grid_para.n_nodes);
J_QR = zeros(Grid_para.n_nodes);
J_QX = zeros(Grid_para.n_nodes);
J_VR = zeros(Grid_para.n_nodes);
J_VX = zeros(Grid_para.n_nodes);
for k = 1:size(K,1)
    if( sum( K{k,1} == idx3.pqac ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3.vscac_vq ))
        J_QR(:,k) = real(K{k,2}{1,1});
        J_QX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3.vscdc_vq ))
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3.pdc ) )
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    end
end


%% Analyse results
disp(mode)
switch mode
case 'P23'
    
    i = Grid_para.n_ac*(Grid_para.n_ph -1) + 23;
    r = complex(J_PR(:,i),J_PX(:,i));
    c = (E_star-E_star2)/real(S_star(i)-S_star2(i));
    
case 'E22'
    
    i = Grid_para.n_ac*(Grid_para.n_ph -1) + 22;
    r = complex(J_VR(:,i),J_VX(:,i));
    c = (E_star-E_star2)/real(E_star(i)-E_star2(i));
    x_rec = [real(c(1:54));imag(c(1:54));real(c(55:end))];
    
case 'Q18'
    
    i = polyphase_indices(18,Grid_para.n_ph);
    r = complex(sum(J_QR(:,i),2),sum(J_QX(:,i),2));
    c = (E_star-E_star2)/imag(S_star(i(1))-S_star2(i(1)));
  
case 'Q9'
    
    i = polyphase_indices(9,Grid_para.n_ph);
    r = complex(sum(J_QR(:,i),2),sum(J_QX(:,i),2));
    c = (E_star-E_star2)/imag(S_star(i(1))-S_star2(i(1)));
    x_rec = [real(c(1:54));imag(c(1:54));real(c(55:end))];
    
case 'P9'
    
    i = polyphase_indices(9,Grid_para.n_ph);
    r = complex(sum(J_PR(:,i),2),sum(J_PX(:,i),2));
    c = (E_star-E_star2)/real(S_star(i(1))-S_star2(i(1)));
    
end

%% make plot

f = Make_scatter_plot(r,c,Grid_para.n_nodes,mode)

% folder = './Plots/figures';
% saveas(f,[folder filesep() strcat('SC_for_X_',mode)],'epsc');

end

