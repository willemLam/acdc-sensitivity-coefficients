%% main script for linear power system state estimation

% clear all;
close all;
% clc;

addpath(genpath(pwd))
addpath '/Users/willem/Documents/phd/Git clones/acdc_power_flow'

%% Base values
A_b = 1e5;
V_b= 400;
Y_b = A_b/V_b^2; 

Vdc_b = 800;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 

%% Set the Grid parameters
Grid_para.n_dc = 8;
Grid_para.n_ac = 18+4;
Grid_para.n_AFE = 4;
Grid_para.n_ph = 3;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;
Grid_para.A_b = A_b;
Grid_para.Vdc_b = Vdc_b;
Grid_para.Ydc_b = Ydc_b;
Grid_para.Adc_b = Adc_b;


%% Set the nodes types
idx1.slack = 1;
idx1.pqac = [2:18]';
idx1.pvac = []';

idx1.pdc = [27:30]';
idx1.vdc = []';

idx1.vscac_pq = []';
idx1.vscac_vq = [19:22]';

idx1.vscdc_pq = []';
idx1.vscdc_vq = [23:26]';


%% Set the filter parameters
Filter_para.R = 0.016*Y_b; %checked 
Filter_para.X = 0.000127325*2*pi*50*Y_b;  %checked
Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Include_losses = 1;

%% Get the EMTP measurements
repeat=1;
ZIN_polyphase = [];
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);
[Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC_augmented(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);


V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
Vdc_LF = transpose(Vdc_LF);

I_complex_LF = transpose(complex(Nodal_I_mag.*cos(Nodal_I_angle), Nodal_I_mag.*sin(Nodal_I_angle)));
Idc_LF = transpose(Idc_inj/2);

I_Flow = complex(Flow_I_mag.*cos(Flow_I_angle), Flow_I_mag.*sin(Flow_I_angle));

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
            z2 = 1000:1045; %1000:1045
        otherwise
            disp('NOPE') 
    end

E_star = mean([V_complex_LF(:,z1); Vdc_LF(:,z1)],2);
S_star = mean([transpose(complex(Nodal_P(z1,:), Nodal_Q(z1,:))); transpose(complex(Pdc_inj(z1,:),0))],2);
I_star = [mean(I_complex_LF(:,z1),2);mean(Idc_LF(:,z1),2)];

E_star2 = mean([V_complex_LF(:,z2); Vdc_LF(:,z2)],2);
S_star2 = mean([transpose(complex(Nodal_P(z2,:), Nodal_Q(z2,:))); transpose(complex(Pdc_inj(z2,:),0))],2);

IFlow = mean(I_Flow(:,z1),2);
% 
% I_star - Grid_para.YY*E_star.*multiplier
%% Augment YY to include the filter and IGBT losses into the admitance matrix


% Text = 'linedata_AC.txt';
% idx1.slack = 1;
% idx1.pqac = [2:14]';
% idx1.pvac = []';
% 
% idx1.pdc = [24:26]';
% idx1.vdc = []';
% 
% idx1.vscac_pq = []';
% idx1.vscac_vq = [15:18]';
% 
% idx1.vscdc_pq = []';
% idx1.vscdc_vq = [19:23]';
% type = 0;
% [Yac, A, Zloss,Zfilter,linedata] = include_losses_filter_in_Y(Text,Grid_para,Filter_para,idx1,E_star([1:54,67:74]),type);
% Yac3 = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
idx3 = Get_multiphase_Node_indices(idx1,Grid_para);
[Grid_para, linedata2] = Get_YY('linedata_AC.txt', 'linedata_DC.txt' , Grid_para, Filter_para, idx3, IFlow, idx1);


Grid_para = Get_Converter_para(idx1,linedata2,Grid_para);

% Grid_para.YY = Y_old;
% E_star = E_old;
S_star = E_star.*conj(Grid_para.YY*E_star);
S_star2 = E_star2.*conj(Grid_para.YY*E_star2);

%% Compute SC
idxCtrl = 1:Grid_para.n_nodes; %estimate all voltage sensitivity coefficients
[K, Time] = SC_voltage_rectangular(E_star,idx3,Grid_para,idxCtrl);

J_PR = zeros(Grid_para.n_nodes);
J_PX = zeros(Grid_para.n_nodes);
J_QR = zeros(Grid_para.n_nodes);
J_QX = zeros(Grid_para.n_nodes);
J_VR = zeros(Grid_para.n_nodes);
J_VX = zeros(Grid_para.n_nodes);
for k = 1:size(K,1)
    if( sum( K{k,1} == idx3.slack))
        continue
    elseif( sum( K{k,1} == idx3.pqac))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3.pvac ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_VR(:,k) = real(K{k,2}{2,1});
        J_VX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3.pdc ) )
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3.vdc ) )
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3.vscac_pq))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3.vscac_vq))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3.vscdc_pq ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3.vscdc_vq ))
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
    else
        warning('somethings off mate')
    end
end

multiplier = [3*ones(length(E_star)-Grid_para.n_dc,1) ; ones(Grid_para.n_dc,1)];
%% Analyse results
disp(mode)
switch mode
case 'P23'
    i = Grid_para.n_ac*(Grid_para.n_ph -1) + 27; %31
    r = complex(J_PR(:,i),J_PX(:,i)).*multiplier;
    c = (E_star -E_star2 )/real(S_star(i)-S_star2(i));
    
case 'E22'
    i  = Grid_para.n_ac*(Grid_para.n_ph -1) + 26; %30
    r = complex(J_VR(:,i ),J_VX(:,i )).*multiplier;
    c = (E_star -E_star2 )/real(E_star(i)-E_star2(i));
    
case 'Q18'
   
    i  = polyphase_indices(22,Grid_para.n_ph); %22
    r = complex(sum(J_QR(:,i ),2),sum(J_QX(:,i ),2)).*multiplier;
    c = (E_star -E_star2 )/imag(S_star(i(1))-S_star2(i(1)));

case 'Q9'
    i = polyphase_indices(9,Grid_para.n_ph);
    r = complex(sum(J_QR(:,i),2),sum(J_QX(:,i),2));
    c = (E_star -E_star2 )/imag(S_star(i(1))-S_star2(i(1)));
%     c = (E_star -E_star 2)/imag(S_star (i(1))-S_star 2(i(1)));

    
case 'P9'
    i = polyphase_indices(9,Grid_para.n_ph);
    r = complex(sum(J_PR(:,i),2),sum(J_PX(:,i),2));
    c = (E_star - E_star2 )/real(S_star(i(1))-S_star2(i(1)));
%     c = (E_star -E_star 2)/real(S_star (i(1))-S_star 2(i(1)));
end

%% Jacobian based
tol = 1e-7;
% n_max = 100;
% E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para .n_dc,1)];
% S_star  = E_star_augmented.*conj(YY_augmented*E_star_augmented);
% Filter_para.Exclude_losses = 1;
% [E,J,n_iter] = NR_rectangularACDC_3ph_general(Grid_para_augmented,Filter_para,S_star_augmented,E_star_augmented,E_0,idx3_augmented,tol,n_max);
% Jinv = inv(J);

% Jinv(1:67,1:51)
% Jinv(68:end,1:51)
% 
% Jinv(1:67,52:102)
% Jinv(68:end,52:102)
% 
% Jinv(1:67,103:114)
% Jinv(68:end,103:114)
% 
% Jinv(1:67,115:126)
% Jinv(68:end,115:126)
% 
% Jinv(1:67,127:130)
% Jinv(68:end,127:130)

%% Make figure
f = Make_scatter_plot(r,c,Grid_para.n_nodes,mode);

%% Save plot
% folder = './Plots/figures';
% saveas(f,[folder filesep() strcat('SC_for_X_',mode)],'epsc');

end

