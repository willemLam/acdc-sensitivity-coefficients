%% main script for linear power system state estimation

% clear all;
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
Grid_para.n_ac = 22;
Grid_para.n_ph = 3;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;
Grid_para.A_b = A_b;

[Yac, YYLac, YLac, YTac, YYTac, I_bac, Ampacitiesac, y_lx, y_tx, Aac, linedata_ac]  = Ymatrix_non_sym_shunt('linedata_AC_LCL.txt',A_b,V_b,[]);
[Ydc, YYLdc, YLdc, YTdc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix_non_sym_shunt('linedata_DC_LCL.txt',Adc_b,Vdc_b,[]);

Ydc = Ydc(23:30,23:30)/2;
YYLdc = YYLdc(23:30,23:30)/2;
YYTdc = YYTdc(23:30,23:30)/2;
Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
YYLac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),YYLac,'UniformOutput',false));
YYTac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),YYTac,'UniformOutput',false));
YY = blkdiag(Yac,Ydc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);
Grid_para.YY = YY;

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

idx3 = Get_multiphase_Node_indices(idx1,Grid_para);
linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para(idx1,linedata,Grid_para);

%% Get the EMTP measurements
repeat=1;
ZIN_polyphase = [];
[Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC_LCL(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC_augmented(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);


V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
Vdc_LF = transpose(Vdc_LF);

I_complex_LF = transpose(complex(Nodal_I_mag.*cos(Nodal_I_angle), Nodal_I_mag.*sin(Nodal_I_angle)));
Idc_LF = transpose(Idc_inj/2);

% balanced
modes = {'P23';'E22';'Q18';'Q9';'P9'};
modes = {'P9';'Q9';'Q18';'E22';'P23'};
% unbalaced
% modes = {'P9';'Q9'};
data = [];

%% LF solution
% Set the filter parameters
Filter_para.R = 0*0.008*Y_b; %checked
Filter_para.X = 0*0.04*Y_b;  %checked
Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Exclude_losses = 1;
% Initialize
E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
% Solve the PF
tol = 1e-7;
n_max = 100;



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

E_star2 = mean([V_complex_LF(:,z2); Vdc_LF(:,z2)],2);
S_star2 = mean([transpose(complex(Nodal_P(z2,:), Nodal_Q(z2,:))); transpose(complex(Pdc_inj(z2,:),0))],2);

[E,J,~] = NR_rectangularACDC_3ph_general(Grid_para,Filter_para,S_star,E_star,E_0,idx3,tol,n_max);
[E2,J2,~] = NR_rectangularACDC_3ph_general(Grid_para,Filter_para,S_star2,E_star2,E_0,idx3,tol,n_max);

Jinv = inv(J);
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

%% Analyse results
disp(mode)
switch mode
case 'P23'
    i = 71; 
    r = complex(J_PR(:,i),J_PX(:,i));
    c = (E_star-E_star2)/real(S_star(i)-S_star2(i));
    lf = (E-E2)/real(S_star(i)-S_star2(i));
    j_re = zeros(n_nodes,1);
    j_re([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq;idx.vscdc_pq;idx.pdc]) = Jinv(1:67,127);
    j_im = zeros(n_nodes,1);
    j_im([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq]) = Jinv(68:end,127);
    j = complex(j_re,j_im);
    
    
case 'E22'
    i = 70; 
    r = complex(J_VR(:,i),J_VX(:,i));
    c = (E_star-E_star2)/real(E_star(i)-E_star2(i));
    lf = (E-E2)/real(E_star(i)-E_star2(i));
    j_re = zeros(n_nodes,1);
    j_re([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq;idx.vscdc_pq;idx.pdc]) = Jinv(1:67,104);
    j_im = zeros(n_nodes,1);
    j_im([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq]) = Jinv(68:end,104);
    j = complex(j_re,j_im);
    
case 'Q18'
    i = polyphase_indices(22,Grid_para.n_ph); %22
    r = complex(sum(J_QR(:,i),2),sum(J_QX(:,i),2));
    c = (E_star-E_star2)/imag(S_star(i(1))-S_star2(i(1)));
    lf = (E-E2)/imag(S_star(i(1))-S_star2(i(1)));
    j_re = zeros(n_nodes,1);
    j_re([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq;idx.vscdc_pq;idx.pdc]) = Jinv(1:67,125);
    j_im = zeros(n_nodes,1);
    j_im([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq]) = Jinv(68:end,125);
    j = complex(j_re,j_im);
    
case 'Q9'
    i = polyphase_indices(9,Grid_para.n_ph);
    r = complex(sum(J_QR(:,i),2),sum(J_QX(:,i),2));
    c = (E_star-E_star2)/imag(mean(S_star(i(:))-S_star2(i(:))));
    lf = (E-E2)/imag(mean(S_star(i(:))-S_star2(i(:))));
    j_re = zeros(n_nodes,1);
    j_re([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq;idx.vscdc_pq;idx.pdc]) = sum(Jinv(1:67,73:75),2);
    j_im = zeros(n_nodes,1);
    j_im([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq]) = sum(Jinv(68:end,73:75),2);
    j = complex(j_re,j_im);
    
    
case 'P9'
    i = polyphase_indices(9,Grid_para.n_ph);
    r = complex(sum(J_PR(:,i),2),sum(J_PX(:,i),2));
    c = (E_star-E_star2)/real(mean(S_star(i(:))-S_star2(i(:))));
    lf = (E-E2)/real(mean(S_star(i(:))-S_star2(i(:))));
    j_re = zeros(n_nodes,1);
    j_re([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq;idx.vscdc_pq;idx.pdc]) = sum(Jinv(1:67,22:24),2);
    j_im = zeros(n_nodes,1);
    j_im([idx.pqac;idx.pvac;idx.vscac_pq;idx.vscac_vq]) = sum(Jinv(68:end,22:24),2);
    j = complex(j_re,j_im);
    
end

%% Make figure
f = Make_scatter_plot(r,c,lf,j,Grid_para.n_nodes,mode);

data(m,:) = [mean((abs(r) - abs(c))), max((abs(r) - abs(c))), min((abs(r) - abs(c))) ] 

%% Save plot
% folder = './Plots/figures';
% saveas(f,[folder filesep() strcat('SC_for_X_',mode)],'epsc');

end


