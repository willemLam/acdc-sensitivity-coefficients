%% main script for linear power system state estimation


%clear all;

% clear all;
close all;
% clc;

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
linedata_dc(:,1:2) = linedata_dc(:,1:2)-4;
Ydc = Ydc(23:30,23:30)/2;
Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
YY = blkdiag(Yac,Ydc);

Grid_para.G = real(YY);
Grid_para.B = imag(YY);
Grid_para.YY = YY;
Grid_para.Yac = Yac;
Grid_para.Ydc = Ydc;

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

%% Set the filter parameters
Filter_para.R = 10e-4 + 0*0.016*Y_b; %checked 
Filter_para.X = 10e-4 + 0*0.000127325*2*pi*50*Y_b;  %checked
% Filter_para.R = 0.04*Y_b; %checked %0.08
% Filter_para.X = 0.04*Y_b;  %checked
Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Include_losses = 0;

%% Get the EMTP measurements
repeat=1;
ZIN_polyphase = [];
[Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);
% [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC_augmented(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);
%[Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc_LF, n_timesteps, Mreal, Mimag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M_LF] = GetEMTPdata_SC_augmented(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN_polyphase,Grid_para.n_ph);

V_complex_LF = transpose(complex(Nodal_V_mag.*cos(Nodal_V_angle), Nodal_V_mag.*sin(Nodal_V_angle)));
Vdc_LF = transpose(Vdc_LF);

I_complex_LF = transpose(complex(Nodal_I_mag.*cos(Nodal_I_angle), Nodal_I_mag.*sin(Nodal_I_angle)));
Idc_LF = transpose(Idc_inj/2);

% balanced
modes = {'P23';'E22';'Q18';'Q9';'P9'}%;'P11a'};

modes = {'P23'}%;'E22';'Q18';'Q9';'P9'};

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
        case 'P11a'
            z2 = 1130:1200;
        otherwise
            disp('NOPE') 
    end

E_star = mean([V_complex_LF(:,z1); Vdc_LF(:,z1)],2);
S_star = mean([transpose(complex(Nodal_P(z1,:), Nodal_Q(z1,:))); transpose(complex(Pdc_inj(z1,:),0))],2);

E_star2 = mean([V_complex_LF(:,z2); Vdc_LF(:,z2)],2);
S_star2 = mean([transpose(complex(Nodal_P(z2,:), Nodal_Q(z2,:))); transpose(complex(Pdc_inj(z2,:),0))],2);

alp = exp(2*pi/3*1i);
A = 1/3 *[1 1     1; 
         1 alp   alp^2; 
         1 alp^2 alp];
ACell =  repmat({A}, 1, Grid_para.n_ac+4);
ICell =  repmat({1}, 1, Grid_para.n_dc); 
Atot = blkdiag(ACell{:},ICell{:});

%% Augment YY to include the filter and IGBT losses into the admitance matrix


% 
type = 0; %equivalent PI with filter and losses
[Yac_augmented, A_augmented, Zloss,Zfilter,linedata1] = include_losses_filter_in_Y('linedata_AC.txt',Grid_para,Filter_para,idx1,E_star,type);
% Yac_augmented = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac_augmented,'UniformOutput',false));
% YY_augmented = blkdiag(Yac_augmented,Ydc);
% 
% Grid_para_augmented = Grid_para;
% 
% Grid_para_augmented.YY = YY_augmented;
% Grid_para_augmented.Yac = Yac_augmented;
% Grid_para_augmented.Ydc = Ydc;
% Grid_para_augmented.n_ac = 18+Grid_para_augmented.n_AFE; %Every AFE gets one extra virtual node
% Grid_para_augmented.n_nodes = Grid_para_augmented.n_ac*Grid_para_augmented.n_ph + Grid_para_augmented.n_dc;
% 
% %% Augment E_star and S_star to include the filter and IGBT losses into the admitance matrix
Zloss = cell2mat(arrayfun(@(x) x*ones(Grid_para.n_ph,1),Zloss,'UniformOutput',false));
Zfilter = cell2mat(arrayfun(@(x) x*ones(Grid_para.n_ph,1),Zfilter,'UniformOutput',false));
% 
% E_augment_loss = E_star(sort([idx3.vscac_pq;idx3.vscac_vq])) + (Zloss+Zfilter).*(YY(sort([idx3.vscac_pq;idx3.vscac_vq]),:)*E_star);
% E_star_augmented = [E_star(1:Grid_para.n_ac*Grid_para.n_ph); E_augment_loss; E_star(Grid_para.n_ac*Grid_para.n_ph+1 : end)  ];
% 
E_2_augment_loss = E_star2(sort([idx3.vscac_pq;idx3.vscac_vq])) + (Zloss+Zfilter).*(YY(sort([idx3.vscac_pq;idx3.vscac_vq]),:)*E_star2);
E_star2_augmented = [E_star2(1:Grid_para.n_ac*Grid_para.n_ph); E_2_augment_loss; E_star2(Grid_para.n_ac*Grid_para.n_ph+1 : end)  ];

[E_star_augmented, YY_augmented, Grid_para_augmented] = Augment_system_with_filter_n_losses(E_star, Grid_para, Filter_para, idx1);

%% Augment indices
idx1_augmented.slack = 1;
idx1_augmented.pqac = [2:18]';
idx1_augmented.pvac = []';

idx1_augmented.pdc = [27:30]';
idx1_augmented.vdc = []';

idx1_augmented.vscac_pq = []';
idx1_augmented.vscac_vq = [19:22]';

idx1_augmented.vscdc_pq = []';
idx1_augmented.vscdc_vq = [23:26]';

idx3_augmented = Get_multiphase_Node_indices(idx1_augmented,Grid_para_augmented);


%% Compute SC
idxCtrl = 1:Grid_para_augmented.n_nodes; %estimate all voltage sensitivity coefficients

[K, Time] = SC_voltage_rectangular(E_star_augmented,idx3_augmented,Grid_para_augmented,idxCtrl);

% time_ar = [];
% for i = 1:100
% [K, Time] = SC_voltage_rectangular(E_star_augmented,idx3_augmented,Grid_para_augmented,idxCtrl);
% time_ar = [time_ar,Time.F];
% end
% mean(time_ar)

J_PR = zeros(Grid_para_augmented.n_nodes);
J_PX = zeros(Grid_para_augmented.n_nodes);
J_QR = zeros(Grid_para_augmented.n_nodes);
J_QX = zeros(Grid_para_augmented.n_nodes);
J_VR = zeros(Grid_para_augmented.n_nodes);
J_VX = zeros(Grid_para_augmented.n_nodes);
for k = 1:size(K,1)
    if( sum( K{k,1} == idx3_augmented.slack))
        continue
    elseif( sum( K{k,1} == idx3_augmented.pqac))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3_augmented.pvac ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_VR(:,k) = real(K{k,2}{2,1});
        J_VX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3_augmented.pdc ) )
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3_augmented.vdc ) )
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3_augmented.vscac_pq))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3_augmented.vscac_vq))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx3_augmented.vscdc_pq ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx3_augmented.vscdc_vq ))
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
    else
        warning('somethings off mate')
    end
end

multiplier = [3*ones(length(E_star_augmented)-Grid_para.n_dc,1) ; ones(Grid_para.n_dc,1)];
%% Analyse results
disp(mode)
switch mode
case 'P23'
    i = Grid_para.n_ac*(Grid_para.n_ph -1) + 23; 
    i_augmented = Grid_para_augmented.n_ac*(Grid_para_augmented.n_ph -1) + 27; %31
    r = complex(J_PR(:,i_augmented),J_PX(:,i_augmented)).*multiplier;
    c = (E_star_augmented-E_star2_augmented)/real(S_star(i)-S_star2(i));
    
case 'E22'
    i = Grid_para.n_ac*(Grid_para.n_ph -1) + 22; 
    i_augmented = Grid_para_augmented.n_ac*(Grid_para_augmented.n_ph -1) + 26; %30
    r = complex(J_VR(:,i_augmented),J_VX(:,i_augmented)).*multiplier;
    c = (E_star_augmented-E_star2_augmented)/real(E_star(i)-E_star2(i));
    
case 'Q18'
    i = polyphase_indices(18,Grid_para.n_ph); %22
    i_augmented = polyphase_indices(22,Grid_para_augmented.n_ph); %22
    r = complex(sum(J_QR(:,i_augmented),2),sum(J_QX(:,i_augmented),2)).*multiplier;
    c = (E_star_augmented-E_star2_augmented)/imag(S_star(i(1))-S_star2(i(1)));

case 'Q9'
    i = polyphase_indices(9,Grid_para_augmented.n_ph);
    r = complex(sum(J_QR(:,i),2),sum(J_QX(:,i),2));
    c = (E_star_augmented-E_star2_augmented)/imag(S_star(i(1))-S_star2(i(1)));

    
case 'P9'
    i = polyphase_indices(9,Grid_para_augmented.n_ph);
    r = complex(sum(J_PR(:,i),2),sum(J_PX(:,i),2));
    c = (E_star_augmented-E_star2_augmented)/real(S_star(i(1))-S_star2(i(1)));
    
case 'P11a'
    i = polyphase_indices(11,Grid_para_augmented.n_ph);

    r = complex(J_PR(:,i),J_PX(:,i));
    i = i(1);
    c = (E_star_augmented-E_star2_augmented)/real(S_star(i)-S_star2(i));
    
end

%% Make figure
f = Make_scatter_plot(r,c,Grid_para_augmented.n_nodes,mode);

mean((abs(r) - abs(c)))*10000
max((abs(r) - abs(c)))*10000
min((abs(r) - abs(c)))*10000
%% Save plot
% folder = './Plots/figures';
% saveas(f,[folder filesep() strcat('SC_for_X_',mode)],'epsc');

end
