%% main script for linear power system state estimation

clear all;
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
Grid_para.n_dc = 2;
Grid_para.n_ac = 3;
Grid_para.n_AFE = 1;
Grid_para.n_ph = 1;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;
Grid_para.A_b = A_b;
Grid_para.Vdc_b = Vdc_b;
Grid_para.Ydc_b = Ydc_b;
Grid_para.Adc_b = Adc_b;


%% Set the nodes types
idx1.slack = 1;
idx1.pqac = []'; %2
idx1.pvac = [2]';

idx1.pdc = [5]';
idx1.vdc = []';

idx1.vscac_pq = []';
idx1.vscac_vq = [3]';

idx1.vscdc_pq = []';
idx1.vscdc_vq = [4]';

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
Filter_para.Include_losses = 0;
% Initialize
E_0 = [repmat([1], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
% Solve the PF
tol = 1e-7;
n_max = 100;

%% Get grid YY
[Yac, YYLac, YLac, YTac, YYTac, I_bac, Ampacitiesac, y_lx, y_tx, Aac, linedata_ac]  = Ymatrix('linedata_AC_small.txt',A_b,V_b,[]);
[Ydc, YYLdc, YLdc, YTdc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC_small.txt',Adc_b,Vdc_b,[]);

Ydc = Ydc(4:end,4:end)/2;
YYLdc = YYLdc(4:end,4:end)/2;
YYTdc = YYTdc(4:end,4:end)/2;
YY = blkdiag(Yac,Ydc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);
Grid_para.YY = YY;
linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para(idx1,linedata,Grid_para);


%% Start
modes = {'P2';'V2';'Q3';'E4';'P5'}; %{'P2';'Q2';'Q3';'E4';'P5'};
for m=1:length(modes)
    
    mode = char(modes(m));
    
    S_star = [0;complex(0.4,0.2);complex(0,-0.4);0;-0.5];
    E_star = [complex(1);complex(0.98,0.1);0;1;0];
    
    [E_star,J,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star,E_star,E_0,idx1,tol,n_max);
    S_star = E_star.*conj(Grid_para.YY*E_star);
    
    
    
   


    %% balanced
    switch mode
        case 'P2'
            S_star2 = [0;complex(0.41,0.2);complex(0,-0.4);0;-0.5];
            E_star2 = [complex(1);complex(0.98,0.1);0;1;0];
        case 'Q2'
            S_star2 = [0;complex(0.4,0.21);complex(0,-0.4);0;-0.5];
            E_star2 = [complex(1);complex(0.98,0.1);0;1;0];
        case 'V2'
            S_star2 = [0;complex(0.4,0.2);complex(0,-0.4);0;-0.5];
            E_star2 = [complex(1);1.001*complex(0.98,0.1);0;1;0];
        case 'Q3'
            S_star2 = [0;complex(0.4,0.2);complex(0,-0.41);0;-0.5];
            E_star2 = [complex(1);complex(0.98,0.1);0;1;0];
        case 'E4'
            S_star2 = [0;complex(0.4,0.2);complex(0,-0.4);0;-0.5];
            E_star2 = [complex(1);complex(0.98,0.1);0;1.01;0];
        case 'P5'
            S_star2 = [0;complex(0.4,0.2);complex(0,-0.4);0;-0.51];
            E_star2 = [complex(1);complex(0.98,0.1);0;1;0];
        otherwise
            disp('NOPE') 
    end


    [E_star2,J2,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star2,E_star2,E_0,idx1,tol,n_max);
    S_star2 = E_star2.*conj(Grid_para.YY*E_star2);

%% Compute SC
idxCtrl = 1:Grid_para.n_nodes; %estimate all voltage sensitivity coefficients
[K, Time] = SC_voltage_rectangular(E_star,idx1,Grid_para,idxCtrl);

J_PR = zeros(Grid_para.n_nodes);
J_PX = zeros(Grid_para.n_nodes);
J_QR = zeros(Grid_para.n_nodes);
J_QX = zeros(Grid_para.n_nodes);
J_VR = zeros(Grid_para.n_nodes);
J_VX = zeros(Grid_para.n_nodes);
for k = 1:size(K,1)
    if( sum( K{k,1} == idx1.slack))
        continue
    elseif( sum( K{k,1} == idx1.pqac))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1.pvac ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_VR(:,k) = real(K{k,2}{2,1});
        J_VX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1.pdc ) )
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx1.vdc ) )
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx1.vscac_pq))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1.vscac_vq))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
        J_QR(:,k) = real(K{k,2}{2,1});
        J_QX(:,k) = imag(K{k,2}{2,1});
    elseif( sum( K{k,1} == idx1.vscdc_pq ))
        J_PR(:,k) = real(K{k,2}{1,1});
        J_PX(:,k) = imag(K{k,2}{1,1});
    elseif( sum( K{k,1} == idx1.vscdc_vq ))
        J_VR(:,k) = real(K{k,2}{1,1});
        J_VX(:,k) = imag(K{k,2}{1,1});
    else
        warning('somethings off mate')
    end
end

%% Analyse results
disp(mode)
switch mode
case 'P2'
    i = 2; %31
    r = complex(J_PR(:,i),J_PX(:,i));
    c = (E_star -E_star2 )/real(S_star(i)-S_star2(i));
    
case 'V2'
    i = 2;
    r = complex(J_VR(:,i ),J_VX(:,i ));
    c = (E_star -E_star2 )/(abs(E_star(i))-abs(E_star2(i)));

case 'Q2'
    i = 2;
    r = complex(sum(J_QR(:,i),2),sum(J_QX(:,i),2));
    c = (E_star -E_star2 )/imag(S_star(i(1))-S_star2(i(1)));
%     c = (E_star -E_star 2)/imag(S_star (i(1))-S_star 2(i(1)));

case 'Q3'
    i = 3; %22
    r = complex(sum(J_QR(:,i ),2),sum(J_QX(:,i ),2));
    c = (E_star -E_star2 )/imag(S_star(i(1))-S_star2(i(1)));

case 'E4'
    i = 4; %30
    r = complex(J_VR(:,i ),J_VX(:,i ));
    c = (E_star -E_star2 )/real(E_star(i)-E_star2(i));

case 'P5'
    i = 5;
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
f = Make_scatter_plot_small(r,c,Grid_para.n_nodes,mode);


data(m,:) = [mean((abs(r) - abs(c))), max((abs(r) - abs(c))), min((abs(r) - abs(c))) ] 

%% Save plot
% folder = './Plots/figures';
% saveas(f,[folder filesep() strcat('SC_for_X_',mode)],'epsc');

end

