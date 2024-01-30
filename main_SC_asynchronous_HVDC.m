%% Valiation of the SC for the high voltage transmission network.
%% The network consists of two asynchronous AC networks that are coupled by two HVDC networks
clear all;
close all;
%clc;

addpath(genpath(pwd))

load('data_HVDC.mat')
idx1 = data.idx;

idx1.slack = [1,64]'; 
idx1.pvac = [2,3,6,8,9,12,66,69,71]';
idx1.pqac = setdiff([2:58,65:77]',idx1.pvac); 

idx1.pdc = [81,82,83,87,88,89,90,91,92,93]';
idx1.vdc = []';

idx1.vscac_pq = [80,61,78,79,62,63]';
idx1.vscac_vq = [59,60]';

idx1.vscdc_pq = [85,86,95,96,97,98]';
idx1.vscdc_vq = [84,94]';



Grid_para = data.Grid_para;
S = data.S_star;
E = data.E_star;

[Yac, YYLac, YLac, YTac, YYTac, I_bac, Ampacitiesac, y_lx, y_tx, Aac, linedata_ac]  = Ymatrix('linedata_AC_HVDC.txt',Grid_para.A_b*0+1,Grid_para.V_b*0+1,[]);
[Ydc, YYLdc, YLdc, YTdc, YYTdc, I_bdc, Ampacitiesdc, y_ihdc, y_idc, Adc, linedata_dc]  = Ymatrix('linedata_DC_HVDC.txt',Grid_para.Adc_b*0+1,Grid_para.Vdc_b*0+1,[]);

Ydc = Ydc(Grid_para.n_ac+1:end,Grid_para.n_ac+1:end)/2;
YYLdc = YYLdc(Grid_para.n_ac+1:end,Grid_para.n_ac+1:end)/2;
YYTdc = YYTdc(Grid_para.n_ac+1:end,Grid_para.n_ac+1:end)/2;
YY = blkdiag(Yac,Ydc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);
Grid_para.YY = YY;
linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para(idx1,linedata,Grid_para);




%% Set the filter parameters
Filter_para.R = 0.008*Grid_para.Y_b; %checked
Filter_para.X = 0.04*Grid_para.Y_b;  %checked
Filter_para.IGBT_piecewise = [  0                   0
                                0.04926559627563	0.7
                                2.30625399327864	0.8
                                15.7793399043317	0.85
                                107.547461516782	0.8999
                                735.837403888342	0.9499
                                1588.01477341768	0.9699];
Filter_para.Exclude_losses = 1;


%% Initialize
if Grid_para.n_ph == 3
    E_0 = [repmat([1; exp(-1i*pi*2/3);  exp(1i*pi*2/3) ], Grid_para.n_ac,1 ) ; ones(Grid_para.n_dc,1)];
elseif Grid_para.n_ph == 1
    E_0 = [ones(Grid_para.n_ac + Grid_para.n_dc,1)];
end

tol = 1e-6;
n_max = 20;

[E_star,J,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S,E,E_0,idx1,tol,n_max);
S_star = E_star.*conj(Grid_para.YY*E_star);

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

%% Start
modes = {'P_pvac', 'V_pvac', 'P_pqac', 'Q_pqac', 'Q_qvic', 'V_qvic', 'P_pqic', 'Q_pqic', 'P_pdc'};
modes_length = [length(idx1.pvac), length(idx1.pvac), length(idx1.pqac), length(idx1.pqac), length(idx1.vscac_vq), length(idx1.vscdc_vq), length(idx1.vscac_pq), length(idx1.vscac_pq),length(idx1.pdc)];

for m=1:length(modes)
    
    mode = char(modes(m));
    
    for i = 1:modes_length(m)
    %% balanced
    
    S_star2 = S_star;
    E_star2 = E_star;
    
    ad = 1e-8;
    switch mode
        
        case 'P_pvac'
            S_star2(idx1.pvac(i)) = complex(real(S_star(idx1.pvac(i)))+ad, imag(S_star(idx1.pvac(i)))) ;
        case 'V_pvac'
            E_star2(idx1.pvac(i)) = abs(E_star(idx1.pvac(i)))*(1+ad) * exp(1i*angle(E_star(idx1.pvac(i))));

        case 'P_pqac'
            S_star2(idx1.pqac(i)) = complex(real(S_star(idx1.pqac(i)))+ad, imag(S_star(idx1.pqac(i)))) ;
        case 'Q_pqac'
            S_star2(idx1.pqac(i)) = complex(real(S_star(idx1.pqac(i))), imag(S_star(idx1.pqac(i)))+ad) ;
            
        case 'Q_qvic'
            S_star2(idx1.vscac_vq(i)) = complex(real(S_star(idx1.vscac_vq(i))), imag(S_star(idx1.vscac_vq(i)))+ad) ;
        case 'V_qvic'
            E_star2(idx1.vscdc_vq(i)) = E_star(idx1.vscdc_vq(i))*(1+ad) ;

        case 'P_pqic'
            S_star2(idx1.vscac_pq(i)) = complex(real(S_star(idx1.vscac_pq(i)))+ad, imag(S_star(idx1.vscac_pq(i)))) ;
        case 'Q_pqic'
            S_star2(idx1.vscac_pq(i)) = complex(real(S_star(idx1.vscac_pq(i))), imag(S_star(idx1.vscac_pq(i)))+ad) ;

        case 'P_pdc'
            S_star2(idx1.pdc(i)) = S_star(idx1.pdc(i))+ad ;
            
        otherwise
            disp('NOPE') 
            
    end

    [E_star2,J2,~] = NR_rectangularACDC_1ph_general(Grid_para,Filter_para,S_star2,E_star2,E_0,idx1,tol,n_max);
    S_star2 = E_star2.*conj(Grid_para.YY*E_star2);



%% Analyse results
disp(mode)
switch mode
    
    
case 'P_pvac'
    l = idx1.pvac(i);
    r = complex(J_PR(:,l),J_PX(:,l));
    c = (E_star -E_star2 )/real(S_star(l)-S_star2(l));
    
%     Make_scatter_plot_small(r,c,Grid_para.n_nodes,mode);
    
    data_real.P_pvac(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.P_pvac(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];
     
case 'V_pvac'
    l = idx1.pvac(i);
    r = complex(J_VR(:,l),J_VX(:,l));
    c = (E_star -E_star2 )/(abs(E_star(l))-abs(E_star2(l)));
%     Make_scatter_plot_small(r,c,Grid_para.n_nodes,mode);
    
    data_real.V_pvac(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.V_pvac(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];
    
case 'P_pqac'
    l = idx1.pqac(i); 
    r = complex(sum(J_PR(:,l),2),sum(J_PX(:,l),2));
    c = (E_star -E_star2 )/real(S_star(l)-S_star2(l));
    
%     Make_scatter_plot_small(r,c,Grid_para.n_nodes,mode);
    
    data_real.P_pqac(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.P_pqac(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];
    
case 'Q_pqac'
    l = idx1.pqac(i);
    r = complex(sum(J_QR(:,l ),2),sum(J_QX(:,l ),2));
    c = (E_star -E_star2 )/imag(S_star(l)-S_star2(l));
    data_real.Q_pqac(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.Q_pqac(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];

case 'Q_qvic'
    l = idx1.vscac_vq(i);
    r = complex(sum(J_QR(:,l ),2),sum(J_QX(:,l ),2));
    c = (E_star -E_star2 )/imag(S_star(l)-S_star2(l));
    data_real.Q_qvic(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.Q_qvic(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];
    
case 'V_qvic'
    l = idx1.vscdc_vq(i);
    r = complex(J_VR(:,l ),J_VX(:,l ));
    c = (E_star -E_star2 )/real(E_star(l)-E_star2(l));
    data_real.V_qvic(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.V_qvic(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];
    
case 'P_pqic'
    l = idx1.vscac_pq(i);
    r = complex(sum(J_PR(:,l),2),sum(J_PX(:,l),2));
    c = (E_star -E_star2 )/real(S_star(l)-S_star2(l));
    data_real.P_pqic(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.P_pqic(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];

case 'Q_pqic'
    l = idx1.vscac_pq(i);
    r = complex(sum(J_QR(:,l ),2),sum(J_QX(:,l ),2));
    c = (E_star -E_star2 )/imag(S_star(l)-S_star2(l));
    data_real.Q_pqic(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.Q_pqic(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];
    
case 'P_pdc'
    l = idx1.pdc(i);
    r = complex(sum(J_PR(:,l),2),sum(J_PX(:,l),2));
    c = (E_star - E_star2 )/real(S_star(l)-S_star2(l));
    data_real.P_pdc(i,:) = [sqrt(mean((real(r) - real(c)).^2)), max(abs(max((real(r) - real(c)))), abs(min((real(r) - real(c))))) ];
    data_imag.P_pdc(i,:) = [sqrt(mean((imag(r) - imag(c)).^2)), max(abs(max((imag(r) - imag(c)))), abs(min((imag(r) - imag(c))))) ];
    
end

%% Make figure
% f = Make_scatter_plot_small(r,c,Grid_para.n_nodes,mode);



%% Save plot
% folder = './Plots/figures';
% saveas(f,[folder filesep() strcat('SC_for_X_',mode)],'epsc');
    end
end


%% post processing
results_PQ =   [mean([data_real.P_pqac(:,1); data_real.Q_pqac(:,1)]),max([data_real.P_pqac(:,2); data_real.Q_pqac(:,2)]); ...
                mean([data_imag.P_pqac(:,1); data_imag.Q_pqac(:,1)]),max([data_imag.P_pqac(:,2); data_imag.Q_pqac(:,2)])];

results_PV =   [mean([data_real.P_pvac(:,1); data_real.V_pvac(:,1)]),max([data_real.P_pvac(:,2); data_real.V_pvac(:,2)]); ...
                mean([data_imag.P_pvac(:,1); data_imag.V_pvac(:,1)]),max([data_imag.P_pvac(:,2); data_imag.V_pvac(:,2)])];

results_ICpq = [mean([data_real.P_pqic(:,1); data_real.Q_pqic(:,1)]),max([data_real.P_pqic(:,2); data_real.Q_pqic(:,2)]); ...
                mean([data_imag.P_pqic(:,1); data_imag.Q_pqic(:,1)]),max([data_imag.P_pqic(:,2); data_imag.Q_pqic(:,2)])];

results_ICvq = [mean([data_real.V_qvic(:,1); data_real.Q_qvic(:,1)]),max([data_real.V_qvic(:,2); data_real.Q_qvic(:,2)]); ...
                mean([data_imag.V_qvic(:,1); data_imag.Q_qvic(:,1)]),max([data_imag.V_qvic(:,2); data_imag.Q_qvic(:,2)])];

results_Pdc =  [mean([data_real.P_pdc(:,1)]),max([data_real.P_pdc(:,2)]); ...
                mean([data_imag.P_pdc(:,1)]),max([data_imag.P_pdc(:,2)])];
            
results_ALL = [results_PQ, results_PV, results_ICpq, results_ICvq, results_Pdc];
T2 = array2table(results_ALL,'VariableNames',{'PQ_rmse','PQ_max','PV_rmse','PV_max','ICpq_rmse','ICpq_max','ICvq_rmse','ICvq_max','Pdc_rmse','Pdc_max'},'RowName',{'SC_real','SC_imag'})


