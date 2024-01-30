
clear all;
close all;
%clc;

addpath(genpath(pwd))

% t = tic;
% mpc = loadcase('FUBM_test_grid');
mpc = loadcase('fubm_case_57_14_2MTDC_ctrls_EPFL_4_4'); %FUBM_m_grid, fubm_case_57_14_2MTDC_ctrls_EPFL, fubm_case_30_2MTDC_ctrls_vt2_pf_EPFL, fubm_case1354pegase_2MTDC_ctrls_pf_qt_dp_EPFL


%% Get power setpoints
S = complex(mpc.bus(:,3), mpc.bus(:,4));

%Get Q setpoint from vsc
ind_branch_Q = find(mpc.branch(:,17));
mpc.branch(ind_branch_Q,2);
mpc.branch(ind_branch_Q,17);

ind_bus_from_Q = get_index(mpc.bus(:,1), mpc.branch(ind_branch_Q,1));
ind_bus_to_Q   = get_index(mpc.bus(:,1), mpc.branch(ind_branch_Q,2));




S(ind_bus_to_Q) = S(ind_bus_to_Q) + complex(-mpc.branch(ind_branch_Q,14),mpc.branch(ind_branch_Q,17));
S(ind_bus_from_Q) = S(ind_bus_from_Q) + mpc.branch(ind_branch_Q,14);
S = -S/mpc.baseMVA; % matpower uses a different convention



%% Get voltage setpoints 
E = mpc.bus(:,8) .* exp(1i* mpc.bus(:,9));


%Get Vdc setpoint from vsc
ind_branch_Vdc = find(mpc.branch(:,22));
mpc.branch(ind_branch_Vdc,1);
mpc.branch(ind_branch_Q,22);

ind_bus_Vdc = [];
for i = 1:length(ind_branch_Vdc)
ind_bus_Vdc = [ind_bus_Vdc ; find(mpc.bus(:,1)== mpc.branch(ind_branch_Vdc(i),1))];
end

E(ind_bus_Vdc) = mpc.branch(ind_branch_Vdc,22);




%% Base values
A_b = mpc.baseMVA*1e6;
V_b= mpc.bus(1,10)*1e3;
Y_b = A_b/V_b^2; 

Vdc_b = mpc.bus(end,10)*1e3;
Adc_b = A_b;
Ydc_b = Adc_b/Vdc_b^2; 




%% Set the nodes types


%% 57 node
bus_idx_ac = [find(mpc.bus(:,7)== 1); find(mpc.bus(:,7)== 2) ];
bus_idx_dc = [find(mpc.bus(:,7)== 3); find(mpc.bus(:,7)== 4) ];
bus_idx_vsc_dc = [find(mpc.bus(:,7)== 5); find(mpc.bus(:,7)== 6) ];
bus_idx_vsc_ac = [find(mpc.bus(:,7)== 7); find(mpc.bus(:,7)== 8) ];

bus_ac = mpc.bus(bus_idx_ac,1);
bus_dc = mpc.bus(bus_idx_dc,1);
bus_vsc_dc = mpc.bus(bus_idx_vsc_dc,1);
bus_vsc_ac = mpc.bus(bus_idx_vsc_ac,1);

branch_idx_ac = [];
branch_idx_dc = [];
branch_idx_vsc_dc = [];
branch_idx_vsc_ac = [];

for i = 1:length(mpc.branch(:,1))
    
    if any(bus_ac(:) == mpc.branch(i,1)) && any(bus_ac(:) == mpc.branch(i,2)) && mpc.branch(i,11)
        branch_idx_ac = [branch_idx_ac; i];
        
    elseif any(bus_dc(:) == mpc.branch(i,1)) && any(bus_dc(:) == mpc.branch(i,2))  && mpc.branch(i,11)
        branch_idx_dc = [branch_idx_dc; i];

    elseif any(bus_dc(:) == mpc.branch(i,1)) && any(bus_vsc_dc(:) == mpc.branch(i,2))  && mpc.branch(i,11)
        branch_idx_vsc_dc = [branch_idx_vsc_dc; i];
        
    elseif any(bus_ac(:) == mpc.branch(i,1)) && any(bus_vsc_ac(:) == mpc.branch(i,2))  && mpc.branch(i,11)
        branch_idx_vsc_ac = [branch_idx_vsc_ac; i];
    else
        mpc.branch(i,1:2)
    end
end

branch_vsc_dc = mpc.branch(branch_idx_vsc_dc,[1,2]);
branch_vsc_dc = branch_vsc_dc(:);

branch_vsc_ac = mpc.branch(branch_idx_vsc_ac,[1,2]);
branch_vsc_ac = branch_vsc_ac(:);



branch_vscVQ = mpc.branch(find(mpc.branch(:,22)),[1,2]);
branch_vscPQ = mpc.branch(find(mpc.branch(:,14)),[1,2]);

%% Set the Grid parameters
Grid_para.n_dc = length(bus_dc) + length(bus_vsc_dc);
Grid_para.n_ac = length(bus_ac) + length(bus_vsc_ac);
Grid_para.n_ph = 1;
Grid_para.n_nodes = Grid_para.n_ac*Grid_para.n_ph + Grid_para.n_dc;
Grid_para.V_b = V_b;
Grid_para.Y_b = Y_b;
Grid_para.Vdc_b = Vdc_b;
Grid_para.Ydc_b = Ydc_b;


[YYac, linedata_ac, YYLac, YYTac]  =  Ymatrix_matpower(mpc.branch([branch_idx_ac;branch_idx_vsc_ac ],1:5));
[YYdc, linedata_dc, YYLdc, YYTdc]  =  Ymatrix_matpower(mpc.branch([branch_idx_dc;branch_idx_vsc_dc ],1:5));
YYac = YYac(sort([bus_ac;bus_vsc_ac]),sort([bus_ac;bus_vsc_ac]));
YYdc = YYdc(sort([bus_dc;bus_vsc_dc]),sort([bus_dc;bus_vsc_dc]));

YYLac = YYLac(sort([bus_ac;bus_vsc_ac]),sort([bus_ac;bus_vsc_ac]));
YYLdc = YYLdc(sort([bus_dc;bus_vsc_dc]),sort([bus_dc;bus_vsc_dc]));

YYTac = YYTac(sort([bus_ac;bus_vsc_ac]),sort([bus_ac;bus_vsc_ac]));
YYTdc = YYTdc(sort([bus_dc;bus_vsc_dc]),sort([bus_dc;bus_vsc_dc]));


bus_index = mpc.bus(:,1);

YY = blkdiag(YYac,YYdc);
Grid_para.G = real(YY);
Grid_para.B = imag(YY);
Grid_para.YY = YY;



% indices
idx1.slack = intersect(mpc.areas(:,2), bus_ac); %
idx1.pqac = [setdiff(bus_ac,[bus_vsc_ac;idx1.slack;mpc.gen(:,1)])];
idx1.pvac = [setdiff(mpc.gen(:,1),idx1.slack)];

idx1.pdc = setdiff(bus_dc,bus_vsc_dc);
idx1.vdc = []';

idx1.vscac_pq = intersect(branch_vscPQ,bus_vsc_ac);
idx1.vscac_vq = intersect(branch_vscVQ,bus_vsc_ac);

idx1.vscdc_pq = intersect(branch_vscPQ,bus_vsc_dc);
idx1.vscdc_vq = intersect(branch_vscVQ,bus_vsc_dc);

linedata = [linedata_ac;linedata_dc];
Grid_para = Get_Converter_para_FUBM(idx1,linedata,Grid_para,mpc);


idx1.slack = get_index(bus_index, idx1.slack);
idx1.pvac = get_index(bus_index, idx1.pvac);
idx1.pqac = get_index(bus_index, idx1.pqac);
idx1.vscac_pq = get_index(bus_index, idx1.vscac_pq);
idx1.vscac_vq = get_index(bus_index, idx1.vscac_vq);
idx1.vscdc_pq = get_index(bus_index, idx1.vscdc_pq);
idx1.vscdc_vq = get_index(bus_index, idx1.vscdc_vq);
idx1.pdc = get_index(bus_index, idx1.pdc);


idx = idx1;
Grid_para.pos_ac3(:,1) = get_index(bus_index, Grid_para.pos_ac3(:,1));
Grid_para.pos_ac3(:,2) = get_index(bus_index, Grid_para.pos_ac3(:,2));


Grid_para.pos_dc3(:,1) = get_index(bus_index, Grid_para.pos_dc3(:,1));
Grid_para.pos_dc3(:,2) = get_index(bus_index, Grid_para.pos_dc3(:,2));


%% update voltage in PV nodes
id_PV = get_index(mpc.gen(:,1), setdiff(mpc.gen(:,1),intersect(mpc.areas(:,2), bus_ac)));
E(idx1.pvac) = mpc.gen(id_PV,6);


%% Set the filter parameters
Filter_para.R = 0.008*Y_b; %checked
Filter_para.X = 0.04*Y_b;  %checked
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

tol = 1e-8;
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
    %Make_scatter_plot_small(r,c,Grid_para.n_nodes,mode);
    
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

