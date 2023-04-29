function [Nodal_V_mag,Nodal_V_angle, Nodal_I_mag, Nodal_I_angle, Flow_I_mag, Flow_I_angle, Idc_flow, Idc_inj, Vdc, n_timesteps, M_real, M_imag, Mabs,Nodal_P,Nodal_Q,Pdc_inj,M] = GetEMTPdata_SC(A_b,V_b,Adc_b, Vdc_b,repeat,ZIN,n_phases) 
    addpath '/Users/willem/Documents/phd/State_estimation/Phasor'
    addpath '/Users/willem/Documents/phd/State_estimation/Phasor/biblio'
%     data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/GridV10.3.mat');
%     data =
%     load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/QV_control_Rf0.1_V5.mat');%
%     step in P09 and P23
%     data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/QV_control_Rf0.1_V6_6.mat');

%% balanced
%     data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/GridV7.2all.mat'); 
%     max_range = 10500;
%% unbalanced
    data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/Full_grid_AVG_noloss_balanced_7.mat'); %  balamced
%     data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/Full_grid_AVG_noloss_unbalanced_8.mat'); % light unblance
%     data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/Full_grid_AVG_noloss_unbalanced_9.mat'); % very strong unblance
%     data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/Full_grid_AVG_noloss_unbalanced_10.mat'); % very strong unblance + filter
    
    data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/data_SC_1.mat'); % very strong unblance + filter
    
    max_range = floor(length(data.B01_Va_mag_control)/10)*10;
   
    I_b = A_b/(V_b.*sqrt(3));
    Idc_b = Adc_b/Vdc_b;
    
    range =  1:max_range;
    if n_phases == 1
        n = 3;
    elseif n_phases == 3
        n = 1;
    else
        error('error')
    end
    
    Vac_bus_names = {'B01_Va_mag_control',...
                 'B01_Vb_mag_control',...
                 'B01_Vc_mag_control',...
                 'B02_Va_mag_control',...
                 'B02_Vb_mag_control',...
                 'B02_Vc_mag_control',...
                 'B03_Va_mag_control',...
                 'B03_Vb_mag_control',...
                 'B03_Vc_mag_control',...
                 'B04_Va_mag_control',...
                 'B04_Vb_mag_control',...
                 'B04_Vc_mag_control',...
                 'B05_Va_mag_control',...
                 'B05_Vb_mag_control',...
                 'B05_Vc_mag_control',...
                 'B06_Va_mag_control',...
                 'B06_Vb_mag_control',...
                 'B06_Vc_mag_control',...
                 'B07_Va_mag_control',...
                 'B07_Vb_mag_control',...
                 'B07_Vc_mag_control',...
                 'B08_Va_mag_control',...
                 'B08_Vb_mag_control',...
                 'B08_Vc_mag_control',...
                 'B09_Va_mag_control',...
                 'B09_Vb_mag_control',...
                 'B09_Vc_mag_control',...
                 'B10_Va_mag_control',...
                 'B10_Vb_mag_control',...
                 'B10_Vc_mag_control',...
                 'B11_Va_mag_control',...
                 'B11_Vb_mag_control',...
                 'B11_Vc_mag_control',...
                 'B12_Va_mag_control',...
                 'B12_Vb_mag_control',...
                 'B12_Vc_mag_control',...
                 'B13_Va_mag_control',...
                 'B13_Vb_mag_control',...
                 'B13_Vc_mag_control',...
                 'B14_Va_mag_control',...
                 'B14_Vb_mag_control',...
                 'B14_Vc_mag_control',...
                 'B15_Va_mag_control',...
                 'B15_Vb_mag_control',...
                 'B15_Vc_mag_control',...
                 'B16_Va_mag_control',...
                 'B16_Vb_mag_control',...
                 'B16_Vc_mag_control',...
                 'B17_Va_mag_control',...
                 'B17_Vb_mag_control',...
                 'B17_Vc_mag_control',...
                 'B18_Va_mag_control',...
                 'B18_Vb_mag_control',...
                 'B18_Vc_mag_control'}; 
Nodal_V_mag = [];
for i = 1:n:length(Vac_bus_names)
    Nodal_V_mag(:,(i-1)/n+1) = mean(reshape(data.(Vac_bus_names{i})(range),10,[]));
end
Nodal_V_mag = Nodal_V_mag/(sqrt(2)/sqrt(3)*V_b);
Nodal_V_mag(:,ZIN) = [];

    Vac_bus_names = {'B01_Va_rad_control',...
                 'B01_Vb_rad_control',...
                 'B01_Vc_rad_control',...
                 'B02_Va_rad_control',...
                 'B02_Vb_rad_control',...
                 'B02_Vc_rad_control',...
                 'B03_Va_rad_control',...
                 'B03_Vb_rad_control',...
                 'B03_Vc_rad_control',...
                 'B04_Va_rad_control',...
                 'B04_Vb_rad_control',...
                 'B04_Vc_rad_control',...
                 'B05_Va_rad_control',...
                 'B05_Vb_rad_control',...
                 'B05_Vc_rad_control',...
                 'B06_Va_rad_control',...
                 'B06_Vb_rad_control',...
                 'B06_Vc_rad_control',...
                 'B07_Va_rad_control',...
                 'B07_Vb_rad_control',...
                 'B07_Vc_rad_control',...
                 'B08_Va_rad_control',...
                 'B08_Vb_rad_control',...
                 'B08_Vc_rad_control',...
                 'B09_Va_rad_control',...
                 'B09_Vb_rad_control',...
                 'B09_Vc_rad_control',...
                 'B10_Va_rad_control',...
                 'B10_Vb_rad_control',...
                 'B10_Vc_rad_control',...
                 'B11_Va_rad_control',...
                 'B11_Vb_rad_control',...
                 'B11_Vc_rad_control',...
                 'B12_Va_rad_control',...
                 'B12_Vb_rad_control',...
                 'B12_Vc_rad_control',...
                 'B13_Va_rad_control',...
                 'B13_Vb_rad_control',...
                 'B13_Vc_rad_control',...
                 'B14_Va_rad_control',...
                 'B14_Vb_rad_control',...
                 'B14_Vc_rad_control',...
                 'B15_Va_rad_control',...
                 'B15_Vb_rad_control',...
                 'B15_Vc_rad_control',...
                 'B16_Va_rad_control',...
                 'B16_Vb_rad_control',...
                 'B16_Vc_rad_control',...
                 'B17_Va_rad_control',...
                 'B17_Vb_rad_control',...
                 'B17_Vc_rad_control',...
                 'B18_Va_rad_control',...
                 'B18_Vb_rad_control',...
                 'B18_Vc_rad_control'}; 


Nodal_V_angle = [];
for i = 1:n:length(Vac_bus_names)
    Nodal_V_angle(:,(i-1)/n+1) = mean(reshape(data.(Vac_bus_names{i})(range),10,[])); %data.(Vac_bus_names{i})(range);
end
Nodal_V_angle(:,ZIN) = [];


Iac_bus_names = {'B01_Ia_mag_control',...
                 'B01_Ib_mag_control',...
                 'B01_Ic_mag_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B03_Ia_mag_control',...
                 'B03_Ib_mag_control',...
                 'B03_Ic_mag_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B05_Ia_mag_control',...
                 'B05_Ib_mag_control',...
                 'B05_Ic_mag_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B09_Ia_mag_control',...
                 'B09_Ib_mag_control',...
                 'B09_Ic_mag_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B11_Ia_mag_control',...
                 'B11_Ib_mag_control',...
                 'B11_Ic_mag_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B13_Ia_mag_control',...
                 'B13_Ib_mag_control',...
                 'B13_Ic_mag_control',...
                 'B14_Ia_mag_control',...
                 'B14_Ib_mag_control',...
                 'B14_Ic_mag_control',...
                 'VSI1_Ia_mag_control',...
                 'VSI1_Ib_mag_control',...
                 'VSI1_Ic_mag_control',...
                 'VSI1_Ia_mag_control',...
                 'VSI1_Ib_mag_control',...
                 'VSI1_Ic_mag_control',...
                 'VSI1_Ia_mag_control',...
                 'VSI1_Ib_mag_control',...
                 'VSI1_Ic_mag_control',...
                 'VSI4_Ia_mag_control',...
                 'VSI4_Ib_mag_control',...
                 'VSI4_Ic_mag_control'}; 
Nodal_I_mag = [];
for i = 1:n:length(Iac_bus_names)
    if Iac_bus_names{i} == "no_inj"
        Nodal_I_mag(:,(i-1)/n+1) = zeros(size(Nodal_I_mag,1),1); 
    else
        Nodal_I_mag(:,(i-1)/n+1) = mean(reshape(data.(Iac_bus_names{i})(range),10,[]));%data.(Iac_bus_names{i})(range);
    end
end
Nodal_I_mag = Nodal_I_mag/sqrt(2)/I_b;
Nodal_I_mag(:,ZIN) = [];

Iac_bus_names = {'B01_Ia_rad_control',...
                 'B01_Ib_rad_control',...
                 'B01_Ic_rad_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B03_Ia_rad_control',...
                 'B03_Ib_rad_control',...
                 'B03_Ic_rad_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B05_Ia_rad_control',...
                 'B05_Ib_rad_control',...
                 'B05_Ic_rad_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B09_Ia_rad_control',...
                 'B09_Ib_rad_control',...
                 'B09_Ic_rad_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B11_Ia_rad_control',...
                 'B11_Ib_rad_control',...
                 'B11_Ic_rad_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B13_Ia_rad_control',...
                 'B13_Ib_rad_control',...
                 'B13_Ic_rad_control',...
                 'B14_Ia_rad_control',...
                 'B14_Ib_rad_control',...
                 'B14_Ic_rad_control',...
                 'VSI1_Ia_rad_control',...
                 'VSI1_Ib_rad_control',...
                 'VSI1_Ic_rad_control',...
                 'VSI1_Ia_rad_control',...
                 'VSI1_Ib_rad_control',...
                 'VSI1_Ic_rad_control',...
                 'VSI1_Ia_rad_control',...
                 'VSI1_Ib_rad_control',...
                 'VSI1_Ic_rad_control',...
                 'VSI4_Ia_rad_control',...
                 'VSI4_Ib_rad_control',...
                 'VSI4_Ic_rad_control'}; 
Nodal_I_angle = [];
for i = 1:n:length(Iac_bus_names)
    if Iac_bus_names{i} == "no_inj"
        Nodal_I_angle(:,(i-1)/n+1) = zeros(size(Nodal_I_angle,1),1); 
    else
        Nodal_I_angle(:,(i-1)/n+1) = mean(reshape(data.(Iac_bus_names{i})(range),10,[]));%data.(Iac_bus_names{i})(range);
    end
end
Nodal_I_angle(:,ZIN) = [];

Pac_bus_names = {'B01_P_control',...
                 'B01_P_control',...
                 'B01_P_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B03_P_control',...
                 'B03_P_control',...
                 'B03_P_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B05_P_control',...
                 'B05_P_control',...
                 'B05_P_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B09_P_control',...
                 'B09_P_control',...
                 'B09_P_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B11_P_control',...
                 'B11_P_control',...
                 'B11_P_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B13_P_control',...
                 'B13_P_control',...
                 'B13_P_control',...
                 'B14_P_control',...
                 'B14_P_control',...
                 'B14_P_control',...
                 'B15_P_control',...
                 'B15_P_control',...
                 'B15_P_control',...
                 'B16_P_control',...
                 'B16_P_control',...
                 'B16_P_control',...
                 'B17_P_control',...
                 'B17_P_control',...
                 'B17_P_control',...
                 'B18_P_control',...
                 'B18_P_control',...
                 'B18_P_control'}; 
Nodal_P = [];
for i = 1:n:length(Pac_bus_names)
    if Pac_bus_names{i} == "no_inj"
        Nodal_P(:,(i-1)/n+1) = zeros(size(Nodal_P,1),1); 
    else
        Nodal_P(:,(i-1)/n+1) = mean(reshape(data.(Pac_bus_names{i})(range),10,[]));%data.(Pac_bus_names{i})(range);
    end
end
Nodal_P = Nodal_P/A_b;

Qac_bus_names = {'B01_Q_control',...
                 'B01_Q_control',...
                 'B01_Q_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B03_Q_control',...
                 'B03_Q_control',...
                 'B03_Q_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B05_Q_control',...
                 'B05_Q_control',...
                 'B05_Q_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B09_Q_control',...
                 'B09_Q_control',...
                 'B09_Q_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B11_Q_control',...
                 'B11_Q_control',...
                 'B11_Q_control',...
                 'no_inj',...
                 'no_inj',...
                 'no_inj',...
                 'B13_Q_control',...
                 'B13_Q_control',...
                 'B13_Q_control',...
                 'B14_Q_control',...
                 'B14_Q_control',...
                 'B14_Q_control',...
                 'B15_Q_control',...
                 'B15_Q_control',...
                 'B15_Q_control',...
                 'B16_Q_control',...
                 'B16_Q_control',...
                 'B16_Q_control',...
                 'B17_Q_control',...
                 'B17_Q_control',...
                 'B17_Q_control',...
                 'B18_Q_control',...
                 'B18_Q_control',...
                 'B18_Q_control'}; 
Nodal_Q = [];
for i = 1:n:length(Qac_bus_names)
    if Qac_bus_names{i} == "no_inj"
        Nodal_Q(:,(i-1)/n+1) = zeros(size(Nodal_Q,1),1); 
    else
        Nodal_Q(:,(i-1)/n+1) = mean(reshape(data.(Qac_bus_names{i})(range),10,[]));%data.(Qac_bus_names{i})(range);
    end
end
Nodal_Q = Nodal_Q/A_b;


Iac_flow_name ={'I_B09_B15_Ia_mag_control',...
                'I_B09_B15_Ib_mag_control',...
                'I_B09_B15_Ic_mag_control',...
                'I_B13_B16_Ia_mag_control',...
                'I_B13_B16_Ib_mag_control',...
                'I_B13_B16_Ic_mag_control',...
                'I_B11_B17_Ia_mag_control',...
                'I_B11_B17_Ib_mag_control',...
                'I_B11_B17_Ic_mag_control',...
                'I_B07_B18_Ia_mag_control',...
                'I_B07_B18_Ib_mag_control',...
                'I_B07_B18_Ic_mag_control'};
            
Flow_I_mag = [];
for i = 1:n:length(Iac_flow_name)
    Flow_I_mag(:,(i-1)/n+1) = mean(reshape(data.(Iac_flow_name{i})(range),10,[]));%data.(Iac_flow_name{i})(range);
end
Flow_I_mag = Flow_I_mag/sqrt(2)/I_b;


Iac_flow_name ={'I_B09_B15_Ia_rad_control',...
                'I_B09_B15_Ib_rad_control',...
                'I_B09_B15_Ic_rad_control',...
                'I_B13_B16_Ia_rad_control',...
                'I_B13_B16_Ib_rad_control',...
                'I_B13_B16_Ic_rad_control',...
                'I_B11_B17_Ia_rad_control',...
                'I_B11_B17_Ib_rad_control',...
                'I_B11_B17_Ic_rad_control',...
                'I_B07_B18_Ia_rad_control',...
                'I_B07_B18_Ib_rad_control',...
                'I_B07_B18_Ic_rad_control'};
Flow_I_angle = [];
for i = 1:n:length(Iac_flow_name)
    Flow_I_angle(:,(i-1)/n+1) = mean(reshape(data.(Iac_flow_name{i})(range),10,[]));%data.(Iac_flow_name{i})(range);
end



Vdc_p_bus_name={'B19_Vdc_flow_p_control',...
                'B20_Vdc_flow_p_control',...
                'B21_Vdc_flow_p_control',...
                'B22_Vdc_flow_p_control',...
                'B23_Vdc_p_control',...
                'B24_Vdc_p_control',...
                'B25_Vdc_p_control',...
                'B26_Vdc_p_control'};
Vdc_p = [];
for i = 1:length(Vdc_p_bus_name)
    Vdc_p(:,i) = mean(reshape(data.(Vdc_p_bus_name{i})(range),10,[]));%data.(Vdc_p_bus_name{i})(range);
end
Vdc_p = Vdc_p/Vdc_b;


Vdc_n_bus_name={'B19_Vdc_flow_p_control',...
                'B20_Vdc_flow_p_control',...
                'B21_Vdc_flow_p_control',...
                'B22_Vdc_flow_p_control',...
                'B23_Vdc_p_control',...
                'B24_Vdc_p_control',...
                'B25_Vdc_p_control',...
                'B26_Vdc_p_control'};
Vdc_n = [];
for i = 1:length(Vdc_n_bus_name)
    Vdc_n(:,i) = mean(reshape(data.(Vdc_n_bus_name{i})(range),10,[]));%data.(Vdc_n_bus_name{i})(range);
end
Vdc_n = Vdc_n/Vdc_b;

Idc_p_bus_name={'B19_Idc_flow_p_control',...
                'B20_Idc_flow_p_control',...
                'B21_Idc_flow_p_control',...
                'B22_Idc_flow_p_control',...
                'B23_Idc_inj_p_control',...
                'B24_Idc_inj_p_control',...
                'B25_Idc_inj_p_control',...
                'B26_Idc_inj_p_control'};
Idc_p = [];
for i = 1:length(Idc_p_bus_name)
    Idc_p(:,i) = mean(reshape(data.(Idc_p_bus_name{i})(range),10,[]));%data.(Idc_p_bus_name{i})(range);
end
Idc_p = Idc_p/Idc_b;


Idc_n_bus_name={'B19_Idc_flow_p_control',...
                'B20_Idc_flow_p_control',...
                'B21_Idc_flow_p_control',...
                'B22_Idc_flow_p_control',...
                'B23_Idc_inj_p_control',...
                'B24_Idc_inj_p_control',...
                'B25_Idc_inj_p_control',...
                'B26_Idc_inj_p_control'};
Idc_n = [];
for i = 1:length(Idc_n_bus_name)
    Idc_n(:,i) = mean(reshape(data.(Idc_n_bus_name{i})(range),10,[]));%data.(Idc_n_bus_name{i})(range);
end
Idc_n = Idc_n/Idc_b;

Pdc_p_bus_name={'B19_P_p_control',...
                'B20_P_p_control',...
                'B21_P_p_control',...
                'B22_P_p_control',...
                'B23_P_p_control',...
                'B24_P_p_control',...
                'B25_P_p_control',...
                'B26_P_p_control'};
Pdc_p = [];
for i = 1:length(Pdc_p_bus_name)
    Pdc_p(:,i) = mean(reshape(data.(Pdc_p_bus_name{i})(range),10,[]));%data.(Pdc_p_bus_name{i})(range);
end
Pdc_p = Pdc_p/Adc_b;

Pdc_n_bus_name={'B19_P_n_control',...
                'B20_P_n_control',...
                'B21_P_n_control',...
                'B22_P_n_control',...
                'B23_P_n_control',...
                'B24_P_n_control',...
                'B25_P_n_control',...
                'B26_P_n_control'};
Pdc_n = [];
for i = 1:length(Pdc_n_bus_name)
    Pdc_n(:,i) = mean(reshape(data.(Pdc_n_bus_name{i})(range),10,[]));%data.(Pdc_n_bus_name{i})(range);
end
Pdc_n = Pdc_n/Adc_b;


Vdc = Vdc_n + Vdc_p;
Idc_inj = (Idc_n + Idc_p);
Pdc_inj = (Pdc_n + Pdc_p);
Idc_flow = Idc_inj(:,1:4);

M_name = {  'VSI1_Mreal_a_2_control',...
            'VSI1_Mreal_b_2_control',...
            'VSI1_Mreal_c_2_control',...
            'VSI1_Mreal_a_2_control',...
            'VSI1_Mreal_b_2_control',...
            'VSI1_Mreal_c_2_control',...
            'VSI1_Mreal_a_2_control',...
            'VSI1_Mreal_b_2_control',...
            'VSI1_Mreal_c_2_control',...
            'VSI4_Mreal_a_2_control',...
            'VSI4_Mreal_b_2_control',...
            'VSI4_Mreal_c_2_control'};
M_real = [];
for i = 1:n:length(M_name)
    M_real(:,(i-1)/n+1) = mean(reshape(data.(M_name{i})(range),10,[]));%data.(M_name{i})(range);
end
M_real = M_real*(Vdc_b/2)/(sqrt(2)/sqrt(3)*V_b)*2/2;%(Vdc_b/2)/(sqrt(2)/sqrt(3)*V_b)*2/2;



M_name = {  'VSI1_Mimag_a_2_control',...
            'VSI1_Mimag_b_2_control',...
            'VSI1_Mimag_c_2_control',...
            'VSI1_Mimag_a_2_control',...
            'VSI1_Mimag_b_2_control',...
            'VSI1_Mimag_c_2_control',...
            'VSI1_Mimag_a_2_control',...
            'VSI1_Mimag_b_2_control',...
            'VSI1_Mimag_c_2_control',...
            'VSI4_Mimag_a_2_control',...
            'VSI4_Mimag_b_2_control',...
            'VSI4_Mimag_c_2_control'};

M_imag = [];
for i = 1:n:length(M_name)
    M_imag(:,(i-1)/n+1) = mean(reshape(data.(M_name{i})(range),10,[]));%data.(M_name{i})(range);
end
M_imag = M_imag*(Vdc_b/2)/(sqrt(2)/sqrt(3)*V_b)*2/2;%(Vdc_b/2)/(sqrt(2)/sqrt(3)*V_b)*2/2;

M = complex(M_real,M_imag);
% M = 2*conj(M./abs(M).^2)/2;
M(isnan(M))=0;

Mabs = 0;
    
   
    Nodal_V_mag = repmat(Nodal_V_mag,repeat,1);
    Nodal_V_angle = repmat(Nodal_V_angle,repeat,1);
    Nodal_I_mag = repmat(Nodal_I_mag,repeat,1);
    Nodal_I_angle = repmat(Nodal_I_angle,repeat,1);
    Nodal_P = repmat(Nodal_P,repeat,1);
    Nodal_Q = repmat(Nodal_Q,repeat,1);
    Flow_I_mag = repmat(Flow_I_mag,repeat,1);
    Flow_I_angle = repmat(Flow_I_angle,repeat,1);
    Vdc = repmat(Vdc,repeat,1);
    Idc_flow = repmat(Idc_flow,repeat,1);
    Idc_inj = repmat(Idc_inj,repeat,1);
    M_real = repmat(M_real,repeat,1);
    M_imag = repmat(M_imag,repeat,1);
    Pdc_inj = repmat(Pdc_inj,repeat,1);
%     Nodal_V_mag = make_nice_step(Nodal_V_mag,repeat);
%     Nodal_V_angle = make_nice_step(Nodal_V_angle,repeat);
%     Nodal_I_mag = make_nice_step(Nodal_I_mag,repeat);
%     Nodal_I_angle = make_nice_step(Nodal_I_angle,repeat);
%     Flow_I_mag = make_nice_step(Flow_I_mag,repeat);
%     Flow_I_angle = make_nice_step(Flow_I_angle,repeat);
%     Vdc = make_nice_step(Vdc,repeat);
%     Idc_flow = make_nice_step(Idc_flow,repeat);
%     Idc_inj = make_nice_step(Idc_inj,repeat);
%     Mreal = make_nice_step(Mreal,repeat);
%     Mimag = make_nice_step(Mimag,repeat);
    
%     Nodal_V_mag = mean(Nodal_V_mag(:,:)).*ones(size(Nodal_V_mag));
%     Nodal_V_angle = mean(Nodal_V_angle(:,:)).*ones(size(Nodal_V_angle));
%     Nodal_I_angle = mean(Nodal_I_angle(:,:)).*ones(size(Nodal_I_angle));
%     Nodal_I_mag = mean(Nodal_I_mag(:,:)).*ones(size(Nodal_I_mag));
%     Nodal_P = mean(Nodal_P(:,:)).*ones(size(Nodal_P));
%     Nodal_Q = mean(Nodal_Q(:,:)).*ones(size(Nodal_Q));
%     Flow_I_mag = mean(Flow_I_mag(:,:)).*ones(size(Flow_I_mag));
%     Flow_I_angle = mean(Flow_I_angle(:,:)).*ones(size(Flow_I_angle));
%     Vdc = mean(Vdc(:,:)).*ones(size(Vdc));
%     Idc_flow = mean(Idc_flow(:,:)).*ones(size(Idc_flow));
%     Idc_inj = mean(Idc_inj(:,:)).*ones(size(Idc_inj));
%     M_real = mean(M_real(:,:)).*ones(size(M_real));
%     M_imag = mean(M_imag(:,:)).*ones(size(M_imag));
%     Pdc_inj = mean(Pdc_inj(:,:)).*ones(size(Pdc_inj));
    

% Nodal_V_mag = [Nodal_V_mag(1,:).*ones(500,1) ; Nodal_V_mag; Nodal_V_mag(end,:).*ones(500,1)];
% Nodal_V_angle = [Nodal_V_angle(1,:).*ones(500,1) ; Nodal_V_angle; Nodal_V_angle(end,:).*ones(500,1)];
% Nodal_I_angle = [Nodal_I_angle(1,:).*ones(500,1) ; Nodal_I_angle; Nodal_I_angle(end,:).*ones(500,1)];
% Nodal_I_mag = [Nodal_I_mag(1,:).*ones(500,1) ; Nodal_I_mag; Nodal_I_mag(end,:).*ones(500,1)];
% Nodal_P = [Nodal_P(1,:).*ones(500,1) ; Nodal_P; Nodal_P(end,:).*ones(500,1)];
% Nodal_Q = [Nodal_Q(1,:).*ones(500,1) ; Nodal_Q; Nodal_Q(end,:).*ones(500,1)];
% Flow_I_mag = [Flow_I_mag(1,:).*ones(500,1) ; Flow_I_mag; Flow_I_mag(end,:).*ones(500,1)];
% Flow_I_angle = [Flow_I_angle(1,:).*ones(500,1) ; Flow_I_angle; Flow_I_angle(end,:).*ones(500,1)];
% Vdc = [Vdc(1,:).*ones(500,1) ; Vdc; Vdc(end,:).*ones(500,1)];
% Idc_flow = [Idc_flow(1,:).*ones(500,1) ; Idc_flow; Idc_flow(end,:).*ones(500,1)];
% Idc_inj = [Idc_inj(1,:).*ones(500,1) ; Idc_inj; Idc_inj(end,:).*ones(500,1)];
% M_real = [M_real(1,:).*ones(500,1) ; M_real; M_real(end,:).*ones(500,1)];
% M_imag = [M_imag(1,:).*ones(500,1) ; M_imag; M_imag(end,:).*ones(500,1)];
% Pdc_inj = [Pdc_inj(1,:).*ones(500,1) ; Pdc_inj; Pdc_inj(end,:).*ones(500,1)];


    n_timesteps = max(size(Idc_flow)) ;
end


function data = make_nice_step(data,repeat)
    data_r = mean(data(1:round(2*end/5),:)).*ones(size(data(1:round(2*end/5),:)));
    data_r = repmat(data_r,3*repeat/4,1);
    
    data_l = mean(data(round(3*end/5):end,:)).*ones(size(data(round(3*end/5):end,:)));
    data_l = repmat(data_l,repeat/4,1);
    
    data_m = [data_r(end-24:end,:) ; data_l(end-24:end,:)];
%     data_m = repmat(data_m,10,1);
%     data = [data_r ; data(2*end/5:3*end/5,:) ; data_l];
    data = [data_r ; data_m; data_l];
end


% data.('VSI1_Mreal_a_2_control')(range).*data.('I_B09_B15_Ia_mag_control')(range).*cos(data.('I_B09_B15_Ia_rad_control')(range)) + data.('VSI1_Mimag_a_2_control')(range).*data.('I_B09_B15_Ia_mag_control')(range).*sin(data.('I_B09_B15_Ia_rad_control')(range))
% data.('B19_Idc_flow_p_control')(range)
% 
% data.('VSI1_Mimag_a_2_control')(range).*data.('I_B09_B15_Ia_mag_control')(range).*cos(data.('I_B09_B15_Ia_rad_control')(range)) - data.('VSI1_Mreal_a_2_control')(range).*data.('I_B09_B15_Ia_mag_control')(range).*sin(data.('I_B09_B15_Ia_rad_control')(range))
% 0
% 
% M = sqrt(data.('VSI1_Mreal_a_2_control')(range).^2 + data.('VSI1_Mimag_a_2_control')(range).^2)*(Vdc_b/2)/(sqrt(2)/sqrt(3)*V_b)*2/2
% 
% data.('B19_P_p_control')(range)*2
% data.('B15_P_control')(range)
% 
% sqrt(data.('B15_P_control')(range).^2 + data.('B15_Q_control')(range).^2)
% data.('B15_Va_mag_control')(range).*data.('I_B09_B15_Ia_mag_control')(range)
% 
% 
% data.('B19_Idc_flow_p_control')(range)*2
% data.('I_B09_B15_Ia_mag_control')(range)*(sqrt(2)/sqrt(3))
% 
% data.('B19_Vdc_flow_p_control')(range)
% data.('B15_Va_mag_control')(range)/(sqrt(2)/sqrt(3))
% 
% V = data.('B15_Va_mag_control')(range).*exp(data.('B15_Va_rad_control')(range)*1i)
% I = data.('I_B09_B15_Ia_mag_control')(range).*exp(data.('I_B09_B15_Ia_rad_control')(range)*1i)
% conj(V).*I
% 
% data.('B19_Idc_flow_p_control')(range).*data.('B19_Vdc_flow_p_control')(range)*2
% 
% 
% 
% data.('VSI1_vac2_r_reconstructed')(range)