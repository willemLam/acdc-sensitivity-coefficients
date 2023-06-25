
A_b = 1e5;
V_b= 400;

Y_b = A_b/V_b^2; 
I_b=A_b/(V_b*sqrt(3));


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



data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/data_SC_test_6.mat'); % balanced + filter




V1 = transpose(complex(data.B2_Va_mag_control.*cos(data.B2_Va_rad_control), data.B2_Va_mag_control.*sin(data.B2_Va_rad_control)))./V_b;
V2 = transpose(complex(data.VSI1_Bp_Va_mag_control.*cos(data.VSI1_Bp_Va_rad_control), data.VSI1_Bp_Va_mag_control.*sin(data.VSI1_Bp_Va_rad_control)))./V_b;
V1p = transpose(complex(data.VSI1_B_Va_mag_control.*cos(data.VSI1_B_Va_rad_control), data.VSI1_B_Va_mag_control.*sin(data.VSI1_B_Va_rad_control)))./V_b;
V1pp = transpose(complex(data.VSI1_Bpp_Va_mag_control.*cos(data.VSI1_Bpp_Va_rad_control), data.VSI1_Bpp_Va_mag_control.*sin(data.VSI1_Bpp_Va_rad_control)))./V_b;

R = data.VSI1_R_control
L = data.VSI1_L_control

R = R(40:end);
L = L(40:end);



I1 = transpose(complex(data.VSI1_Ia_mag_control.*cos(data.VSI1_Ia_rad_control), data.VSI1_Ia_mag_control.*sin(data.VSI1_Ia_rad_control)))./I_b;

V1 = V1(40:end);
V1p = V1p(40:end);
V1pp = V1pp(40:end);
V2 = V2(40:end);
I1 = I1(40:end);

 Imag = mean(abs(I1))
    R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b/sqrt(3))./V_b./Imag;
    R_eq = (R_eq_ctu*4/pi); 
    R_eq(isnan(R_eq))=0;
    
    
mean((V2 - V1)./I1)
mean((V2 - V1p)./I1)

mean(complex(R,L*2*pi*50))*Y_b + R_eq
complex(Filter_para.R + R_eq,Filter_para.X)

%% full system


data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/data_SC_4.mat'); % balanced + filter



V1 = transpose(complex(data.B15_Va_mag_control.*cos(data.B15_Va_rad_control), data.B15_Va_mag_control.*sin(data.B15_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V1ppp = transpose(complex(data.VSI1_Bppp_Va_mag_control.*cos(data.VSI1_Bppp_Va_rad_control), data.VSI1_Bppp_Va_mag_control.*sin(data.VSI1_Bppp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V1pp = transpose(complex(data.VSI1_Bpp_Va_mag_control.*cos(data.VSI1_Bpp_Va_rad_control), data.VSI1_Bpp_Va_mag_control.*sin(data.VSI1_Bpp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V1p = transpose(complex(data.VSI1_Bp_Va_mag_control.*cos(data.VSI1_Bp_Va_rad_control), data.VSI1_Bp_Va_mag_control.*sin(data.VSI1_Bp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
I1 = transpose(complex(data.VSI1_Ia_mag_control.*cos(data.VSI1_Ia_rad_control), data.VSI1_Ia_mag_control.*sin(data.VSI1_Ia_rad_control)))/sqrt(2)./I_b;

V1 = V1(1000:end);
V1ppp = V1ppp(1000:end);
V1pp = V1pp(1000:end);
V1p = V1p(1000:end);
I1 = I1(1000:end);
R_eq1 = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),abs(I1)*sqrt(2)*I_b)./V_b./(abs(I1)/sqrt(3)*sqrt(2))*4/pi;



V3 = transpose(complex(data.B17_Va_mag_control.*cos(data.B17_Va_rad_control), data.B17_Va_mag_control.*sin(data.B17_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V3ppp = transpose(complex(data.VSI3_Bppp_Va_mag_control.*cos(data.VSI3_Bppp_Va_rad_control), data.VSI3_Bppp_Va_mag_control.*sin(data.VSI3_Bppp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V3pp = transpose(complex(data.VSI3_Bpp_Va_mag_control.*cos(data.VSI3_Bpp_Va_rad_control), data.VSI3_Bpp_Va_mag_control.*sin(data.VSI3_Bpp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
V3p = transpose(complex(data.VSI3_Bp_Va_mag_control.*cos(data.VSI3_Bp_Va_rad_control), data.VSI3_Bp_Va_mag_control.*sin(data.VSI3_Bp_Va_rad_control)))/(sqrt(2)/sqrt(3)*V_b);
I3 = transpose(complex(data.VSI3_Ia_mag_control.*cos(data.VSI3_Ia_rad_control), data.VSI3_Ia_mag_control.*sin(data.VSI3_Ia_rad_control)))/sqrt(2)./I_b;

V3 = V3(1000:end);
V3ppp = V3ppp(1000:end);
V3pp = V3pp(1000:end);
V3p = V3p(1000:end);
I3 = I3(1000:end);
R_eq3 = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),abs(I3)*sqrt(2)*I_b)./V_b./(abs(I3)/sqrt(3)*sqrt(2))*4/pi;


% 
% V1 = transpose(complex(data.B17_Va_mag_control.*cos(data.B17_Va_rad_control), data.B17_Va_mag_control.*sin(data.B17_Va_rad_control)))./V_b;
% V2 = transpose(complex(data.VSI3_Bp_Va_mag_control.*cos(data.VSI3_Bp_Va_rad_control), data.VSI3_Bp_Va_mag_control.*sin(data.VSI3_Bp_Va_rad_control)))./V_b;
% I = transpose(complex(data.VSI3_Ia_mag_control.*cos(data.VSI3_Ia_rad_control), data.VSI3_Ia_mag_control.*sin(data.VSI3_Ia_rad_control)))./I_b;
% 
% 
% V1 = transpose(complex(data.B18_Va_mag_control.*cos(data.B18_Va_rad_control), data.B18_Va_mag_control.*sin(data.B18_Va_rad_control)))./V_b;
% V2 = transpose(complex(data.VSI4_Bp_Va_mag_control.*cos(data.VSI4_Bp_Va_rad_control), data.VSI4_Bp_Va_mag_control.*sin(data.VSI4_Bp_Va_rad_control)))./V_b;
% I = transpose(complex(data.VSI4_Ia_mag_control.*cos(data.VSI4_Ia_rad_control), data.VSI4_Ia_mag_control.*sin(data.VSI4_Ia_rad_control)))./I_b;





    
(V1(1) - V1ppp(1))./I1(1)
(V1(1) - V1pp(1))./I1(1)
(V1pp(1) - V1p(1))./I1(1)

%% 1
figure
subplot(2,1,1)
hold on
plot(-real(((V1pp - V1p)./I1)))
plot(real(R_eq1))
subplot(2,1,2)
hold on
plot(imag(((V1pp - V1p)./I1)))
plot(imag(R_eq1))

%% 3
figure
subplot(2,1,1)
hold on
plot(-real(((V3pp - V3p)./I3)))
plot(real(R_eq3))
subplot(2,1,2)
hold on
plot(imag(((V3pp - V3p)./I3)))
plot(imag(R_eq3))
