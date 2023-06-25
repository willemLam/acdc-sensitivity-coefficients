
A_b = 1e3;
V_b= 10;
Y_b = A_b/V_b^2; 
I_b=A_b/(V_b);


Filter_para.R = 0.04*Y_b; %checked %0.08
Filter_para.X = 0.04*Y_b;  %checked

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



data = load('/Users/willem/Documents/phd/State_estimation/EMTP/Experiments/TEST_Loss3.mat'); % balanced + filter

V = transpose(complex(data.V_mag_control.*cos(data.V_rad_control), data.V_mag_control.*sin(data.V_rad_control)))./V_b;
I = transpose(complex(data.I_mag_control.*cos(data.I_rad_control), data.I_mag_control.*sin(data.I_rad_control)))./I_b;

V = V(2:end);
I = I(2:end);



    alp = exp(2*pi/3*1i);
    A = 1/3*[1 1     1; 
                 1 alp   alp^2; 
                 1 alp^2 alp];
    ACell =  repmat({A}, 1, n_AFE);
    Atot = blkdiag(ACell{:});
    
    
    Imag = abs(I);

    R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*I_b)./V_b./Imag;
    R_eq = (R_eq_ctu*4/pi); 
    R_eq(isnan(R_eq))=0;
    R_eq(isinf(R_eq))=0;
    
    Zloss = R_eq*Filter_para.Include_losses;
    Zfilter = repmat(complex(Filter_para.R,Filter_para.X),n_AFE,1);

    
    
    Zloss_model = Zloss + complex(Filter_para.R,Filter_para.X);
    Zloss_emtp = V./I;
    
    
    figure
    subplot(2,1,1)
    hold on
    plot(real(Zloss_model.*I),abs(I).*sign(angle(I)))
    plot(real(Zloss_emtp.*I),abs(I).*sign(angle(I)))
    xlabel('V loss')
    ylabel('Current')
    
    subplot(2,1,2)
    hold on
    plot(real(((Zloss_model - Zloss_emtp)./Zloss_emtp).*I),abs(I).*sign(angle(I)))
    xlabel('Delta V loss')
    ylabel('Current')
    
    