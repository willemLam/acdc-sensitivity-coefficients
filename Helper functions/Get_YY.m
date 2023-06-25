function [Grid_para, linedata] = Get_YY(linedata_ac, linedata_dc , Grid_para, Filter_para, idx3, IFlow, idx1)


    V_b = Grid_para.V_b;
    A_b = Grid_para.A_b;
    Y_b = Grid_para.Y_b;
    Vdc_b = Grid_para.Vdc_b;
    Adc_b = Grid_para.Adc_b;
    n_AFE = Grid_para.n_AFE;
    n_ph = Grid_para.n_ph;
    linedata_ac = 'linedata_AC.txt';
    linedata_dc = 'linedata_DC.txt';
    %% Get Linedata ofthe AC and DC grid
    linedata_ac = load (linedata_ac);
    linedata_dc = load (linedata_dc);
    
    %% Make linedata for the virtual nodes
    alp = exp(2*pi/3*1i);
   
     A = 1/3*[1 1     1; 
                 1 alp   alp^2; 
                 1 alp^2 alp];
    ACell =  repmat({A}, 1, n_AFE);
    Atot = blkdiag(ACell{:});
    
    I_b=A_b/(V_b);
    Imag = abs(Atot*IFlow);
    Imag = abs(IFlow);%Imag(2:3:end);
    R_eq_ctu = interp1(Filter_para.IGBT_piecewise(:,1),Filter_para.IGBT_piecewise(:,2),Imag*sqrt(2)*I_b)./Grid_para.V_b./(Imag/sqrt(3)*sqrt(2));
  
    R_eq = (R_eq_ctu*4/pi); 
    R_eq(isnan(R_eq))=0;
    R_eq(isinf(R_eq))=0;
    
    
%     Zloss = (R_eq*Filter_para.Include_losses + (-0.0007 + 0.0020*1i)).*[ones(3,1);zeros(3,1);ones(6,1)];
    Zloss = [0.0314127577579024 + 0.00201268704979953i;0.0314132737817996 + 0.00201259264730807i;0.0314129829837461 + 0.00201218724945477i;2.48699626903600e-05 - 7.99232005513698e-06i;2.46678567656297e-05 - 8.01692419096639e-06i;2.47558925818084e-05 - 7.82762596440029e-06i;0.0299418883701312 + 0.00147479768736150i;0.0299419948526609 + 0.00147481623400933i;0.0299419348549692 + 0.00147469968908606i;0.0301938990969787 + 0.00201491015900307i;0.0301935952104912 + 0.00201498357831528i;0.0301938350902403 + 0.00201519151943597i];
    Zfilter = repmat(complex(Filter_para.R,Filter_para.X),n_AFE*n_ph,1);

    idx_ac3 = [(idx3.vscac_vq);(idx3.vscac_pq)];
    linedata_afe3 = repmat([0 0 0 0 0 1 100 0 1],n_AFE*n_ph,1); 
        for i = 1:n_AFE*n_ph
            linedata_afe3(i,1) = idx_ac3(i)-n_AFE*n_ph;
            linedata_afe3(i,2) = idx_ac3(i);
            linedata_afe3(i,3) = (real(Zloss(i) + Zfilter(i)))./Y_b;
            linedata_afe3(i,4) = (imag(Zloss(i) + Zfilter(i)))./Y_b;
            linedata_afe3(i,5) = (1E-9).*Y_b; % for numeric stability
        end   

    idx_ac = [(idx1.vscac_vq);(idx1.vscac_pq)];
    linedata_afe = repmat([0 0 0 0 0 1 100 0 1],n_AFE,1); 
        for i = 1:n_AFE
            linedata_afe(i,1) = idx_ac(i)-n_AFE;
            linedata_afe(i,2) = idx_ac(i);
            linedata_afe(i,3) = (real(Zloss(i*3) + Zfilter(i*3)))./Y_b;
            linedata_afe(i,4) = (imag(Zloss(i*3) + Zfilter(i*3)))./Y_b;
            linedata_afe(i,5) = (1E-9).*Y_b; % for numeric stability
        end   
    %% Make YY
    
    linedata_ac3 = cell2mat(arrayfun(@(x) x*ones(Grid_para.n_ph,1),linedata_ac,'UniformOutput',false));
    linedata_ac3(:,1:2) = polyphase_indices(linedata_ac(:,1:2),n_ph);
    
    
    linedata_ac3_afe = [linedata_ac3 ; linedata_afe3]   ; 
    [Yac, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~]  = Ymatrix(linedata_ac3_afe,A_b,V_b,[]);

    linedata_ac_afe = [linedata_ac ; linedata_afe]   ; 
    [ ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, linedata2]  = Ymatrix(linedata_ac_afe,A_b,V_b,[]);
    
    [Ydc, ~, ~, ~, ~, ~, ~, ~, ~, ~, linedata2dc]  = Ymatrix(linedata_dc,Adc_b,Vdc_b,[]);


%     Yac = cell2mat(arrayfun(@(x) x*eye(Grid_para.n_ph),Yac,'UniformOutput',false));
    Yac = Yac/3;
    Ydc = Ydc(23:30,23:30)/2;

    YY = blkdiag(Yac,Ydc);
    linedata = [linedata2;linedata2dc];

    Grid_para.G = real(YY);
    Grid_para.B = imag(YY);
    Grid_para.YY = YY;
    Grid_para.Yac = Yac;
    Grid_para.Ydc = Ydc;
end
