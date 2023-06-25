function [YY, YYL, YL, YT, YYT, I_b, Ampacities, y_ih, y_i, A, linedata]  = make_Y_augmented(Text,Grid_para,idx1,Zloss,Zfilter,type)

    V_b = Grid_para.V_b;
    A_b = Grid_para.A_b;
    n_ac = Grid_para.n_ac;
    n_dc = Grid_para.n_dc;
    linedata = load(Text);


    I_b=A_b/(V_b.*sqrt(3));
    Yb=A_b/V_b^2;

    %%%The line parameters are expressed here in p.u.
    line_lengths=linedata(:,6);
    linedata(:,3)=(linedata(:,3).*Yb).*line_lengths;
    linedata(:,4)=(linedata(:,4).*Yb).*line_lengths;
    linedata(:,5)=(linedata(:,5)./Yb).*line_lengths;


    %% give nodes a new number
    idx_to_augment = sort([idx1.vscac_pq; idx1.vscac_vq]);
% 
%     for i = 1:length(idx_to_change)
%         idx_tmp = find(linedata(:,2)==idx_to_change(i));
%         linedata(idx_tmp,2) = n_ac + i;
%     end
    if type == 0
    %% add the IGBT losses
        linedata_template = [0 0 0 0 0 1 100 0 1];
        for i = 1:length(idx_to_augment)
            linedata_append = linedata_template;
            linedata_append(1) = idx_to_augment(i);
            linedata_append(2) = n_ac + i;
            linedata_append(3) = real(Zloss(i) + Zfilter(i));
            linedata_append(4) = imag(Zloss(i) + Zfilter(i));
            linedata_append(5) = 1E-9; % for numeric stability

            linedata = [linedata ; linedata_append];
        end   
    
    elseif type == 1
        
    %% add the IGBT losses
        linedata_template = [0 0 0 0 0 1 100 0 1];
        for i = 1:length(idx_to_augment)
            linedata_append = linedata_template;
            linedata_append(1) = idx_to_augment(i);
            linedata_append(2) = n_ac + i;
            linedata_append(3) = real(Zloss(i));
            linedata_append(4) = imag(Zloss(i));
            linedata_append(5) = 1E-9; % for numeric stability

            linedata = [linedata ; linedata_append];
        end   

        %% add the filter
        linedata_template = [0 0 0 0 0 1 100 0 1];
        for i = 1:length(idx_to_augment)
            linedata_append = linedata_template;
            linedata_append(1) = n_ac + i;
            linedata_append(2) = n_ac + length(idx_to_augment) + i;
            linedata_append(3) = real(Zfilter(i));
            linedata_append(4) = imag(Zfilter(i));
            linedata_append(5) = 1E-9; % for numeric stability

            linedata = [linedata ; linedata_append];
        end    

    end
    
    n_lines=length(linedata(:,1));                          %%%number of lines
    n_nodes=max(max(linedata(:,1)),max(linedata(:,2))); %%%number of nodes

    %% Third step: Building the primitive branch and shunt admittance matrices

    % disp('STEP 3: Building the network primitive branch and shunt admittance matrices....')
    % disp('                                                            ')

    y_i_ih=zeros(n_nodes,n_nodes);
    y_i=zeros(1,n_nodes); 

    y_ih=1./(linedata(:,3)+1i*linedata(:,4)); %%%the longitudinal admittance of each line
    absy_ih=abs(y_ih);
    angley_ih=angle(y_ih);
    g_ih=real(y_ih);
    b_ih=imag(y_ih);


    for k=1:n_lines
        y_i_ih(linedata(k,1),linedata(k,2))=linedata(k,5)*1i/2;
        y_i_ih(linedata(k,2),linedata(k,1))=linedata(k,5)*1i/2;
    end

    for k=1:n_nodes
        y_i(k)=sum(y_i_ih(k,:)); %%%the transversal admittance elements
    end

    g_i_ih=real(y_i_ih);
    b_i_ih=imag(y_i_ih);

    % Extra
    for k=1:n_lines
        % This is not 100% Correct, I am assuming that both shunt elements are
        % equal (reasonable...)
        y_i_cur(k)=y_i(linedata(k,1)) ;%+ y_i(linedata(k,2)); 
    end

    YL = diag(y_ih);
    YT = diag(y_i);
    YT_cur = diag(y_i_cur);

    %----------------------------------------------------------------------------------------------%
    %% Fourth step: compute incidence matrix of the network

    % disp('STEP 4: Building the network incidence matrix....')
    % disp('                                                 ')


    %%% The branch-to-node incidence matrix A is computed here
    %%% A is of size (number of branches, number of nodes)
    %%% A_ij=0 if branch i is not connected to node j
    %%% A_ij=1 if current in branch i is directed away from node j
    %%% A_ij=-1 if current in branch i is directed towards node j

    A=zeros(n_lines,n_nodes); 

    for k=1:n_lines
        A(k,linedata(k,1))=1;
        A(k,linedata(k,2))=-1;
    end


    %% Fifth step: compute [Y] matrix and visualize network connectivity

    % disp('STEP 5: Computing the [Y] matrix....')
    % disp('                                                 ')

    Y=A.'*YL*A + YT;

    % For Current Sensitivites
    YY = Y;
    YYL = -(YY - diag(diag(YY)));
    YYT = -(A.'*YT_cur*A - diag(diag(A.'*YT_cur*A)));
    Ampacities = (-(A.'*diag(linedata(:,7))*A - diag(diag(A.'*diag(linedata(:,7))*A))))/I_b;

    %% kron reduce
% 
%     N_kron = 1:n_ac;
%     YY_prime = reduce(YY,N_kron',1);
    
   
end