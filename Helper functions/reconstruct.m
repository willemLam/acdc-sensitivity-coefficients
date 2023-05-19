function V = reconstruct(Y,N_kron,V_kron,n_phases)
% V = reconstruct(Y,N_kron,V_kron)
%
% INPUT
% - Y       Original admittance matrix.
% - N_kron  Nodes which were not reduced.
% - V_kron  Corresponding phase-to-ground voltage phasors.
% - f       Fundamental frequency.
%
% OUTPUT
% - V       Complete phase-to-ground voltage phasors.

%% Check

if(~isa(Y,'numeric'))
    error('type of input "Y"');
elseif(~isa(N_kron,'numeric'))
    error('type of input "N_kron"');
elseif(~isa(V_kron,'numeric'))
    error('type of input "V_kron"');
elseif(size(Y,1)~=size(Y,2))
    error('size of input "Y"');
elseif(size(N_kron,2)~=1)
    error('size of input "N_kron"');
elseif(size(V_kron,1)~=n_phases*size(N_kron,1))
    error('size of input "V_kron"');
elseif(~all(mod(N_kron,1)==0))
    error('value of input "N_kron"');
elseif(~all(ismember(N_kron,1:size(Y,1)))) % does this make sense? Y is in 3 phases as well
    error('inconsistency of input "Y" and "N_kron"');
end

%% Calculate
n_nodes = size(Y,1)/n_phases;
n_frequencies = size(Y,3);

idx_given = N_kron;
idx_wanted = setdiff((1:n_nodes)',idx_given);

Y = mat2cell(Y,n_phases*ones(n_nodes,1),n_phases*ones(n_nodes,1),n_frequencies);

Y_22 = cell2mat(Y(idx_wanted,idx_wanted,:));
Y_21 = cell2mat(Y(idx_wanted,idx_given,:));

V_1 = V_kron;
V_2 = zeros(size(Y_22,1),size(V_1,2));
for h = 1 : size(V_kron,2)
    V_2(:,h) = - (Y_22(:,:,h) \ (Y_21(:,:,h)) * V_1(:,h));
end

V = zeros(n_nodes*n_phases,n_frequencies);
V = mat2cell(V,n_phases*ones(n_nodes,1),n_frequencies);

V_1 = mat2cell(V_1,n_phases*ones(length(idx_given),1),n_frequencies);
V_2 = mat2cell(V_2,n_phases*ones(length(idx_wanted),1),n_frequencies);

V(idx_given) = V_1;
V(idx_wanted) = V_2;

V = cell2mat(V);
end