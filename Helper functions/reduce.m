function Y_kron = reduce(Y,N_kron,n_phases)
% Y_kron = reduce(Y,N_kron) performs Kron reduction.
%
% INPUT
% - Y       Original admittance matrix.
% - N_kron  Nodes remaining after the reduction.
%
% OUTPUT
% - Y_kron  Reduced admittance matrix.

%% Check

if(~isa(Y,'numeric'))
    error('type of input "Y"');
elseif(~isa(N_kron,'numeric'))
    error('type of input "N_kron"');
elseif(size(Y,1)~=size(Y,2))
    error('size of input "Y"');
elseif(size(N_kron,2)~=1)
    error('size of input "N_kron"');
elseif(~all(mod(N_kron,1)==0))
    error('value of input "N_kron"');
elseif(~all(ismember(N_kron,1:size(Y,1))))
    error('inconsistency of input "Y" and "N_kron"');
end

%% Calculate
n_nodes = size(Y,1)/n_phases;
n_frequencies = size(Y,3);

idx_consider = N_kron;
idx_remove = setdiff((1:n_nodes)',idx_consider);

Y = mat2cell(Y,n_phases*ones(n_nodes,1),n_phases*ones(n_nodes,1),n_frequencies);

Y_11 = cell2mat(Y(idx_consider,idx_consider,:));
Y_12 = cell2mat(Y(idx_consider,idx_remove,:));
Y_21 = cell2mat(Y(idx_remove,idx_consider,:));
Y_22 = cell2mat(Y(idx_remove,idx_remove,:));

if isempty(idx_remove)
    Y_kron = Y_11;
else
    Y_kron = zeros(size(Y_11));
    for h = 1:size(Y_kron,3)
        Y_kron(:,:,h) = Y_11(:,:,h) - Y_12(:,:,h) * (Y_22(:,:,h) \ Y_21(:,:,h)); % Schur complement
    end
end

end