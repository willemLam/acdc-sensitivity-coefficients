function [Kp, Kq, Kv] = transform_K_complex(K,Grid_para, idx)

    J_PR = zeros(Grid_para.n_nodes);
    J_QR = zeros(Grid_para.n_nodes);
    J_VR = zeros(Grid_para.n_nodes);
    
    for k = 1:size(K,1)
        if( sum( K{k,1} == idx.slack))
            continue
        elseif( sum( K{k,1} == idx.pqac))
            J_PR(:,k) = (K{k,2}{1,1});
            J_QR(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.pvac ))
            J_PR(:,k) = (K{k,2}{1,1});
            J_VR(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.pdc ) )
            J_PR(:,k) = (K{k,2}{1,1});
         
        elseif( sum( K{k,1} == idx.vdc ) )
            J_VR(:,k) = (K{k,2}{1,1});

        elseif( sum( K{k,1} == idx.vscac_pq))
            J_PR(:,k) = (K{k,2}{1,1});    
            J_QR(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.vscac_vq))
            J_PR(:,k) = (K{k,2}{1,1});
            J_QR(:,k) = (K{k,2}{2,1});

        elseif( sum( K{k,1} == idx.vscdc_pq ))
            J_PR(:,k) = (K{k,2}{1,1});

        elseif( sum( K{k,1} == idx.vscdc_vq ))
            J_VR(:,k) = (K{k,2}{1,1});
        else
            warning('somethings off mate')
        end
    end

    Kp = J_PR;%(4:end,4:end);
    Kq = J_QR;%(4:end,4:end);
    Kv = J_VR;%(4:end,4:end);
end