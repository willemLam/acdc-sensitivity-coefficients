function [Kp_abs, Kq_abs, Kv_abs,Kp_ang, Kq_ang, Kv_ang] = transform_K(K,Grid_para_augmented, idx3_augmented)

    J_PR_abs = zeros(Grid_para_augmented.n_nodes);
    J_QR_abs = zeros(Grid_para_augmented.n_nodes);
    J_VR_abs = zeros(Grid_para_augmented.n_nodes);
    J_PR_ang = zeros(Grid_para_augmented.n_nodes);
    J_QR_ang = zeros(Grid_para_augmented.n_nodes);
    J_VR_ang = zeros(Grid_para_augmented.n_nodes);
    for k = 1:size(K,1)
        if( sum( K{k,1} == idx3_augmented.slack))
            continue
        elseif( sum( K{k,1} == idx3_augmented.pqac))
            J_PR_abs(:,k) = real(K{k,5}{1,1});
            J_QR_abs(:,k) = real(K{k,5}{2,1});
            J_PR_ang(:,k) = real(K{k,6}{1,1});
            J_QR_ang(:,k) = real(K{k,6}{2,1});

        elseif( sum( K{k,1} == idx3_augmented.pvac ))
            J_PR_abs(:,k) = real(K{k,5}{1,1});
            J_VR_abs(:,k) = real(K{k,5}{2,1});
            J_PR_ang(:,k) = real(K{k,6}{1,1});
            J_VR_ang(:,k) = real(K{k,6}{2,1});

        elseif( sum( K{k,1} == idx3_augmented.pdc ) )
            J_PR_abs(:,k) = real(K{k,5}{1,1});
            J_PR_ang(:,k) = real(K{k,6}{1,1});

        elseif( sum( K{k,1} == idx3_augmented.vdc ) )
            J_VR_abs(:,k) = real(K{k,5}{1,1});
            J_VR_ang(:,k) = real(K{k,6}{1,1});

        elseif( sum( K{k,1} == idx3_augmented.vscac_pq))
            J_PR_abs(:,k) = real(K{k,5}{1,1});    
            J_QR_abs(:,k) = real(K{k,5}{2,1});
            J_PR_ang(:,k) = real(K{k,6}{1,1});    
            J_QR_ang(:,k) = real(K{k,6}{2,1});

        elseif( sum( K{k,1} == idx3_augmented.vscac_vq))
            J_PR_abs(:,k) = real(K{k,5}{1,1});
            J_QR_abs(:,k) = real(K{k,5}{2,1});
            J_PR_ang(:,k) = real(K{k,6}{1,1});
            J_QR_ang(:,k) = real(K{k,6}{2,1});

        elseif( sum( K{k,1} == idx3_augmented.vscdc_pq ))
            J_PR_abs(:,k) = real(K{k,5}{1,1});
            J_PR_ang(:,k) = real(K{k,6}{1,1});

        elseif( sum( K{k,1} == idx3_augmented.vscdc_vq ))
            J_VR_abs(:,k) = real(K{k,5}{1,1});
            J_VR_ang(:,k) = real(K{k,6}{1,1});

        else
            warning('somethings off mate')
        end
    end

    Kp_abs = J_PR_abs;%(4:end,4:end);
    Kq_abs = J_QR_abs;%(4:end,4:end);
    Kv_abs = J_VR_abs;%(4:end,4:end);
    Kp_ang = J_PR_abs;%(4:end,4:end);
    Kq_ang = J_QR_abs;%(4:end,4:end);
    Kv_ang = J_VR_abs;%(4:end,4:end);
end