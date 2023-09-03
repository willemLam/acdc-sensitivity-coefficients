function [Kp_abs, Kq_abs, Kv_abs,Kp_ang, Kq_ang, Kv_ang] = transform_K_polar(K,Grid_para, idx)

    J_PR_abs = zeros(Grid_para.n_nodes);
    J_QR_abs = zeros(Grid_para.n_nodes);
    J_VR_abs = zeros(Grid_para.n_nodes);
    J_PR_ang = zeros(Grid_para.n_nodes);
    J_QR_ang = zeros(Grid_para.n_nodes);
    J_VR_ang = zeros(Grid_para.n_nodes);
    for k = 1:size(K,1)
        if( sum( K{k,1} == idx.slack))
            continue
        elseif( sum( K{k,1} == idx.pqac))
            J_PR_abs(:,k) = (K{k,5}{1,1});
            J_QR_abs(:,k) = (K{k,5}{2,1});
            J_PR_ang(:,k) = (K{k,6}{1,1});
            J_QR_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.pvac ))
            J_PR_abs(:,k) = (K{k,5}{1,1});
            J_VR_abs(:,k) = (K{k,5}{2,1});
            J_PR_ang(:,k) = (K{k,6}{1,1});
            J_VR_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.pdc ) )
            J_PR_abs(:,k) = (K{k,5}{1,1});
            J_PR_ang(:,k) = (K{k,6}{1,1});

        elseif( sum( K{k,1} == idx.vdc ) )
            J_VR_abs(:,k) = (K{k,5}{1,1});
            J_VR_ang(:,k) = (K{k,6}{1,1});

        elseif( sum( K{k,1} == idx.vscac_pq))
            J_PR_abs(:,k) = (K{k,5}{1,1});    
            J_QR_abs(:,k) = (K{k,5}{2,1});
            J_PR_ang(:,k) = (K{k,6}{1,1});    
            J_QR_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.vscac_vq))
            J_PR_abs(:,k) = (K{k,5}{1,1});
            J_QR_abs(:,k) = (K{k,5}{2,1});
            J_PR_ang(:,k) = (K{k,6}{1,1});
            J_QR_ang(:,k) = (K{k,6}{2,1});

        elseif( sum( K{k,1} == idx.vscdc_pq ))
            J_PR_abs(:,k) = (K{k,5}{1,1});
            J_PR_ang(:,k) = (K{k,6}{1,1});

        elseif( sum( K{k,1} == idx.vscdc_vq ))
            J_VR_abs(:,k) = (K{k,5}{1,1});
            J_VR_ang(:,k) = (K{k,6}{1,1});

        else
            warning('somethings off mate')
        end
    end

    Kp_abs = J_PR_abs;%(4:end,4:end);
    Kq_abs = J_QR_abs;%(4:end,4:end);
    Kv_abs = J_VR_abs;%(4:end,4:end);
    Kp_ang = J_PR_ang;%(4:end,4:end);
    Kq_ang = J_QR_ang;%(4:end,4:end);
    Kv_ang = J_VR_ang;%(4:end,4:end);
end