function I_flow = get_Current_flow(E,Grid_para)

    YYL = Grid_para.YYL;
    YYT = Grid_para.YYT;
    I_flow = zeros(size( Grid_para.Ampacities));

    for i = 1:size(YYL,1)
        for j = 1:size(YYL,2)
                I_flow(ceil(i),ceil(j)) = YYL(i,j)*(E(i) - E(j)) + YYT(i,j)*E(i);
        end
    end 

end