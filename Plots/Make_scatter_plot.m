function f = Make_scatter_plot(r,c,n_nodes,mode)
    L3 = [];

    L1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
    for i = 1:length(L1)
        L3 =  [ L3; convertCharsToStrings(['$' num2str(L1(i)) '_a$']);
                        convertCharsToStrings(['$' num2str(L1(i)) '_b$']);
                        convertCharsToStrings(['$' num2str(L1(i)) '_c$'])];
    end

    LDC = [];
    for i = 19:26
        LDC =  [ LDC; convertCharsToStrings(['$' num2str(i)  '$'])];
    end

    n_L = [L3;LDC];

    n_name = {'$|$','$V_{AC}$','$|$','$V_{DC}$','$|$'}; 
    n_index = [0.5, 24, 54.5, 58, 62];
    n_line = [0.75 54.25;54.75 61.75];

    f = figure('Renderer', 'painters', 'Position', [10 10 1500 500]);
    subplot(2,1,1)
    title(strcat('SC dEi/dX for X = ',mode)) % -  max abs error = ',string(max(abs(r-c))))
    hold on
    scatter(1:n_nodes,real(r),40,'fill')
    scatter(1:n_nodes,real(c),40,'fill')
    ylabel('real part')
    legend({'analytical sens coef','ground truth'}) %'analytic senc coef','inverse jacob',
    set(gca,'FontSize',15) 
    xticks(1:62)
    xticklabels(n_L)
    set(gca,'FontSize',15) 
     pos = get(gca, 'Position');
        pos(1) = 0.055;
        pos(3) = 0.9;
        set(gca, 'Position', pos)


    ylim([min(ylim(gca)),max(ylim(gca))])
    xlim([0.5,62.5])

    hold off

    subplot(2,1,2)
    hold on
    scatter(1:n_nodes,imag(r),40,'fill')
    scatter(1:n_nodes,imag(c),40,'fill')
    ylabel('imaginary part')
    xticks(1:62)
    xticklabels(n_L)

    set(gca,'FontSize',15) 
     pos = get(gca, 'Position');
        pos(1) = 0.055;
        pos(3) = 0.9;
        set(gca, 'Position', pos)

    lowerLevel = min(ylim(gca))-range(ylim(gca))*.2; 
    lowerLabels = n_name; 
    text(n_index, repmat(lowerLevel,1,numel(lowerLabels)), lowerLabels,...
        'VerticalAlignment','Top','HorizontalAlignment','Center','Interpreter','Latex', 'FontWeight','bold')
    ylim([min(ylim(gca)),max(ylim(gca))])
    xlim([0.5,62.5])

    for i = 1:length(n_line)
        [figx, figy] = axxy2figxy(gca, n_line(i,:), repmat(lowerLevel,1,numel(n_line(i,:))));
        annotation('line',figx,figy)
    end

    hold off

end