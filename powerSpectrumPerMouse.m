%%
C = categorical(T.Line);
ncolor = grp2idx(C);
uniqueC = length(unique(ncolor));
colors = cbrewer('qual','Set1',uniqueC);
text1 = {'V1 point1','V1_point2','PPC_point1','PPC_point2','PFC_point1','PFC_point2'};
for kk = 1:mouseN
    figure('position',[500,500,900,600]);
    count1 = 1;
    for m = 1:6
        for n = 1:2
            subplot(6,2,count1)
            plot(squeeze(freq(:,kk,n)),squeeze(psdx(:,m,kk,n)),'color',colors(ncolor(kk),:),'LineWidth',2)
            hold on
            xlabel('Frequency (Hz)');
            ylabel({'Normalized Power', '(dF/F)^2'})   
            ylim([0 8*10.^(-3)]);
            text(15,4*10.^(-3),text1{m}, 'Interpreter', 'none');
            if count1 == 1
                % legend([p{1},p{5},p{7},p{10},p{13},p{14}],{'Triple','Ai210','Ai195','Ai162','g7f','Ai203'})
                title('RF Mapping')
            elseif count1 == 2
                title('Linearity test');
            end
            count1 = count1+1;
        end
    end
    title1 = sprintf(T.MouseID{kk});
    title2 = sprintf(T.Line{kk});
    sgtitle({title1,title2}, 'Interpreter', 'none')
end
