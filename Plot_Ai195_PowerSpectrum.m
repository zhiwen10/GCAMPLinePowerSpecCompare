%%
T = readtable('mouseID.xlsx');
T = T(1:end-1,:);
mouseN = size(T,1);
%%
C = categorical(T.Line);
ncolor = grp2idx(C);
uniqueC = unique(ncolor);
uniqueCN = length(uniqueC);
colors = cbrewer('qual','Set1',uniqueCN);
count1 = 1;
text1 = {'V1 point1','V1_point2','PPC_point1','PPC_point2','PFC_point1','PFC_point2'};
%%
for kk = 1:uniqueCN
    indx = find(ncolor == uniqueC(kk));  
    psdxTemp = psdx(:,:,indx,:);  
    psdxAll{kk} = (psdxTemp);
    lineName{kk} = T.Line{indx(1)};
end
%%
% lineName1 = lineName(1,[2,3,7]);
% psdxAll1 = psdxAll(1,[2,3,7]);
lineName1 = lineName(1,[6,7,8]);
psdxAll1 = psdxAll(1,[6,7,8]);
lineN = 3;
colors = cbrewer('qual','Set1',lineN);
text1 = {'V1 point1','V1_point2','PPC_point1','PPC_point2','PFC_point1','PFC_point2'};
figure('position',[500,200,900,800]);   
count1 = 1;
for m = 1:6
    for n = 1:2   
        subplot(6,2,count1)
        count2 = 1;
        for kk = 1:lineN
            psdxM1 = psdxAll1{kk}; 
            for qq = 1:size(psdxM1,3)
                p{count2} = plot(freq(:,1,1),psdxM1(:,m,qq,n),'color',colors(kk,:),'LineWidth',2);
                xlabel('Frequency (Hz)');
                ylabel({'(dF/F)^2'})   
                % ylabel('dB/Hz')
                ylim([0 6*10.^(-3)]);               
                % ylim([0 2*10.^(-4)]);
                % xlim([0,12])
                hold on
                count2 = count2+1;
            end
        end
        % text(15,3*10.^(-3),text1{m}, 'Interpreter', 'none');
        text(15,1*10.^(-4),text1{m}, 'Interpreter', 'none');
        if count1 == 1
            %legend([p{1},p{3},p{4}],{'Ai195.Cre','Ai195.flpO.','Triple'})
            legend([p{1},p{3},p{4}],{'g7f','g7s','Triple'})
            title('RF Mapping')
        elseif count1 == 2
            title('Linearity test')
        end
        count1 = count1+1;
    end
end
title1 = 'Ai195_PowerSpectrum_Normalisation';
sgtitle(title1, 'Interpreter', 'none')
savefig(['Figures/' title1])
saveas(gcf, ['Figures/' title1 '.png'])