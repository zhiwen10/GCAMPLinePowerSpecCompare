%%
T = readtable('mouseID.xlsx');
T = T(1:end-1,:);
mouseN = size(T,1);
%%
T.Line(strcmp('Ai210',string(T.Line))| strcmp('Ai195',string(T.Line))) = {'PHP inj'};
%%
% Check = {'Ai162','triple','PHP inj'};
% indx1 = cellfun(@(x) strcmp(x,string(T.Line)),Check,'UniformOutput',false);
% indx2 = cat(2,indx1{:});
% T = T(any(indx2,2),:);
%%
C = categorical(T.Line);
ncolor = grp2idx(C);
uniqueC = unique(ncolor);
uniqueCN = length(uniqueC);
text1 = {'V1 point1','V1_point2','PPC_point1','PPC_point2','PFC_point1','PFC_point2'};
%%
for kk = 1:uniqueCN
    indx = find(ncolor == uniqueC(kk));  
    psdxTemp = squeeze(psdx(:,:,indx,:));  
    psdxMean{kk} = squeeze(mean(psdxTemp,3));
    psdxStd{kk} = squeeze(std(psdxTemp,0,3));
    psdxStdev{kk} = squeeze(std(psdxTemp,0,3))/sqrt(size(psdxTemp,3));
    lineName{kk} = T.Line{indx(1)};
end
%%
lineName1 = lineName(1,[1,3,5]);
psdxMean1 = psdxMean(1,[1,3,5]);
psdxStd1 = psdxStd(1,[1,3,5]);
psdxStdev1 = psdxStdev(1,[1,3,5]);
lineN = 3;
colors = cbrewer('qual','Set1',lineN);
text1 = {'V1 point1','V1_point2','PPC_point1','PPC_point2','PFC_point1','PFC_point2'};
figure('position',[500,500,900,800]);   
count1 = 1;
for m = 1:6
    for n = 1:2   
        subplot(6,2,count1)
        for kk = 1:lineN
            psdxM1 = psdxMean1{kk};
            psdxS1 = psdxStdev1{kk};   
            shadedErrorBar(freq(:,1,1),psdxM1(:,m,n),psdxS1(:,m,n),'lineprops', {'color',colors(kk,:)},'transparent',true);
            xlabel('Frequency (Hz)');
            ylabel({'(dF/F)^2'})   
            ylim([0 6*10.^(-3)]);
            hold on
        end
        text(15,3*10.^(-3),text1{m}, 'Interpreter', 'none');
        if count1 == 1
            legend({'Ai162','PHP inj.','Triple'})
            title('RF Mapping')
        elseif count1 == 2
            title('Linearity test')
        end
        count1 = count1+1;
    end
end
title1 = 'meanPowerSpectrum';
sgtitle(title1, 'Interpreter', 'none')
savefig(['Figures/SecondHalf_' title1])

            