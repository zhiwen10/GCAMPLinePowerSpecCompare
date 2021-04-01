%%
T = readtable('mouseID.xlsx');
T = T(1:end-1,:);
mouseN = size(T,1);
%%
% T.Line(strcmp('Ai210',string(T.Line))| strcmp('Ai195',string(T.Line))) = {'PHP inj'};
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
    sizeIndx = length(indx);
    colors = cbrewer('qual','Set1',sizeIndx);
    figure('position',[500,500,900,800]);
    count1 = 1;
    for m = 1:6
        for n = 1:2
            subplot(6,2,count1)
            for j = 1:sizeIndx                
                plot(squeeze(freq(:,indx(j),n)),squeeze(psdx(:,m,indx(j),n)),'color',colors(j,:),'LineWidth',2)
                hold on
                xlabel('Frequency (Hz)');
                ylabel({'Normalized Power', '(dF/F)^2'})   
                ylim([0 8*10.^(-3)]);
                text(15,4*10.^(-3),text1{m}, 'Interpreter', 'none');
            end
            if count1 == 1
                % legend([p{1},p{5},p{7},p{10},p{13},p{14}],{'Triple','Ai210','Ai195','Ai162','g7f','Ai203'})
                title('RF Mapping')
            elseif count1 == 2
                title('Linearity test');
            end
            count1 = count1+1;
        end
    end
    title1 = sprintf(T.Line{indx(1)});
    sgtitle(title1, 'Interpreter', 'none')
    savefig(['Figures/SecondHalf_' title1])
end
