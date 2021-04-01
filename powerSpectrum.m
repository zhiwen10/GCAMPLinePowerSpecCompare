% add path for dependencies 
githubDir = 'C:\Users\Steinmetz lab\Documents\git';

addpath(genpath(fullfile(githubDir, 'widefield'))) % cortex-lab/widefield
addpath(genpath(fullfile(githubDir, 'Pipelines')))
addpath(genpath(fullfile(githubDir, 'wheelAnalysis'))) % cortex-lab/wheelAnalysis
%%
T = readtable('mouseID.xlsx');
T = T(1:end-1,:);
mouseN = size(T,1);
%%
for n = 1:2        
    % set pixel data
    points(:,:,1) = [T.v1X,T.v1Y];
    points(:,:,2) = [T.v2X,T.v2Y];
    points(:,:,3) = repmat([330,330],[mouseN,1]);
    points(:,:,4) = repmat([200,300],[mouseN,1]);
    points(:,:,5) = repmat([120,330],[mouseN,1]);
    points(:,:,6) = repmat([80,330],[mouseN,1]);
    % set session 
    formatOut = 'yyyy-mm-dd';
    if n == 1
    tdAll = datestr(T.rfDate,formatOut); enAll = T.rfMain; enexpAll = T.rfFolder; nSVAll = T.nSV2;
    else 
    tdAll = datestr(T.ExpDate,formatOut); enAll = T.mainFolder; enexpAll = T.expFolder; nSVAll = T.nSV1;
    end
    % loop through mouse ID for powerRatio and PowerSpectrum(psdx)
    for i = 1:mouseN
    mn = char(T.MouseID(i));
    td = tdAll(i,:);
    en = enAll(i);
    enexp = enexpAll(i);
    if not(isnan(en))
        serverRoot = expPath(mn, td, en);
        serverRoot2 = expPath(mn, td, enexp);
        nSV = nSVAll(i);
        point = squeeze(points(i,:,:))';
        [pRatio(:,i,n),freq(:,i,n),psdx(:,:,i,n)] = powerRatio(nSV,serverRoot,serverRoot2,point);
    else
        pRatio(:,i,n) = nan;
    end
    end
end
%%
pRatioRF = squeeze(pRatio(:,:,1))';
pRatioLinear = squeeze(pRatio(:,:,2))';
T.pRatioRF = pRatioRF;
T.pRatioLinear = pRatioLinear;
%%
C = categorical(T.Line);
ncolor = grp2idx(C);
uniqueC = length(unique(ncolor));
colors = cbrewer('qual','Set1',uniqueC);
figure('position',[500,100,900,800]);
count1 = 1;
text1 = {'V1 point1','V1_point2','PPC_point1','PPC_point2','PFC_point1','PFC_point2'};
for m = 1:6
    for n = 1:2
        subplot(6,2,count1)
        for kk = 1:mouseN
            p{kk} = plot(squeeze(freq(:,kk,n)),squeeze(psdx(:,m,kk,n)),'color',colors(ncolor(kk),:),'LineWidth',2)
            hold on
            xlabel('Frequency (Hz)');
            ylabel({'Normalized Power', '(dF/F)^2'})   
            %ylim([0 8*10.^(-3)]);
            % ylim([0 2*10.^(-4)]);
        end
        % text(15,4*10.^(-3),text1{m}, 'Interpreter', 'none');
        if count1 == 1
            legend([p{1},p{5},p{7},p{10},p{13},p{14}],{'Triple','Ai210','Ai195','Ai162','g7f','Ai203'})
            title('RF Mapping')
        elseif count1 == 2
            title('Linearity test')
        end
        count1 = count1+1;
    end
end

%%
function [pRatio,freqQuery,psdxInterp] = powerRatio(nSV,serverRoot,serverRoot2,points)
    %plot PCA components and traces for blue and purple channels
    
    corrPath = fullfile(serverRoot, 'corr', 'svdTemporalComponents_corr.npy');
    if ~exist(corrPath, 'file')
        colors = {'blue', 'violet'};
        computeWidefieldTimestamps(serverRoot, colors); % preprocess video
        [U, V, t, mimg] = hemoCorrect(serverRoot, nSV); % process hemodynamic correction
    else
        [U, V, t, mimg] = loadUVt(serverRoot, nSV);
    end
    if length(t) > size(V,2)
      %t = t(1:end-1);
      t = t(1:size(V,2));
    end
    %%
    pixelCorrelationViewerSVD(U,V)
    dV = [zeros(size(V,1),1) diff(V,[],2)];
    ddV = [zeros(size(dV,1),1) diff(dV,[],2)];
    close all
    %%
    tFile = fullfile(serverRoot2, ['photodiodeFlips.timestamps.npy']); 
    ts = readNPY(tFile); 
    ts(isnan(ts)) = [];
    %%
    for i = 1:size(points,1)
        point_v1 = points(i,:);
        px_v1 = squeeze(U(point_v1(1),point_v1(2),:))'*dV;
        %%
        tlFile = fullfile(serverRoot, 'blue\meanImage.npy');
        meanImage = readNPY(tlFile);
        meanF = meanImage(point_v1(1),point_v1(2));
        %%
        df_v1 = px_v1/meanF;
        % qx = ts(1,1):1/35:ts(end,1);
        qx = (ts(1,1)+(ts(end,1)-ts(1,1))/2):1/35:ts(end,1);
        df2_v1 = interp1(t,df_v1,qx);
        fr = 35;
        [freq,psdx] = powerSpectrumPlot(df2_v1,fr);
        % plot(freq,10*log10(psdx))
        frange = [2,4];
        pBand = bandpower(df2_v1,fr,frange);
        pTot = bandpower(df2_v1,fr,[0.5,17]);
        pRatio(i,1) = pBand./pTot;    
        %%
        gw = myGaussWin(0.5, 35); 
        freqTemp = freq; psdxTemp = psdx;
        freqQuery = 0.5:0.01:17.4;
        psdxTemp1 = interp1(freqTemp,psdxTemp,freqQuery);
        psdxSum = sum(psdxTemp1);
        psdxTemp1 = psdxTemp1/psdxSum;
        psdxInterp(:,i) = conv(psdxTemp1,gw,'same');
    end
end

function [freq,psdx] = powerSpectrumPlot(trace,fr)
    N = length(trace);
    xdft = fft(trace);
    xdft = xdft(1:N/2+1);
    psdx = (1/(fr*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fr/length(trace):fr/2;
end