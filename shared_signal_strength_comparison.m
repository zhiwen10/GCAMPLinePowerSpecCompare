startup_linearity;

% This script loads all of the data and splits it by the line used to
% express GCaMP. It then plots (1) comparison plots of the different
% conditions against line, and (2) computes SNR

%% Load the probe data and the widefield data
preprocDir = fullfile(githubDir,'LinearityTest','preproc');

files = dir(fullfile(preprocDir,'*widefield.mat'));
%%
gcampLines = struct; 

gcampLines.AD_0011 = 0;
gcampLines.DPR_0005 = 0;
gcampLines.ZYE_0003 = 0;
line_names = {'6s - ai162'};
gcampLines.DPR_0008 = 1;
gcampLines.ZYE_0009 = 1;
gcampLines.ZYE_0010 = 1;
gcampLines.DPR_0012 = 1;
line_names{end+1} = '7s - ai195 triple';
gcampLines.DPR_0009 = 2;
gcampLines.ZYE_0011 = 2;
line_names{end+1} = '7f - ai210 double cre/flpo inj';
gcampLines.ZYE_0016 = 3;
gcampLines.ZYE_0012 = 3;
gcampLines.DPR_0011 = 3;
line_names{end+1} = '7f - ai210 flpo inj';
gcampLines.ZYE_0015 = 4;
line_names{end+1} = '6s - ai162 cre inj';
gcampLines.ZYE_0019 = 5;
gcampLines.ZYE_0022 = 5;
line_names{end+1} = '7f - flex inj';
gcampLines.ZYE_0014 = 6;
gcampLines.ZYE_0025 = 6;
line_names{end+1} = '7s - ai195 single flpo or cre inj';
gcampLines.none = 7;
line_names{end+1} = '7s - flex inj';
gcampLines.DO_0004 = 8;
gcampLines.ZYE_0018 = 8;
gcampLines.DO_0007 = 8;
line_names{end+1} = 'soma - ignore';
%%
ignore_lines = [2 8];

gcamp_lineID = fields(gcampLines);
gcamp_line = cellfun(@(x) gcampLines.(x),gcamp_lineID);
gcamp_lineID{end+1} = '00';
gcamp_line(end+1) = 0;

widefield = {};
for fi = 1:length(files)
    fname = fullfile(files(fi).folder,files(fi).name);
    temp = load(fname);
    temp.sessionData.typenum = gcampLines.(temp.sessionData.mn);

    widefield{end+1} = temp.sessionData;
    widefield{end}.fname = fname;
end
disp(sprintf('Found %i GCaMP sessions',length(widefield)));

%% Drop lines that are failed for any reason
% lines whose data we will drop
% line 2 - failed because double injections don't work well 
keep = zeros(size(widefield));
for wi = 1:length(widefield)
    keep(wi) = all(widefield{wi}.typenum~=ignore_lines);
end
widefield = widefield(logical(keep));
disp(sprintf('Keeping %i GCaMP sessions',length(widefield)));

%% Get types
widefield_line = [];
for wi = 1:length(widefield)
    idx = cellfun(@(wide) find(cellfun(@(line) ~isempty(strfind(wide.fname,line)),gcamp_lineID),1),widefield);
    widefield_line = gcamp_line(idx);
end

%% Remove frame counts that are annoying
for wi = 1:length(widefield)
    widefield{wi} = lt_mergeFrameCounts(widefield{wi});
end

%% Get all the adata and figure out what contrasts and duratinos are in the dataset
alldata = [];
for wi = 1:length(widefield)
    alldata = [alldata ; widefield{wi}.adata];
end

contrasts = unique(alldata(:,2));
durations = unique(alldata(:,7)); % in frames
durations(durations==0) = [];

%% Sort responses for probe and widefield according to contrast/duration
for wi = 1:length(widefield)
    widefield{wi} = lt_sortVars(widefield{wi},contrasts,durations);
end

% resample
time = widefield{end}.time;
for wi = 1:length(widefield)
    
    % resample
    widefield{wi} = lt_resampleToProbe(widefield{wi},time,contrasts,durations);
    widefield{wi}.time = time;
end

%% Separate by line
uline = unique(widefield_line);
lineNames = line_names(uline+1);
lines = {};
for ui = 1:length(uline)
    cline = uline(ui);
    
    idxs = find(widefield_line==cline);
    clines = widefield(idxs);
    
    lines{end+1} = clines;
end

%% Plot each dataset separately (mean over repeats) to look for issues
close all
for ui = 1:length(lines)
    clines = lines{ui};
    for wi = 1:length(clines)
        lt_quickPlot(clines{wi},contrasts,durations,line_names{uline(ui)+1});
    end
end 
%% Plot each dataset separately to look for issues
close all
for ui = 1:length(lines)
    clines = lines{ui};
    for wi = 1:length(clines)
        lt_quickPlot_nomean(clines{wi},contrasts,durations,line_names{uline(ui)+1});
    end
end 


%% Pull the baseline variance
for ui = 1:length(lines)
    clines = lines{ui};
    
    cBaseStd = [];
    for li = 1:length(clines)
        data = clines{li};
        
        baseIdxs = data.time<0;
        
        sortedResp = data.sortedResp;
        cBaseline = [];
        figure;
        for ci = 1:size(sortedResp,1)
            for di = 1:size(sortedResp,2)
                resp = sortedResp{ci,di};
                temp = resp(:,baseIdxs);
                temp = temp - repmat(mean(temp,2),1,sum(baseIdxs));
                temp = temp';
                cBaseline = [cBaseline ; temp(:)];
                subplot(6,7,(ci-1)*7+di);
                plot(temp(:));
                a = axis;
                axis([a(1) a(2) -0.05 0.05]);
            end
        end
        title(lineNames{ui});
        
        cBaseStd(li) = nanstd(cBaseline);
    end
    
    lineStd{ui} = cBaseStd;
end

% figure
h = figure(4); clf; hold on
for li = 1:length(lines)
    plot(li,lineStd{li},'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','w','MarkerSize',5);
    plot(li,mean(lineStd{li}),'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','w','MarkerSize',8);
end

set(gca,'XTick',1:length(lines),'XTickLabel',lineNames);
ylabel('Baseline period signal variability (s.d. of % signal change)');

drawPublishAxis('figSize=[22,8]');

savepdf(h,fullfile('C:\proj\LinearityTest\figures\baseline.pdf'));

%% Plot the different lines on top of each other
% first, average by line
line_avgs = {};
line_avgsz = {};
for ui = 1:length(lines)
    clines = lines{ui};
    cWideData = lt_pullSorted(clines,contrasts,durations);
    % z-score
    cWidez = lt_zscore(cWideData,time);
    % squash repeats
    cWideMeanz = lt_squashSorted(cWidez,false,length(time));
    line_avgsz{end+1} = squeeze(nanmean(cWideMeanz,1));
    % non-zscored
    
%     cWidez = lt_dff(cWideData,time);
    cWideMean = lt_squashSorted(cWideData,false,length(time));
    line_avgs{end+1} = squeeze(nanmean(cWideMean,1));
end


%% Plot all widefield on top of each other (to look at snr)
cmap = brewermap(8,'Dark2');
cmap = [0 0 0; cmap]; % put black first, so that 6s is a black line
h = figure(1); clf; hold on
for ci = 1:length(contrasts)
    for di = 1:length(durations)
        for ui = 1:length(line_avgs)
            if strfind(lineNames{ui},'7f')
                lt = '--';
            else
                lt = '-';
            end
            p(ui) = plot(di*2.5+time,squeeze(line_avgsz{ui}(ci,di,:)) + ci*6,lt,'Color',cmap(ui,:),'LineWidth',2);
        end        
    end
end
plot([1 2],[5 5],'-k');
text(1.1,4,'1 s');
plot([1 1],[5 6],'-k');
t = text(0.7,4,'1 s.d.');
set(t,'Rotation',90);
axis([0 20 0 48]);

set(gca,'XTick',[],'YTick',[]);

legend(p,line_names(uline+1),'Location','NorthWest');

savepdf(h,fullfile('C:\proj\LinearityTest\figures','snr_alllinesz.pdf'));

%% Do it again, but just the top 4 responses

h = figure(1); clf; hold on
for ci = 5:length(contrasts)
    for di = 6:length(durations)
        for ui = 1:length(line_avgs)
            if strfind(lineNames{ui},'7f')
                lt = '--';
            else
                lt = '-';
            end
            p(ui) = plot((di-5)*3+time,squeeze(line_avgsz{ui}(ci,di,:)) + (ci-4)*10,lt,'Color',cmap(ui,:),'LineWidth',3);
        end        
    end
end
plot([1 2],[5 5],'-k');
text(1.1,4,'1 s');
plot([1 1],[5 6],'-k');
t = text(0.7,4,'1 s.d.');
set(t,'Rotation',90);
axis([0 10 0 30]);

set(gca,'XTick',[3 6],'YTick',[10 20]);
set(gca,'XTickLabel',{'0.5 s','1 s'},'YTickLabel',{'50%','100%'});

% get the sort order
mval = [];
for ui = 1:length(line_avgs)
    mval(ui) = mean(line_avgsz{ui}(end,end,time>0.1 & time<0.4));
end
[~,sortOrder] = sort(mval,'descend');

% sort the legend order by the response strength
sortedP = p(sortOrder);
names = line_names(uline+1);
sortedN = names(sortOrder);
legend(sortedP,sortedN);

savepdf(h,fullfile('C:\proj\LinearityTest\figures','snr_alllines_smlz.pdf'));

%% Do it again, but just the top 4 responses

h = figure(2); clf; hold on
for ci = 5:length(contrasts)
    for di = 6:length(durations)
        for ui = 1:length(line_avgs)
            if strfind(lineNames{ui},'7f')
                lt = '--';
            else
                lt = '-';
            end
            p(ui) = plot((di-5)*3+time,squeeze(line_avgs{ui}(ci,di,:)) + (ci-4)*0.15,lt,'Color',cmap(ui,:),'LineWidth',3);
        end        
    end
end
% get the sort order
mval = [];
for ui = 1:length(line_avgs)
    mval(ui) = mean(line_avgs{ui}(end,end,time>0.1 & time<0.4));
end
plot([2 3],[0 0],'-k');
text(2.2,-0.025,'1 s');
plot([2 2],[0 0.1],'-k');
t = text(1.8,0.05,'0.1 dF/F');
set(t,'Rotation',90);
axis([0 10 -0.1 0.45]);
% 
set(gca,'XTick',[3 6],'YTick',[0.15 0.3]);
set(gca,'XTickLabel',{'1 s','2 s'},'YTickLabel',{'50%','100%'});

legend(p,line_names(uline+1),'Location','SouthEast');

savepdf(h,fullfile('C:\proj\LinearityTest\figures','snr_alllines_sml.pdf'));