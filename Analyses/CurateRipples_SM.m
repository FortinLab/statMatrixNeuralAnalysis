%% CurateRipples(rips, trialRips)
global listSel
listSel = 1;        % Used to keep track of which list is being selected from for ripple viewing
%%
rips = RippleDetection_SM;
[trialRips] = ExtractTrialEventRips_SM(rips, [500 700]);
allTrialRips = sortrows(cell2mat([trialRips.Events(:,1); trialRips.Events(:,2); trialRips.Events(:,3)]));
%% Create Figure
ripCure = figure;
rawAxes = axes(ripCure, 'position', [0.1, 0.75, 0.7, 0.2]);
bpfAxes = axes(ripCure, 'position', [0.1, 0.55, 0.7, 0.2]);
spkAxes = axes(ripCure, 'position', [0.1, 0.2, 0.7, 0.3]);
linkaxes([spkAxes, bpfAxes, rawAxes], 'x');

% Overall Ripple List
ssnRipTitleAx = axes(ripCure, 'position', [0.825, 0.925, 0.11, 0.05]);
set(ssnRipTitleAx, 'xlim', [-0.5 0.5], 'ylim', [0 0.5]);
text(ssnRipTitleAx, 0,0, 'Session Rips', 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
axis(ssnRipTitleAx, 'off');
ssnRipList = uicontrol(ripCure, 'Units', 'Normalized', 'Style', 'listbox', 'String', 1:size(rips.Ripples.Events,1),...
    'Position', [0.825, 0.55, 0.15, 0.375], 'Tag', 'ssnRip_Lst', 'Callback', @SelectSsnRip, 'userData', rips.TimeStamps(rips.Ripples.Events));

% Trial Ripple List
trlRipTitleAx = axes(ripCure, 'position', [0.825, 0.5, 0.11, 0.05]);
set(trlRipTitleAx, 'xlim', [-0.5 0.5], 'ylim', [0 0.5]);
text(trlRipTitleAx, 0,0, 'Trial Rips', 'horizontalalignment', 'center', 'verticalalignment', 'bottom');
axis(trlRipTitleAx, 'off');
trlRipList = uicontrol(ripCure, 'Units', 'Normalized', 'Style', 'listbox', 'String', 1:size(allTrialRips,1),...
    'Position', [0.825, 0.2, 0.15, 0.300], 'Tag', 'trlRip_Lst', 'Callback', @SelectTrlRip, 'userData', rips.TimeStamps(allTrialRips));

% Zoom In
zoomOutBtn = uicontrol(ripCure, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Zm-',...
    'Position', [0.35, 0.06, 0.05, 0.075], 'Callback', @ZoomOut);
% Zoom Out
zoomInBtn = uicontrol(ripCure, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Zm+',...
    'Position', [0.65, 0.06, 0.05, 0.075], 'Callback', @ZoomIn);

% Previous Ripple
prevRipBtn = uicontrol(ripCure, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', '<< Previous Ripple',...
    'Position', [0.05, 0.05, 0.2, 0.1], 'Callback', @PrevRip);
% Next Ripple
nextRipBtn = uicontrol(ripCure, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Next Ripple >>',...
    'Position', [0.75, 0.05, 0.2, 0.1], 'Callback', @NextRip);

% Remove Ripple
removeRipBtn = uicontrol(ripCure, 'Units', 'Normalized', 'Style', 'pushbutton', 'String', 'Remove Ripple',...
    'Position', [0.425, 0.075, 0.2, 0.05], 'Callback', @RmvRip);

%% Plot Stuff
% Plot Raw LFP Traces
raw = plot(rawAxes, rips.TimeStamps, rips.SessionData.RawLFP, 'color', 'k');
for r = 1:length(raw)
    raw(r).Color(4) = 0.2;
end
set(rawAxes, 'Tag', 'Raw_Axes');
% Plot Ripple Band LFP Traces
bpf = plot(bpfAxes, rips.TimeStamps, rips.SessionData.RipBPF, 'color', 'k');
for b = 1:length(bpf)
    bpf(b).Color(4) = 0.2;
end
set(bpfAxes, 'Tag', 'Bpf_Axes');
% Plot Spiking Activity
[spkX, spkY] = find(rips.SessionData.Spikes~=0);
scatter(spkAxes, rips.TimeStamps(spkX),spkY, '*k');
% line(spkAxes, rips.TimeStamps(spkX),spkY, 'marker', '*', 'linestyle', 'none', 'color', 'k');
xlabel(spkAxes, 'Time (m)');
set(spkAxes, 'Tag', 'Spk_Axes');

% Fiddle with Axes
set(rawAxes, 'color', 'none', 'ycolor', 'none', 'xcolor', 'none', 'xticklabel', [], 'box', 'off');
set(bpfAxes,'color', 'none', 'ycolor', 'none', 'xticklabel', [], 'box', 'off');
set(spkAxes, 'color', 'none', 'ycolor', 'none', 'box', 'off');
axis([rawAxes, bpfAxes,spkAxes], 'tight');
rawLim = repmat(get(rawAxes,'ylim'),[2,1]);
bpfLim = repmat(get(bpfAxes,'ylim'),[2,1]);
spkLim = repmat(get(spkAxes,'ylim'),[2,1]);
spkLim(:,1) = spkLim(:,1)-1;
spkLim(:,2) = spkLim(:,2)+1;
% Plot Trial Periods
for trl = 1:size(rips.TrialInfo.TrialPokes,1)
    switch rips.TrialInfo.OdorVect(trl)
        case 1
            patchColor = [44/255 168/255 224/255];
        case 2
            patchColor = [154/255 133/255 122/255];
        case 3
            patchColor = [9/255 161/255 74/255];
        case 4 
            patchColor = [128/255 66/255 151/255];
        case 5
            patchColor = [241/255 103/255 36/255];
    end
    patch(rawAxes, 'XData', [rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,1)), rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,2)), rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,2)),rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,1))],...
        'YData', rawLim(:),...
        'FaceColor', patchColor, 'FaceAlpha', 0.5,...
        'EdgeColor', patchColor, 'EdgeAlpha', 0.5);
    patch(bpfAxes, 'XData', [rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,1)), rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,2)), rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,2)),rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,1))],...
        'YData', bpfLim(:),...
        'FaceColor', patchColor, 'FaceAlpha', 0.5,...
        'EdgeColor', patchColor, 'EdgeAlpha', 0.5);
    patch(spkAxes, 'XData', [rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,1)), rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,2)), rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,2)),rips.TimeStamps(rips.TrialInfo.TrialPokes(trl,1))],...
        'YData', spkLim(:),...
        'FaceColor', patchColor, 'FaceAlpha', 0.5,...
        'EdgeColor', patchColor, 'EdgeAlpha', 0.5);    
end

for rip = 1:size(rips.Ripples.Events,1)
    curRipX = [rips.TimeStamps(rips.Ripples.Events(rip,:)); flipud(rips.TimeStamps(rips.Ripples.Events(rip,:)))]';
    patch(rawAxes, 'XData', curRipX,...
        'YData', rawLim(:),...
        'FaceColor', 'y', 'FaceAlpha', 0.15,...
        'EdgeColor', 'y', 'EdgeAlpha', 0.5);
    patch(bpfAxes, 'XData', curRipX,...
        'YData', bpfLim(:),...
        'FaceColor', 'y', 'FaceAlpha', 0.15,...
        'EdgeColor', 'y', 'EdgeAlpha', 0.5);
    patch(spkAxes, 'XData', curRipX,...
        'YData', spkLim(:),...
        'FaceColor', 'y', 'FaceAlpha', 0.15,...
        'EdgeColor', 'y', 'EdgeAlpha', 0.5);    
end

%% Callbacks
function SelectSsnRip(source,event)
global listSel
listSel = 1;
kids = get(get(source, 'Parent'), 'Children');
tags = arrayfun(@(a)a.Tag, kids, 'uniformoutput',0);
axTag = strcmp(tags, 'Raw_Axes');
set(kids(axTag), 'xLim', [source.UserData(source.Value,1)-(50/60000),...
    source.UserData(source.Value,2)+(50/60000)]);
end

function SelectTrlRip(source,event)
global listSel
listSel = 2;
kids = get(get(source, 'Parent'), 'Children');
tags = arrayfun(@(a)a.Tag, kids, 'uniformoutput',0);
axTag = strcmp(tags, 'Raw_Axes');
set(kids(axTag), 'xLim', [source.UserData(source.Value,1)-(50/60000),...
    source.UserData(source.Value,2)+(50/60000)]);
end

function NextRip(source,event)
global listSel
if listSel == 1
    lstTarg = 'ssnRip_Lst';
else
    lstTarg = 'trlRip_Lst';
end
kids = get(get(source, 'Parent'), 'Children');
tags = arrayfun(@(a)a.Tag, kids, 'uniformoutput',0);
lstTag = strcmp(tags, lstTarg);
kids(lstTag).Value = kids(lstTag).Value+1;
axTag = strcmp(tags, 'Raw_Axes');
set(kids(axTag), 'xLim', [kids(lstTag).UserData(kids(lstTag).Value,1)-(50),...
    kids(lstTag).UserData(kids(lstTag).Value,2)+(50)]);
end

function PrevRip(source,event)
global listSel
if listSel == 1
    lstTarg = 'ssnRip_Lst';
else
    lstTarg = 'trlRip_Lst';
end
kids = get(get(source, 'Parent'), 'Children');
tags = arrayfun(@(a)a.Tag, kids, 'uniformoutput',0);
lstTag = strcmp(tags, lstTarg);
kids(lstTag).Value = kids(lstTag).Value-1;
axTag = strcmp(tags, 'Raw_Axes');
set(kids(axTag), 'xLim', [kids(lstTag).UserData(kids(lstTag).Value,1)-(50),...
    kids(lstTag).UserData(kids(lstTag).Value,2)+(50)]);
end

function ZoomOut(source,event)
kids = get(get(source, 'Parent'), 'Children');
tags = arrayfun(@(a)a.Tag, kids, 'uniformoutput',0);
axTag = strcmp(tags, 'Raw_Axes');
curX = get(kids(axTag), 'xLim');
set(kids(axTag), 'xLim', [curX(1)-(10),...
    curX(2)+(10)]);
end

function ZoomIn(source,event)
kids = get(get(source, 'Parent'), 'Children');
tags = arrayfun(@(a)a.Tag, kids, 'uniformoutput',0);
axTag = strcmp(tags, 'Raw_Axes');
curX = get(kids(axTag), 'xLim');
set(kids(axTag), 'xLim', [curX(1)+(10),...
    curX(2)-(10)]);
end