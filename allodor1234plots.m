%% Plot ALL four odors (1,2,3,4) together — Combined & Aligned (c_al)
% This script reconstructs aligned matrices for Odor 1 and 3 using the
% already-computed alignment indices for Odor 2 and 4 (NoGo) in your files.
% Odor pairing: (1↔2) and (3↔4).

clear; clc;

% --------------------- Settings ---------------------
condToPlot = 'c_al';       % we reconstruct c_al for all four odors
nNan       = 5000;         % canvas width for aligned matrices (matches your pipeline)
minSessionsForPlot = 1;    % require at least N sessions contributing

% Colors: keep user’s preferred for 2 (NoStim) and 4 (Stim), add distinct hues for 1 & 3
col_o2 = [91 141 184]/255;   % Odor 2 (NoStim)  — user preference
col_o4 = [229 115 115]/255;  % Odor 4 (Stim)    — user preference
col_o1 = [100 180 170]/255;  % Odor 1 (Go; paired with 2)
col_o3 = [200 170  90]/255;  % Odor 3 (Go; paired with 4)

rt_ylim   = [0 2];         % reaction-time axis range (s)
lc_ylim   = [0 1.0];       % learning-curve axis range (fraction correct)
lw_mean   = 3;             % line width for mean
lw_shaded = 1;             % shadedErrorBar line width
% ---------------------------------------------------

% Check for shadedErrorBar; if not available, we fall back to simple mean±SEM lines.
useShaded = exist('shadedErrorBar','file') == 2;

dataFiles = dir('*_learning_Stim_Feedback.mat');
if isempty(dataFiles)
    error('No *_learning_Stim_Feedback.mat files found in this folder.');
end

for fileIdx = 1:numel(dataFiles)
    fname = dataFiles(fileIdx).name;
    S = load(fname);

    % Get animalID from filename prefix
    [~, base, ~] = fileparts(fname);
    parts = strsplit(base, '_');
    animalID = parts{1};

    % Pull allodors (has learning_data (col 9), ncond (col 10), rt (col 11), pass flag (col 13))
    eval(['allodors = S.' animalID '_allodors;']);

    % Pre-allocate aligned matrices for all 4 odors
    lc_al_1 = nan(size(allodors,1), nNan);
    lc_al_2 = nan(size(allodors,1), nNan);
    lc_al_3 = nan(size(allodors,1), nNan);
    lc_al_4 = nan(size(allodors,1), nNan);

    rt_al_1 = nan(size(allodors,1), nNan);
    rt_al_2 = nan(size(allodors,1), nNan);
    rt_al_3 = nan(size(allodors,1), nNan);
    rt_al_4 = nan(size(allodors,1), nNan);

    % Rebuild aligned matrices using the final learned indices nval2/nval4
    for i = 1:size(allodors,1)
        % Must have passed combined learning criteria (column 13 == 1)
        if ~(isfield(allodors, 'data') || true) %#ok<*NOPRT>
            % no-op; to quiet MATLAB code analyzer in some editors
        end
        if isempty(allodors{i,13}) || allodors{i,13} ~= 1
            continue
        end

        % Learning traces (sliding mean) for odors 1..4
        lc1 = allodors{i,9}{1}(:,2);  lc2 = allodors{i,9}{2}(:,2);
        lc3 = allodors{i,9}{3}(:,2);  lc4 = allodors{i,9}{4}(:,2);

        % Reaction times for odors 1..4 (seconds)
        rt1 = allodors{i,11}{1}(:,2); rt2 = allodors{i,11}{2}(:,2);
        rt3 = allodors{i,11}{3}(:,2); rt4 = allodors{i,11}{4}(:,2);

        % Learned indices for NoGo odors (computed & stored by your pipeline)
        nval2 = allodors{i,10}{2};    % alignment for pair (1,2)
        nval4 = allodors{i,10}{4};    % alignment for pair (3,4)

        % Require both to exist to keep consistent alignment canvas
        if isempty(nval2) || isempty(nval4)
            continue
        end

        % Place pairs (1,2) aligned to nval2
        s2 = nNan/2 - nval2 + 1;
        e2 = s2 + numel(lc2) - 1;
        if s2 >= 1 && e2 <= nNan
            lc_al_2(i, s2:e2) = lc2;
            rt_al_2(i, s2:e2) = rt2;

            % Align Odor 1 using the SAME anchor as Odor 2 (pair)
            s1 = s2; e1 = s1 + numel(lc1) - 1;
            if s1 >= 1 && e1 <= nNan
                lc_al_1(i, s1:e1) = lc1;
                rt_al_1(i, s1:e1) = rt1;
            end
        end

        % Place pairs (3,4) aligned to nval4
        s4 = nNan/2 - nval4 + 1;
        e4 = s4 + numel(lc4) - 1;
        if s4 >= 1 && e4 <= nNan
            lc_al_4(i, s4:e4) = lc4;
            rt_al_4(i, s4:e4) = rt4;

            % Align Odor 3 using the SAME anchor as Odor 4 (pair)
            s3 = s4; e3 = s3 + numel(lc3) - 1;
            if s3 >= 1 && e3 <= nNan
                lc_al_3(i, s3:e3) = lc3;
                rt_al_3(i, s3:e3) = rt3;
            end
        end
    end

    % Determine common plotting window across all four odors
    cond1 = sum(~isnan(lc_al_1), 1);
    cond2 = sum(~isnan(lc_al_2), 1);
    cond3 = sum(~isnan(lc_al_3), 1);
    cond4 = sum(~isnan(lc_al_4), 1);
    condAll = cond1 + cond2 + cond3 + cond4;

    ninit = find(condAll >= minSessionsForPlot, 1, 'first');
    nfin  = find(condAll >= minSessionsForPlot, 1, 'last');

    if isempty(ninit) || isempty(nfin) || nfin <= ninit
        warning('No overlapping window for %s — skipping.', animalID);
        continue
    end

    % --------- Plot Learning Curves (1,2,3,4) ----------
    figure('Name', [animalID ' — Learning Curves (1–4), c\_al'], 'Color', 'w');

    hold on;
    if useShaded
        shadedErrorBar(ninit:nfin, nanmean(lc_al_1(:,ninit:nfin),1), nanstd(lc_al_1(:,ninit:nfin),0,1)./sqrt(max(1,cond1(ninit:nfin))), ...
            'lineProps', {'Color', col_o1, 'LineWidth', lw_shaded});
        shadedErrorBar(ninit:nfin, nanmean(lc_al_2(:,ninit:nfin),1), nanstd(lc_al_2(:,ninit:nfin),0,1)./sqrt(max(1,cond2(ninit:nfin))), ...
            'lineProps', {'Color', col_o2, 'LineWidth', lw_shaded});
        shadedErrorBar(ninit:nfin, nanmean(lc_al_3(:,ninit:nfin),1), nanstd(lc_al_3(:,ninit:nfin),0,1)./sqrt(max(1,cond3(ninit:nfin))), ...
            'lineProps', {'Color', col_o3, 'LineWidth', lw_shaded});
        shadedErrorBar(ninit:nfin, nanmean(lc_al_4(:,ninit:nfin),1), nanstd(lc_al_4(:,ninit:nfin),0,1)./sqrt(max(1,cond4(ninit:nfin))), ...
            'lineProps', {'Color', col_o4, 'LineWidth', lw_shaded});
    else
        % Fallback: mean lines (no shadedErrorBar)
        plot(ninit:nfin, nanmean(lc_al_1(:,ninit:nfin),1), '-', 'Color', col_o1, 'LineWidth', lw_mean);
        plot(ninit:nfin, nanmean(lc_al_2(:,ninit:nfin),1), '-', 'Color', col_o2, 'LineWidth', lw_mean);
        plot(ninit:nfin, nanmean(lc_al_3(:,ninit:nfin),1), '-', 'Color', col_o3, 'LineWidth', lw_mean);
        plot(ninit:nfin, nanmean(lc_al_4(:,ninit:nfin),1), '-', 'Color', col_o4, 'LineWidth', lw_mean);
    end
    hold off;

    ylim([-0.2 1.2]);
    xlim([ninit nfin]);
    xlabel('Aligned trial index');
    ylabel('Performance (sliding mean)');
    title(['Learning Curves (c\_al) — ' animalID ' — Odors 1–4']);
    legend({'Odor 1 (aligned to Odor 2)', 'Odor 2 (NoStim)', ...
            'Odor 3 (aligned to Odor 4)', 'Odor 4 (Stim)'}, ...
            'Location', 'southeast');

    % --------- Plot Reaction Times (1,2,3,4) ----------
    figure('Name', [animalID ' — Reaction Times (1–4), c\_al'], 'Color', 'w');

    hold on;
    if useShaded
        shadedErrorBar(ninit:nfin, nanmean(rt_al_1(:,ninit:nfin),1), nanstd(rt_al_1(:,ninit:nfin),0,1)./sqrt(max(1,cond1(ninit:nfin))), ...
            'lineProps', {'Color', col_o1, 'LineWidth', lw_shaded});
        shadedErrorBar(ninit:nfin, nanmean(rt_al_2(:,ninit:nfin),1), nanstd(rt_al_2(:,ninit:nfin),0,1)./sqrt(max(1,cond2(ninit:nfin))), ...
            'lineProps', {'Color', col_o2, 'LineWidth', lw_shaded});
        shadedErrorBar(ninit:nfin, nanmean(rt_al_3(:,ninit:nfin),1), nanstd(rt_al_3(:,ninit:nfin),0,1)./sqrt(max(1,cond3(ninit:nfin))), ...
            'lineProps', {'Color', col_o3, 'LineWidth', lw_shaded});
        shadedErrorBar(ninit:nfin, nanmean(rt_al_4(:,ninit:nfin),1), nanstd(rt_al_4(:,ninit:nfin),0,1)./sqrt(max(1,cond4(ninit:nfin))), ...
            'lineProps', {'Color', col_o4, 'LineWidth', lw_shaded});
    else
        plot(ninit:nfin, nanmean(rt_al_1(:,ninit:nfin),1), '-', 'Color', col_o1, 'LineWidth', lw_mean);
        plot(ninit:nfin, nanmean(rt_al_2(:,ninit:nfin),1), '-', 'Color', col_o2, 'LineWidth', lw_mean);
        plot(ninit:nfin, nanmean(rt_al_3(:,ninit:nfin),1), '-', 'Color', col_o3, 'LineWidth', lw_mean);
        plot(ninit:nfin, nanmean(rt_al_4(:,ninit:nfin),1), '-', 'Color', col_o4, 'LineWidth', lw_mean);
    end
    hold off;

    ylim(rt_ylim);
    xlim([ninit nfin]);
    xlabel('Aligned trial index');
    ylabel('Reaction time (s)');
    title(['Reaction Times (c\_al) — ' animalID ' — Odors 1–4']);
    legend({'Odor 1 (aligned to Odor 2)', 'Odor 2 (NoStim)', ...
            'Odor 3 (aligned to Odor 4)', 'Odor 4 (Stim)'}, ...
            'Location', 'northeast');
end

disp('Done: Plotted Odors 1–4 (combined & aligned) for all animals.');

%%

%% Plot ALL individual sessions (no averaging) for each animal
% Odors 1–4 plotted separately for each session, aligned as before.

for fileIdx = 1:numel(dataFiles)
    fname = dataFiles(fileIdx).name;
    S = load(fname);

    % Extract animalID
    [~, base, ~] = fileparts(fname);
    parts = strsplit(base, '_');
    animalID = parts{1};

    eval(['allodors = S.' animalID '_allodors;']);

    figure('Name', [animalID ' — Individual Sessions (Odors 1–4)'], 'Color', 'w');
    hold on;

    for i = 1:size(allodors,1)
        % Only plot if this session passed criteria
        if isempty(allodors{i,13}) || allodors{i,13} ~= 1
            continue
        end

        % Grab learning curves
        lc1 = allodors{i,9}{1}(:,2);
        lc2 = allodors{i,9}{2}(:,2);
        lc3 = allodors{i,9}{3}(:,2);
        lc4 = allodors{i,9}{4}(:,2);

        % Alignment anchors
        nval2 = allodors{i,10}{2};
        nval4 = allodors{i,10}{4};

        if isempty(nval2) || isempty(nval4), continue; end

        % Align Odor 1/2
        s2 = nNan/2 - nval2 + 1;
        e2 = s2 + numel(lc2) - 1;
        plot(s2:e2, lc2, '-', 'Color', col_o2, 'LineWidth', 1);
        s1 = s2; e1 = s1 + numel(lc1) - 1;
        plot(s1:e1, lc1, '-', 'Color', col_o1, 'LineWidth', 1);

        % Align Odor 3/4
        s4 = nNan/2 - nval4 + 1;
        e4 = s4 + numel(lc4) - 1;
        plot(s4:e4, lc4, '-', 'Color', col_o4, 'LineWidth', 1);
        s3 = s4; e3 = s3 + numel(lc3) - 1;
        plot(s3:e3, lc3, '-', 'Color', col_o3, 'LineWidth', 1);
    end

    ylim([-0.2 1.2]);
    xlabel('Aligned trial index');
    ylabel('Performance (sliding mean)');
    title(['All Individual Sessions — ' animalID]);
    legend({'Odor 1','Odor 2','Odor 3','Odor 4'}, 'Location','southeast');
    hold off;
end

