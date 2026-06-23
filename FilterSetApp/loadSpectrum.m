function S = loadSpectrum(filename, lambda)
%LOADSPECTRUM  Robust importer for optical spectra of mixed formats.
%   S = LOADSPECTRUM(FILENAME) reads a spectrum file and returns a struct.
%   S = LOADSPECTRUM(FILENAME, LAMBDA) resamples onto wavelength axis LAMBDA
%       (default 350:1:850 nm).
%
%   Auto-detects:
%     * delimiter (tab, comma, semicolon, whitespace)
%     * header row (present / absent / blank-first-cell)
%     * number of data columns (2 = lambda+T ; 3 = lambda+exc+em)
%     * percent (0-100) vs fractional (0-1) scaling for transmission-type data
%     * arbitrary wavelength step (0.2, 0.5, 1 nm ...) with duplicate handling
%
%   Output struct S:
%     .name    char, file base name
%     .lambda  resampled wavelength axis (column vector)
%     .ex      excitation / transmission column on lambda (col vector), or []
%     .em      emission column on lambda (col vector), or []
%     .kind    'fluorophore' (has ex & em) or 'filter' (single curve in .ex)
%     .file    full path
%
%   The single-curve case (filters, dichroics, lasers) is stored in .ex so
%   downstream code can always read S.ex as "the transmission/reflection".

    if nargin < 2 || isempty(lambda); lambda = (350:1:850)'; end
    lambda = lambda(:);

    raw = fileread(filename);

    lines = regexp(raw, '\r\n|\n|\r', 'split');
    lines = lines(~cellfun(@(s) all(isspace(s)) || isempty(s), lines));

    % --- choose the delimiter that actually parses the DATA, not the header ---
    % Detecting the delimiter from the first few lines is fooled by prose
    % headers (e.g. Semrock files whose preamble contains commas), so instead
    % parse with each candidate delimiter and keep whichever yields the most
    % numeric rows. A data row is one whose first cell is a finite number.
    delims = {'\t', ',', ';', '\s+'};
    data = []; bestRows = 0;
    for d = 1:numel(delims)
        D = []; rows = 0;
        for i = 1:numel(lines)
            toks = regexp(strtrim(lines{i}), delims{d}, 'split');
            vals = str2double(toks);
            if numel(vals) >= 2 && isfinite(vals(1))
                rows = rows + 1;
                D(rows, 1:numel(vals)) = vals; %#ok<AGROW>
            end
        end
        if rows > bestRows; bestRows = rows; data = D; end
    end
    if isempty(data)
        error('loadSpectrum:noData', 'No numeric data parsed from %s', filename);
    end

    wl = data(:,1);
    cols = data(:, 2:end);

    % Drop columns that are entirely NaN (e.g. empty emission column)
    cols(:, all(isnan(cols), 1)) = [];

    % Collapse duplicate / non-monotonic wavelengths (average)
    [wl, ~, ic] = unique(round(wl, 4), 'stable');
    if numel(wl) < size(cols,1)
        cols = splitapply(@(x) mean(x,1,'omitnan'), cols, ic);
    end
    [wl, order] = sort(wl);
    cols = cols(order, :);

    % Resample each column onto target axis
    res = nan(numel(lambda), size(cols,2));
    for c = 1:size(cols,2)
        good = isfinite(cols(:,c));
        if nnz(good) >= 2
            res(:,c) = interp1(wl(good), cols(good,c), lambda, 'linear', 0);
        end
    end
    res(~isfinite(res)) = 0;

    [~, base] = fileparts(filename);
    S.name   = base;
    S.lambda = lambda;
    S.file   = filename;

    if size(res,2) >= 2
        % fluorophore: excitation + emission, normalise each to peak 1
        S.kind = 'fluorophore';
        S.ex = normPeak(res(:,1));
        S.em = normPeak(res(:,2));
        S.floor = NaN;
    else
        % single transmission/reflection curve (filter/dichroic/laser)
        S.kind = 'filter';
        v = res(:,1);
        if max(v) > 1.5      % looks like percent
            v = v/100;
        end
        S.ex = min(max(v,0),1);   % clamp to [0,1]
        S.em = [];
        S.floor = measuredFloor(S.ex);
    end
end

function fl = measuredFloor(v)
%MEASUREDFLOOR Deepest blocking the curve actually demonstrates, used as a
%   realistic out-of-band floor for back-reflection estimates. Returns NaN
%   ("unknown") when the data shows no genuine deep blocking (e.g. linear-
%   scale spectra that round to 0), so the caller falls back to a global OD.
    pos = v(v > 0);
    if isempty(pos); fl = NaN; return; end
    mp = min(pos);
    if mp < 1e-4            % real high-dynamic-range OD data present
        fl = max(mp, 1e-7); % don't credit beyond OD7 from instrument noise
    else
        fl = NaN;           % no evidence of deep blocking -> unknown
    end
end

function v = normPeak(v)
    m = max(v);
    if m > 0; v = v/m; end
end
