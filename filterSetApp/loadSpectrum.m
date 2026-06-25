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

    % --- locate the data block: first line whose first token is a number ---
    % The delimiter is detected from the DATA lines only (not the prose header,
    % whose commas otherwise fool the detector), then the whole block is parsed
    % in one fast textscan call (str2double per token was ~100x slower on the
    % large theoretical-spectrum files of tens of thousands of rows).
    isData = ~cellfun('isempty', regexp(lines, '^\s*[+-]?(\d|\.\d)', 'once'));
    firstData = find(isData, 1);
    if isempty(firstData)
        error('loadSpectrum:noData', 'No numeric data parsed from %s', filename);
    end
    dlines = lines(firstData:end);

    samp = strjoin(dlines(1:min(40,numel(dlines))), newline);
    counts = [numel(strfind(samp, sprintf('\t'))), ...
              numel(strfind(samp, ',')), numel(strfind(samp, ';'))];
    [mx, di] = max(counts);
    delimRe = {'\t', ',', ';'};
    if mx == 0
        delimChar = {' ','\t'}; multi = true; splitRe = '\s+';
    else
        chars = {sprintf('\t'), ',', ';'};
        delimChar = chars{di}; multi = false; splitRe = delimRe{di};
    end
    ncol = numel(regexp(strtrim(dlines{1}), splitRe, 'split'));

    data = [];
    try
        C = textscan(strjoin(dlines, newline), repmat('%f',1,ncol), ...
            'CollectOutput',true, 'Delimiter',delimChar, ...
            'MultipleDelimsAsOne',multi, 'EndOfLine','\n');
        data = C{1};
    catch
        data = [];
    end
    if isempty(data) || size(data,1) < 0.5*numel(dlines)
        % robust fallback for ragged files: per-line numeric parse
        data = [];
        for i = 1:numel(dlines)
            vals = str2double(regexp(strtrim(dlines{i}), splitRe, 'split'));
            if numel(vals) >= 2 && isfinite(vals(1))
                data(end+1, 1:numel(vals)) = vals; %#ok<AGROW>
            end
        end
    end
    % keep only rows with a finite wavelength (drops any trailing footer)
    if ~isempty(data); data = data(isfinite(data(:,1)), :); end
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
