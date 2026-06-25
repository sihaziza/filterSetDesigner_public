function harmonize_proteins
%HARMONIZE_PROTEINS Convert spectra/Proteins/ into one file per protein, 3 columns.

appDir = fileparts(mfilename('fullpath'));
root = fileparts(appDir);
specRoot = fullfile(root, 'spectra');
if ~isfolder(specRoot)
    specRoot = fullfile(root, 'Spectra');
end
protDir = fullfile(specRoot, 'Proteins');
lambda = (350:850)';
backupDir = fullfile(root, ['Proteins_raw_backup_' datestr(now,'yyyymmdd_HHMMSS')]);
mkdir(backupDir);

fprintf('Backup: %s\n', backupDir);

% First convert split-folder spectra. These are preferred when present because
% they are the raw Abs/Em pair the user intentionally added.
dirs = dir(protDir);
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.','..'}));
for k = 1:numel(dirs)
    srcDir = fullfile(dirs(k).folder, dirs(k).name);
    files = dir(fullfile(srcDir, '*.txt'));
    names = lower({files.name});
    absIdx = find(contains(names, 'abs') | contains(names, 'ex'), 1);
    emIdx = find(contains(names, 'em'), 1);
    if isempty(absIdx) || isempty(emIdx)
        warning('Skipping split folder without Abs/Em pair: %s', srcDir);
        continue;
    end
    ex = readSingleCurve(fullfile(files(absIdx).folder, files(absIdx).name), lambda);
    em = readSingleCurve(fullfile(files(emIdx).folder, files(emIdx).name), lambda);
    outFile = fullfile(protDir, [dirs(k).name '.txt']);
    copyExisting(outFile, backupDir);
    writeProtein(outFile, lambda, normPeak(ex), normPeak(em));
    movefile(srcDir, fullfile(backupDir, dirs(k).name));
    fprintf('folder -> %s\n', outFile);
end

% Then normalize remaining top-level text files that were not produced from a
% split folder. CSV duplicates are archived so the app sees one file per entry.
txtFiles = dir(fullfile(protDir, '*.txt'));
for k = 1:numel(txtFiles)
    fp = fullfile(txtFiles(k).folder, txtFiles(k).name);
    copyExisting(fp, backupDir);
    try
        S = loadSpectrum(fp, lambda);
        ex = S.ex;
        em = S.em;
        if isempty(em)
            em = zeros(size(lambda));
        end
        writeProtein(fp, lambda, normPeak(ex), normPeak(em));
        fprintf('normalized %s\n', fp);
    catch ME
        warning('Could not normalize %s: %s', fp, ME.message);
    end
end

csvFiles = dir(fullfile(protDir, '*.csv'));
for k = 1:numel(csvFiles)
    movefile(fullfile(csvFiles(k).folder, csvFiles(k).name), ...
        fullfile(backupDir, csvFiles(k).name));
    fprintf('archived duplicate CSV %s\n', csvFiles(k).name);
end

fprintf('Done. Proteins now uses one top-level TXT per entry on 350:850 nm.\n');
end

function copyExisting(fp, backupDir)
if exist(fp, 'file')
    [~, base, ext] = fileparts(fp);
    copyfile(fp, fullfile(backupDir, [base ext]));
end
end

function y = readSingleCurve(fp, lambda)
raw = fileread(fp);
lines = regexp(raw, '\r\n|\n|\r', 'split');
data = [];
for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line) || isempty(regexp(line, '^[-+]?\d', 'once'))
        continue;
    end
    toks = regexp(line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
    vals = str2double(toks);
    if numel(vals) >= 2 && isfinite(vals(1)) && isfinite(vals(2))
        data(end+1,:) = vals(1:2); %#ok<AGROW>
    end
end
if isempty(data)
    error('No numeric wavelength/value rows parsed from %s', fp);
end
wl = data(:,1);
v = data(:,2);
[wl, ~, ic] = unique(round(wl, 4), 'stable');
if numel(wl) < numel(v)
    v = splitapply(@(x) mean(x, 'omitnan'), v, ic);
end
[wl, order] = sort(wl);
v = v(order);
good = isfinite(wl) & isfinite(v);
y = interp1(wl(good), v(good), lambda, 'linear', 0);
y(~isfinite(y) | y < 0) = 0;
end

function writeProtein(fp, lambda, ex, em)
fid = fopen(fp, 'w');
if fid < 0
    error('Could not write %s', fp);
end
c = onCleanup(@() fclose(fid));
fprintf(fid, 'Wavelength_nm\tExcitation\tEmission\n');
for i = 1:numel(lambda)
    fprintf(fid, '%d\t%.8g\t%.8g\n', lambda(i), ex(i), em(i));
end
end

function v = normPeak(v)
v = v(:);
v(~isfinite(v) | v < 0) = 0;
m = max(v);
if m > 0
    v = v / m;
end
end
