function selftest
% Headless sanity check of loadSpectrum + FilterEngine against real data.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
spectraRoot = fullfile(root, 'Spectra');
addpath(here);
lambda = (350:1:850)';

fprintf('--- loadSpectrum on varied formats ---\n');
tests = {
    'Proteins/GFP.txt'
    'Proteins/mNeonGreen.txt'
    'Proteins/iRFP670.txt'
    'Filters/FF01-520_44.txt'
    'Dichroics/Di01-R488_561.TXT'
    'Filters/T_600_lpxr.txt'
    'Illumations/UHP-F-470_Prizmatix.txt'};
for i = 1:numel(tests)
    S = loadSpectrum(fullfile(spectraRoot, tests{i}), lambda);
    fprintf('%-28s kind=%-11s max(ex)=%.2f max(em)=%.2f\n', ...
        tests{i}, S.kind, max(S.ex), max([S.em;0]));
end

fprintf('\n--- FilterEngine end-to-end (uSMAART-like) ---\n');
gfp = loadSpectrum(fullfile(spectraRoot,'Proteins/GFP.txt'),lambda);
ofp = loadSpectrum(fullfile(spectraRoot,'Proteins/cyOFP.txt'),lambda);
rub = loadSpectrum(fullfile(spectraRoot,'Proteins/mRuby3.txt'),lambda);
fluors = struct('name',{'GFP','cyOFP','mRuby3'}, ...
    'ex',{gfp.ex,ofp.ex,rub.ex},'em',{gfp.em,ofp.em,rub.em}, ...
    'brightness',{33.5,30.4,42.9});
lasers = struct('name',{'488','561'},'wl',{488,561},'power',{1,1});
prim = loadSpectrum(fullfile(spectraRoot,'Dichroics/Di01-R488_561.TXT'),lambda);
spl  = loadSpectrum(fullfile(spectraRoot,'Dichroics/FF562-Di03.TXT'),lambda);
fG   = loadSpectrum(fullfile(spectraRoot,'Filters/FF02-520_28.txt'),lambda);
fR   = loadSpectrum(fullfile(spectraRoot,'Filters/FF01_630_92.txt'),lambda);
ch(1) = struct('name','Green','emFilter',fG.ex, ...
    'path',struct('T',{prim.ex,spl.ex},'mode',{'T','R'}));
ch(2) = struct('name','Red','emFilter',fR.ex, ...
    'path',struct('T',{prim.ex,spl.ex},'mode',{'T','T'}));
assign = [1 3];
[S,eff] = FilterEngine.signalMatrix(fluors,lasers,ch,lambda);
CT = FilterEngine.crosstalkMatrix(S,assign);
disp('Signal matrix (rows=fluor, cols=Green/Red):'); disp(round(S,3));
disp('Crosstalk % :'); disp(round(100*CT,2));
fprintf('System score = %.3f\n', FilterEngine.systemScore(S,assign,5));

fprintf('\n--- detector QE weighting ---\n');
qe = 0.5*ones(numel(lambda),1);          % flat 50%% QE
S2 = FilterEngine.signalMatrix(fluors,lasers,ch,lambda,qe);
fprintf('signal ratio with 50%% QE (should be ~0.5): %.3f\n', S2(1,1)/S(1,1));

fprintf('\n--- laser back-reflection bleed (Rback=0.5%%) ---\n');
bleed = FilterEngine.laserBleed(ch, lasers, lambda, 0.005);
for k=1:2
  for j=1:2
    bg = bleed(k,j); sig = S(assign(k),k);
    od = -log10(max(bleed(k,j)/0.005, 1e-12));   % detection-path OD at laser line
    fprintf('  ch%d laser%dnm: bleed/signal=%.2e  path OD=%.1f\n', ...
        k, lasers(j).wl, bg/sig, od);
  end
end

fprintf('\n--- joint optimizer (makeChannels + optimizeJoint) ---\n');
spl2 = loadSpectrum(fullfile(spectraRoot,'Dichroics/FF564-Di01.txt'),lambda);
splitterCands = struct('name',{'FF562-Di03','FF564-Di01'},'ex',{spl.ex,spl2.ex});
fG2 = loadSpectrum(fullfile(spectraRoot,'Filters/FF02-529_24.txt'),lambda);
fGcands = struct('name',{'FF02-520_28','FF02-529_24'},'ex',{fG.ex,fG2.ex});
fRcands = struct('name',{'FF01_630_92'},'ex',{fR.ex});
R = FilterEngine.optimizeJoint(fluors,lasers,{'Green','Red'},prim.ex, ...
    splitterCands,{'R','T'},{fGcands,fRcands},assign,lambda,5);
fprintf('combinations ranked: %d\n', numel(R));
fprintf('best: splitter=%s  green=%s  red=%s  score=%.3f\n', ...
    R(1).dichroic, R(1).filters{1}, R(1).filters{2}, R(1).score);
fprintf('\nALL TESTS RAN OK\n');
end
