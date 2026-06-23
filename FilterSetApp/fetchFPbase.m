function S = fetchFPbase(query, lambda)
%FETCHFPBASE  Download a fluorophore spectrum from the FPbase web API.
%   S = FETCHFPBASE(QUERY) searches FPbase (https://www.fpbase.org) for a
%   fluorescent protein / dye by name and returns a spectrum struct in the
%   same format as LOADSPECTRUM, including web-sourced brightness.
%
%   S = FETCHFPBASE(QUERY, LAMBDA) resamples onto wavelength axis LAMBDA.
%
%   Requires an internet connection. Uses FPbase's public REST endpoints:
%     /api/proteins/?name__icontains=...   (search + photophysics)
%     /spectra/?q=...                       (spectra as csv/json)
%
%   Fields added beyond loadSpectrum:
%     .exMax, .emMax  peak wavelengths (nm)
%     .ec             extinction coefficient (M^-1 cm^-1)
%     .qy             quantum yield
%     .brightness     ec*qy/1000  (relative units, like FPbase "brightness")

    if nargin < 2 || isempty(lambda); lambda = (350:1:850)'; end
    lambda = lambda(:);

    opts = weboptions('Timeout', 20, 'ContentType', 'json');

    % --- 1) find the protein record ---
    base = 'https://www.fpbase.org';
    url  = sprintf('%s/api/proteins/?name__icontains=%s&format=json', ...
                   base, urlencode(query));
    recs = webread(url, opts);
    if isempty(recs)
        error('fetchFPbase:notFound', 'No FPbase match for "%s".', query);
    end
    if iscell(recs); recs = [recs{:}]; end
    names = {recs.name};
    exact = find(strcmpi(names, query), 1);
    if isempty(exact); exact = 1; end
    rec = recs(exact);

    S.name   = rec.name;
    S.lambda = lambda;
    S.file   = url;
    S.kind   = 'fluorophore';
    state = firstState(rec);
    S.exMax  = getfieldsafe(state, 'ex_max', getfieldsafe(rec, 'ex_max', NaN));
    S.emMax  = getfieldsafe(state, 'em_max', getfieldsafe(rec, 'em_max', NaN));
    S.ec     = getfieldsafe(state, 'ext_coeff', getfieldsafe(rec, 'ext_coeff', NaN));
    S.qy     = getfieldsafe(state, 'qy', getfieldsafe(rec, 'qy', NaN));
    if isfinite(S.ec) && isfinite(S.qy)
        S.brightness = S.ec * S.qy / 1000;
    else
        S.brightness = NaN;
    end

    % --- 2) fetch the spectra for this fluorophore ---
    slug = getfieldsafe(rec, 'slug', lower(strrep(rec.name,' ','-')));
    surl = sprintf('%s/api/proteins/spectra/?slug=%s&format=json', base, urlencode(slug));
    try
        sp = webread(surl, opts);
    catch
        sp = [];
    end
    sp = matchingSpectrumRecord(sp, slug);

    ex = zeros(numel(lambda),1); em = zeros(numel(lambda),1);
    if ~isempty(sp)
        if iscell(sp); sp = sp{1}; end
        if isfield(sp,'spectra')
            spectra = sp.spectra;
            if ~iscell(spectra); spectra = num2cell(spectra); end
            for k = 1:numel(spectra)
                s = spectra{k};
                d = spectrumData(s);
                w = d(:,1); v = d(:,2);
                vi = interp1(w, v, lambda, 'linear', 0);
                vi = vi / max([max(vi) eps]);
                switch spectrumSubtype(s)
                    case {'ex','ab','abs'}; ex = max(ex, vi);
                    case {'em'};            em = max(em, vi);
                end
                if isfield(s,'max') && isfinite(s.max)
                    switch spectrumSubtype(s)
                        case {'ex','ab','abs'}; S.exMax = s.max;
                        case {'em'};            S.emMax = s.max;
                    end
                end
                if isfield(s,'ec') && isfinite(s.ec); S.ec = s.ec; end
                if isfield(s,'qy') && isfinite(s.qy); S.qy = s.qy; end
            end
        end
    end
    S.ex = ex; S.em = em;
    if isfinite(S.ec) && isfinite(S.qy)
        S.brightness = S.ec * S.qy / 1000;
    end

    if ~any(em) && ~any(ex)
        warning('fetchFPbase:noSpectra', ...
            'Record found for "%s" but no spectral data returned.', S.name);
    end
end

function v = getfieldsafe(s, f, d)
    if isstruct(s) && isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

function st = firstState(rec)
    st = struct();
    if ~isstruct(rec) || ~isfield(rec,'states') || isempty(rec.states); return; end
    states = rec.states;
    if iscell(states); states = [states{:}]; end
    st = states(1);
end

function sp = matchingSpectrumRecord(sp, slug)
    if isempty(sp); return; end
    if iscell(sp); sp = [sp{:}]; end
    if numel(sp) > 1 && isfield(sp,'slug')
        idx = find(strcmpi({sp.slug}, slug), 1);
        if ~isempty(idx); sp = sp(idx); else; sp = sp(1); end
    end
end

function d = spectrumData(s)
    d = s.data;
    if iscell(d); d = cell2mat(d); end
    d = squeeze(d);
    if size(d,2) < 2 && size(d,1) >= 2; d = d'; end
    if size(d,2) < 2; error('FPbase spectrum data is not wavelength/value pairs.'); end
    d = d(:,1:2);
end

function sub = spectrumSubtype(s)
    sub = lower(char(string(getfieldsafe(s,'subtype',''))));
    if ~isempty(sub); return; end
    state = lower(char(string(getfieldsafe(s,'state',''))));
    if contains(state, '_em') || endsWith(state, 'em')
        sub = 'em';
    elseif contains(state, '_ex') || contains(state, '_ab') || endsWith(state, 'ex')
        sub = 'ex';
    else
        sub = '';
    end
end
