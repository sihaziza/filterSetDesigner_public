function s = autofluorPreset(name, lambda)
%AUTOFLUORPRESET  Parametric autofluorescence sources for snrModel.
%   Each source has an excitation-wavelength-dependent absorption (.absorb,
%   peak 1) and a broad emission spectrum (.em, peak 1), plus a default
%   strength (peak absorbed fraction ~ ε·c·d) and quantum yield.
%
%   Presets capture the well-known behaviour that endogenous/optical
%   autofluorescence is excited far more strongly at short wavelengths and
%   emits as a broad, red-tailed band:
%     'Brain tissue'  - NADH/flavin-like, strong 400-500 nm exc, em ~ 500 nm
%     'Silica fiber'  - patch-cord/fiber AF, weak broad blue-green
%     'Flavin (FAD)'  - peaked ~450 nm exc, em ~ 530 nm
%     'Lipofuscin'    - very broad, red-shifted, excitable to ~560 nm
%   Use 'list' to get the available names.

    lambda = lambda(:);
    if nargin >= 1 && strcmpi(name,'list')
        s = {'Brain tissue','Silica fiber','Flavin (FAD)','Lipofuscin'}; return;
    end

    switch lower(name)
        case 'brain tissue'
            absorb = expFalloff(lambda, 400, 110);          % strong blue, decays to red
            em     = lognormalBand(lambda, 510, 90);        % broad, red-tailed
            s = pack('Brain tissue', absorb, em, 2e-3, 0.10);
        case 'silica fiber'
            absorb = expFalloff(lambda, 400, 80);
            em     = lognormalBand(lambda, 470, 80);
            s = pack('Silica fiber', absorb, em, 5e-4, 0.05);
        case 'flavin (fad)'
            absorb = gaussBand(lambda, 450, 45);
            em     = lognormalBand(lambda, 530, 70);
            s = pack('Flavin (FAD)', absorb, em, 1e-3, 0.10);
        case 'lipofuscin'
            absorb = expFalloff(lambda, 400, 160);          % excitable well into green
            em     = lognormalBand(lambda, 600, 110);
            s = pack('Lipofuscin', absorb, em, 1e-3, 0.08);
        otherwise
            % flat default custom source
            s = pack(name, gaussBand(lambda,480,40), lognormalBand(lambda,540,90), 1e-3, 0.1);
    end
end

function s = pack(name, absorb, em, strength, qy)
    s = struct('name',name,'absorb',absorb,'em',em,'strength',strength,'qy',qy);
end

function v = expFalloff(lambda, lam0, tau)
    v = exp(-(lambda - lam0)/tau); v(lambda < lam0) = 1; v = v/max(v);
end

function v = gaussBand(lambda, mu, sig)
    v = exp(-0.5*((lambda-mu)/sig).^2); v = v/max(v);
end

function v = lognormalBand(lambda, peak, width)
    % skewed band with a red tail, peaked at 'peak' nm
    b = 0.25;                                   % skew
    x = 1 + 2*b*(lambda-peak)/width;
    v = zeros(size(lambda)); ok = x > 0;
    v(ok) = exp(-log(x(ok)).^2 / (2*b^2));
    v = v/max([v;eps]);
end
