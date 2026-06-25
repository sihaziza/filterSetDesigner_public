function qe = detectorPreset(name, lambda)
%DETECTORPRESET  Built-in detector quantum-efficiency / responsivity curves
%   (0-1 vs wavelength) so the schematic detector dropdown is usable without a
%   measured detector file. Approximate, typical-datasheet shapes.
%
%   Names:  'Si APD (typical)'   - silicon APD, peak ~650-700 nm (e.g. APD440A2)
%           'sCMOS (typical)'    - back-illuminated sCMOS, flat ~85% 500-650 nm
%           'EMCCD (typical)'    - back-illuminated EMCCD, ~90% 550-650 nm
%           'GaAsP PMT (typical)'- peak ~40% near 500 nm
%           'Bialkali PMT (typical)' - peak ~28% near 400 nm
%   Use 'list' to get the available names.

    lambda = lambda(:);
    if nargin >= 1 && strcmpi(name,'list')
        qe = {'Si APD (typical)','sCMOS (typical)','EMCCD (typical)', ...
              'GaAsP PMT (typical)','Bialkali PMT (typical)'};
        return;
    end
    switch lower(name)
        case 'si apd (typical)'
            qe = 0.80 * gauss(lambda, 680, 170);
            qe = qe .* logistic(lambda, 360, 18);          % UV roll-off
        case 'scmos (typical)'
            qe = 0.85 * superg(lambda, 560, 170, 4);
        case 'emccd (typical)'
            qe = 0.92 * superg(lambda, 600, 150, 4);
        case 'gaasp pmt (typical)'
            qe = 0.45 * gauss(lambda, 500, 95) .* logistic(lambda, 350, 15);
        case 'bialkali pmt (typical)'
            qe = 0.28 * gauss(lambda, 400, 80) .* logistic(lambda, 300, 15);
        otherwise
            qe = ones(numel(lambda),1);                    % unknown -> ideal
    end
    qe = min(1, max(0, qe));
end

function v = gauss(lam, mu, sig);  v = exp(-0.5*((lam-mu)/sig).^2); end
function v = superg(lam, mu, sig, p); v = exp(-((lam-mu)/sig).^(2*p)); end
function v = logistic(lam, x0, k);  v = 1./(1+exp(-(lam-x0)/k)); end
