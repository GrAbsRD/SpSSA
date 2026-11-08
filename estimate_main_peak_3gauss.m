function [fmax, delta_f, model, info] = estimate_main_peak_3gauss(x, fs)
%ESTIMATE_MAIN_PEAK_3GAUSS
% Robust estimation of dominant frequency and bandwidth.
% Fixed to prevent NaN/Inf errors in lsqcurvefit.

    %% -------- 1. PSD estimation --------
    x = x(:) - mean(x);
    
    % Use pwelch or periodogram. Ensure vector length is sufficient.
    nfft = max(256, 2^nextpow2(length(x)));
    [PSD, f] = pwelch(x, [], [], nfft, fs);
    
    % Safety check: Zero Signal
    max_psd = max(PSD);
    if max_psd == 0 || isnan(max_psd)
        % Fallback for silence/error
        fmax = 0;
        delta_f = 1.0; 
        model = struct('A1',0, 'sigma1',1, 'mu1',0);
        info = struct('fail', true);
        return;
    end
    
    PSD = PSD / max_psd;   % Normalize

    %% -------- 2. Find dominant and secondary peaks --------
    % Find peaks with a minimum height to ignore noise
    [pks, locs] = findpeaks(PSD, f, 'SortStr','descend', 'NPeaks', 2);

    if isempty(locs)
        % If findpeaks fails (e.g., pure decay or DC), use max value
        [~, idx_max] = max(PSD);
        fmax = f(idx_max);
        pks = PSD(idx_max);
    else
        fmax = locs(1);
    end

    if length(locs) >= 2
        f2 = locs(2);
    else
        % Fallback f2 logic
        if fmax > 0
            f2 = fmax * 1.5; 
        else
            f2 = fs * 0.1; % Arbitrary if fmax is DC
        end
    end
    
    % Ensure f2 is not exactly fmax to avoid sigma3=0
    if abs(f2 - fmax) < 1e-6
         f2 = fmax + (f(2)-f(1))*2;
    end

    f3 = (fmax + f2) / 2;
    
    % Frequency resolution for lower bounds
    df = f(2) - f(1);

    %% -------- 3. Initial parameter guesses --------
    % Helper to safely get amplitude
    get_A = @(freq) max(1e-6, interp1(f, PSD, freq, 'linear', 0));
    
    A1 = 0.5 * get_A(fmax);
    A2 = 0.5 * get_A(f2);
    A3 = 0.25 * get_A(f3);

    % Robust Sigma estimation
    sigma1 = find_sigma_halfpower(f, PSD, fmax, df);
    sigma2 = find_sigma_halfpower(f, PSD, f2, df);
    sigma3 = max(df, 4 * abs(fmax - f2)); % Ensure positive

    p0 = [A1 sigma1 A2 sigma2 A3 sigma3];
    
    % Final NaN check for p0
    if any(isnan(p0)) || any(isinf(p0))
        p0(isnan(p0)) = 0.1;
        p0(isinf(p0)) = 10;
    end

    %% -------- 4. Nonlinear least squares fitting --------
    % Gaussian model
    gauss3 = @(p, f_val) ...
        p(1)*exp(-(f_val-fmax).^2./(2*p(2)^2 + 1e-9)) + ... % Add epsilon to denom
        p(3)*exp(-(f_val-f2  ).^2./(2*p(4)^2 + 1e-9)) + ...
        p(5)*exp(-(f_val-f3  ).^2./(2*p(6)^2 + 1e-9));

    % Set Constraints: Amplitudes > 0, Sigmas >= df/2
    lb = [0,    df/2, 0,    df/2, 0,    df/2];
    ub = [inf,  inf,  inf,  inf,  inf,  inf];

    opts = optimoptions('lsqcurvefit', ...
                        'Display','off', ...
                        'MaxIterations', 200, ...
                        'Algorithm', 'trust-region-reflective'); % Safer than Levenberg-Marquardt for bounds

    try
        p = lsqcurvefit(gauss3, p0, f, PSD, lb, ub, opts);
    catch
        % If fit fails, fall back to initial guess
        p = p0;
    end

    %% -------- 5. Extract bandwidth --------
    sigma1_hat = p(2);
    delta_f = 2.5 * sigma1_hat;

    %% -------- 6. Output structures --------
    model = struct();
    model.A1 = p(1); model.sigma1 = p(2); model.mu1 = fmax;
    % ... fill others if needed ...

    info = struct();
    info.f = f;
    info.PSD = PSD;
    info.p = p;

end

function sigma = find_sigma_halfpower(f, PSD, f0, df)
% Estimate sigma from half-power width with safety clamping
    
    try
        P0 = interp1(f, PSD, f0, 'linear', 0);
    catch
        P0 = 0;
    end
    
    if P0 == 0
        sigma = max(df, f0 * 0.1);
        return;
    end

    target = 2/3 * P0;
    idx = find(PSD >= target);

    if isempty(idx)
        sigma = max(df, f0 * 0.1);
        return
    end

    f_low = f(min(idx));
    f_high = f(max(idx));
    
    width = f_high - f_low;
    
    % CRITICAL FIX: If width is 0 (single point peak), enforce min width
    if width < df
        width = df; 
    end

    sigma = width / 2.355; 
end