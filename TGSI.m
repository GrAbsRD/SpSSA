function tgsi = TGSI(x, y, H, sigma)
% TGSI  Time-series to image similarity (luminance * structure)
% x,y   : real time series (vector)
% H     : image height (default 128)
% sigma : vertical spread in pixels (default 3)

    if nargin < 3 || isempty(H),     H = 128; end
    if nargin < 4 || isempty(sigma), sigma = 3; end

    % --- ensure column & same length ---
    x = x(:); y = y(:);
    N = min(numel(x), numel(y));
    x = x(1:N); y = y(1:N);

    % --- build expanded images (vectorized, Gaussian) ---
    Ix = ts2img_gauss_fast(x, H, sigma);
    Iy = ts2img_gauss_fast(y, H, sigma);

    % --- luminance & structure components (fast, global) ---
    L = 1;              % since images are normalized to [0,1]
    [lxy, sxy] = lum_struct_fast(Ix, Iy, L);

    tgsi = lxy * sxy;
end


function I = ts2img_gauss_fast(y, H, sigma)
% Vectorized: normalize -> row index r -> Gaussian vertical expansion

    y = double(y(:));
    ymin = min(y); ymax = max(y);
    y = (y - ymin) ./ (ymax - ymin + eps);     % [0,1]

    r = round(y * (H-1)) + 1;                  % T×1 in [1,H]

    v = (1:H)';                                % H×1
    D = v - r.';                               % H×T (implicit expansion)
    I = exp(-(D.^2) ./ (2*sigma^2));           % Gaussian

    % normalize to [0,1]
    I = I ./ (max(I(:)) + eps);
end


function [lxy, sxy] = lum_struct_fast(X, Y, L)
% Fast global luminance & structure (covariance) components

    x = double(X(:));
    y = double(Y(:));

    K1 = 0.01; K2 = 0.03;
    C1 = (K1*L)^2;
    C2 = (K2*L)^2 / 2;

    mux = mean(x);
    muy = mean(y);

    % Centered vectors
    xc = x - mux;
    yc = y - muy;

    % Use population moments (divide by N) for speed/consistency
    % (If you want unbiased, divide by (N-1) everywhere; ratio changes little.)
    N = numel(x);
    varx = (xc' * xc) / N;
    vary = (yc' * yc) / N;
    covxy = (xc' * yc) / N;

    sigx = sqrt(varx);
    sigy = sqrt(vary);

    lxy = (2*mux*muy + C1) / (mux^2 + muy^2 + C1);
    sxy = (covxy + C2) / (sigx*sigy + C2);
end