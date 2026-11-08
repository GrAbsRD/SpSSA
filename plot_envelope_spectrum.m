function plot_envelope_spectrum(x, Fs)

x = x(:);
x = x - mean(x);

% Hilbert 包络
env = abs(hilbert(x));

% FFT
N = length(env);
Y = abs(fft(env));
f = (0:N-1) * Fs / N;

% 单边谱
idx = 1:floor(N/2);

plot(f(idx), Y(idx), 'LineWidth', 1.2);
% xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([5, 50]);
grid on;

end