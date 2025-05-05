% Вектор частот для анализа
f_values = 0.1:0.1:2; % Частота от 0.1 до 2 Гц с шагом 0.1
mse_values = zeros(size(f_values)); % Для хранения MSE

% Параметры
dt = 0.01; % Шаг дискретизации
t = 0:dt:100; % Временной вектор (100 секунд)
fs = 1/dt; % Частота дискретизации
window_size = 10; % Размер окна для усреднения

% Цикл по частотам
for i = 1:length(f_values)
    f = f_values(i); % Текущая частота

    % 1. Генерация случайного сигнала
    signal = randn(1, length(t)); % Белый шум

    % 2. Интерполяция
    dt_interp = dt/2; % Удвоенная частота дискретизации
    t_interp = 0:dt_interp:100;
    signal_interp = interp1(t, signal, t_interp, 'linear');

    % 3. Усреднение
    impact = conv(signal_interp, ones(1, window_size)/window_size, 'valid');

    % 4. Измерение MSE
    mse_values(i) = mean((signal_interp(1:length(impact)) - impact).^2);
end

% Спектральный анализ для одной частоты (например, f = 0.5 Гц)
f = 0.5;
signal = randn(1, length(t));
t_interp = 0:dt_interp:100;
signal_interp = interp1(t, signal, t_interp, 'linear');
impact = conv(signal_interp, ones(1, window_size)/window_size, 'valid');

% FFT
N = length(signal);
N_interp = length(signal_interp);
N_impact = length(impact);
fs_interp = 1/dt_interp;

signal_fft = fft(signal);
freq = (0:N-1)*(fs/N);
signal_psd = abs(signal_fft(1:floor(N/2)+1)).^2 / N;
freq = freq(1:floor(N/2)+1);

signal_interp_fft = fft(signal_interp);
freq_interp = (0:N_interp-1)*(fs_interp/N_interp);
signal_interp_psd = abs(signal_interp_fft(1:floor(N_interp/2)+1)).^2 / N_interp;
freq_interp = freq_interp(1:floor(N_interp/2)+1);

impact_fft = fft(impact, N_interp);
impact_psd = abs(impact_fft(1:floor(N_interp/2)+1)).^2 / N_interp;
freq_impact = freq_interp(1:floor(N_interp/2)+1);

% Визуализация
figure;
subplot(2,2,1);
plot(f_values, mse_values, 'b-o');
title('MSE в зависимости от частоты f');
xlabel('Частота f (Гц)'); ylabel('MSE');

subplot(2,2,2);
plot(freq, 10*log10(signal_psd), 'b');
title('Спектр исходного сигнала (f = 0.5 Гц)');
xlabel('Частота (Гц)'); ylabel('PSD (дБ)');

subplot(2,2,3);
plot(freq_interp, 10*log10(signal_interp_psd), 'r');
title('Спектр интерполированного сигнала');
xlabel('Частота (Гц)'); ylabel('PSD (дБ)');

subplot(2,2,4);
plot(freq_impact, 10*log10(impact_psd), 'g');
title('Спектр усредненного сигнала');
xlabel('Частота (Гц)'); ylabel('PSD (дБ)');
