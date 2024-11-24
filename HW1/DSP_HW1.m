% 1. Checking the Audio Data
% (a) Read the given original audio data into Matlab
[audioData, sampleRate] = audioread('PianoSound-44.1kHz.wav');
% (b) Listen to the sound
sound(audioData, sampleRate);
% (c) Plot the waveform with its x-axis in time with unit of seconds
time = (0:length(audioData)-1) / sampleRate;
% 파형 그리기
figure;
plot(time, audioData);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of Audio Data');
xlim("tight")
ylim([-1 1])


% 2. Playing with the sampling frequency of the audio data by resampling
% (a) Resample the original audio data at half of the original sampling frequency.
resampleRate_half = sampleRate / 2;
resampledAudio_half = resample(audioData, resampleRate_half, sampleRate);
% 리샘플된 오디오 데이터를 파일로 저장
audiowrite('ResampledAudio_half.wav', resampledAudio_half, resampleRate_half);
% 파형 그리기
time_half = (0:length(resampledAudio_half)-1) / resampleRate_half;
figure;
plot(time_half, resampledAudio_half);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of Audio Data Resampled at Half the Original Sampling Frequency');
xlim("tight")
ylim([-1 1])

% (b) Resample the original audio data at 1/10 of the original sampling frequency.
resampleRate_10th = sampleRate / 10;
resampledAudio_10th = resample(audioData, resampleRate_10th, sampleRate);
% 리샘플된 오디오 데이터를 파일로 저장
audiowrite('ResampledAudio_10th.wav', resampledAudio_10th, resampleRate_10th);
% 파형 그리기
time_10th = (0:length(resampledAudio_10th)-1) / resampleRate_10th;
figure;
plot(time_10th, resampledAudio_10th);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of Audio Data Resampled at 1/10th the Original Sampling Frequency');
xlim("tight")
ylim([-1 1])

% (c) Resample the original audio data at 1/50 of the original sampling frequency
resampleRate_50th = sampleRate / 50;
resampledAudio_50th = resample(audioData, resampleRate_50th, sampleRate);
% 리샘플된 오디오 데이터를 파일로 저장
audiowrite('ResampledAudio_50th.wav', resampledAudio_50th, resampleRate_50th);
% 파형 그리기
time_50th = (0:length(resampledAudio_50th)-1) / resampleRate_50th;
figure;
plot(time_50th, resampledAudio_50th);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Waveform of Audio Data Resampled at 1/50th the Original Sampling Frequency');
xlim("tight")
ylim([-1 1])


% 3. Analyzing the audio data in the frequency domain
% (a) conversion of original audio data and drawing waveforms in the frequency domain
% FFT 수행
n = length(audioData); % 데이터 길이
fftOriginal = fft(audioData); % FFT 계산
f = (0:n-1)*(sampleRate/n); % 주파수 벡터 생성
magnitudeOriginal = abs(fftOriginal/n); % FFT의 크기 스펙트럼
% 0에서 1000Hz 범위의 주파수 대역만 표시
idx = f <= 1000;
figure;
plot(f(idx), magnitudeOriginal(idx));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain of Original Audio (0 to 1000 Hz)');

% (b) Resampled audio data and drawing waveform in the frequency domain
% 각 리샘플된 오디오 데이터에 대해 FFT 수행
% 절반 샘플링 주파수 데이터
n_half = length(resampledAudio_half);
fft_half = fft(resampledAudio_half);
f_half = (0:n_half-1)*(resampleRate_half/n_half);
magnitude_half = abs(fft_half/n_half);
% 0에서 1000Hz 범위의 주파수 대역만 표시
idx_half = f_half <= 1000;
figure;
plot(f_half(idx_half), magnitude_half(idx_half));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain of Audio Resampled at Half (0 to 1000 Hz)');

% 1/10 샘플링 주파수 데이터
n_10th = length(resampledAudio_10th);
fft_10th = fft(resampledAudio_10th);
f_10th = (0:n_10th-1)*(resampleRate_10th/n_10th);
magnitude_10th = abs(fft_10th/n_10th);
% 0에서 1000Hz 범위의 주파수 대역만 표시
idx_10th = f_10th <= 1000;
figure;
plot(f_10th(idx_10th), magnitude_10th(idx_10th));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain of Audio Resampled at 1/10th (0 to 1000 Hz)');

% 1/50 샘플링 주파수 데이터
n_50th = length(resampledAudio_50th);
fft_50th = fft(resampledAudio_50th);
f_50th = (0:n_50th-1)*(resampleRate_50th/n_50th);
magnitude_50th = abs(fft_50th/n_50th);
% 0에서 1000Hz 범위의 주파수 대역만 표시
idx_50th = f_50th <= 1000;
figure;
plot(f_50th(idx_50th), magnitude_50th(idx_50th));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency Domain of Audio Resampled at 1/50th (0 to 1000 Hz)');
xlim([0 1000]);