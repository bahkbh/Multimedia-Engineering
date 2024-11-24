clc
clear

%%%%%
%%%%%
%%%%%
%%%%%
%%%%%

% < Experiment with Entropy and Huffman coding >
% Understanding Entropy
% 1. Generate three grayscale images of size 100x100 having Gaussian distributed pixel values for the conditions (1),(2),(3)
% Parameters for each condition
params = [
    128, 20; % (1) mean = 128, sigma = 20
    64, 20;  % (2) mean = 64, sigma = 20
    128, 30; % (3) mean = 128, sigma = 30
];
entropies = zeros(1, 3);
for i = 1:3
    mean_val = params(i, 1);   
    sigma_val = params(i, 2);    
    image = normrnd(mean_val, sigma_val, [100, 100]);
    image = max(0, min(255, round(image)));
    
    % (a) Show each image in condition (1),(2),(3).
    figure;
    subplot(1, 2, 1);
    imshow(uint8(image));
    title(['Condition (' num2str(i) ') Image']);
    
    % (b) Show histogram of each image in condition (1),(2),(3).
    subplot(1, 2, 2);
    histogram(image(:), 0:255);
    title(['Condition (' num2str(i) ') Histogram']);
    xlabel('Pixel Intensity');
    ylabel('Frequency');
    
    % (c) Calculate entropy of each image in condition (1),(2),(3).
    entropies(i) = entropy(uint8(image));
end
% (c) Show entropy values for each condition
fprintf('\n');
fprintf('1. (c) Entropy values for each condition:\n');
for j = 1:3
    fprintf("Condition %d Entropy: %f \n",j,entropies(j))
end


%%%%%
%%%%%
%%%%%
%%%%%
%%%%%

% 2. Import the image files into computer and display at your computer using MATLAB.

high_entropy_image_file = 'highent.jpeg';
low_entropy_image_file = 'lowent.jpeg'; 

% (a) Show the two images as color, R-channel only, G-channel only, B-channel only, and grayscale image.
high_entropy_image = imread(high_entropy_image_file);
low_entropy_image = imread(low_entropy_image_file);
figure;
for i = 1:2
    if i == 1
        image = high_entropy_image;
        title_prefix = 'High Entropy Image';
    else
        image = low_entropy_image;
        title_prefix = 'Low Entropy Image';
    end
    gray_image = rgb2gray(image);
    R_channel = image(:,:,1);
    G_channel = image(:,:,2);
    B_channel = image(:,:,3);

    subplot(5, 2, i);
    imshow(image);
    title([title_prefix, ' - Color']);

    subplot(5, 2, i + 2);
    imshow(R_channel);
    title([title_prefix, ' - R Channel']);

    subplot(5, 2, i + 4);
    imshow(G_channel);
    title([title_prefix, ' - G Channel']);

    subplot(5, 2, i + 6);
    imshow(B_channel);
    title([title_prefix, ' - B Channel']);

    subplot(5, 2, i + 8);
    imshow(gray_image);
    title([title_prefix, ' - Grayscale']);
end

% (b) Investigate the RGB channels to find out what channel (among RGB) is most similar to its gray scale image and explain why.
fprintf('\n')
fprintf('2. (b) Investigating RGB channels similarity to grayscale:');
for i = 1:2
    if i == 1
        image = high_entropy_image;
        image_name = 'High Entropy Image';
    else
        image = low_entropy_image;
        image_name = 'Low Entropy Image';
    end
    gray_image = rgb2gray(image);

    mse_R = mean((double(image(:,:,1)) - double(gray_image)).^2, 'all');
    mse_G = mean((double(image(:,:,2)) - double(gray_image)).^2, 'all');
    mse_B = mean((double(image(:,:,3)) - double(gray_image)).^2, 'all');
    fprintf('\n')
    fprintf(['Similarity for ', image_name, ':\n']);
    fprintf(['  R Channel: ', num2str(mse_R)]);
    fprintf('\n')
    fprintf(['  G Channel: ', num2str(mse_G)]);
    fprintf('\n')
    fprintf(['  B Channel: ', num2str(mse_B)]);
    [~, most_similar_channel] = min([mse_R, mse_G, mse_B]);
    channel_names = {'R', 'G', 'B'};
    fprintf('\n')
    fprintf(['Most similar channel: ', channel_names{most_similar_channel}]);
    fprintf('\n')
end

% (c) Show histograms of the two images (only for grayscale image).
figure;
for i = 1:2
    if i == 1
        image = high_entropy_image;
        title_prefix = 'High Entropy Image';
    else
        image = low_entropy_image;
        title_prefix = 'Low Entropy Image';
    end

    gray_image = rgb2gray(image);

    subplot(2, 1, i);
    histogram(gray_image, 0:255);
    title([title_prefix, ' - Grayscale Histogram']);
    xlabel('Pixel Intensity');
    ylabel('Frequency');
end

% (d) Calculate entropies of the two images (only for grayscale image).
fprintf('\n')
fprintf('\n')
fprintf('2. (d) Entropy values for grayscale images:');
high_entropy_gray = rgb2gray(high_entropy_image);
low_entropy_gray = rgb2gray(low_entropy_image);
high_entropy_value = entropy(high_entropy_gray);
low_entropy_value = entropy(low_entropy_gray);
fprintf('\n')
fprintf(['  High Entropy Image: ', num2str(high_entropy_value)]);
fprintf('\n')
fprintf(['  Low Entropy Image: ', num2str(low_entropy_value)]);

%%%%%
%%%%%
%%%%%
%%%%%
%%%%%

%3. Perform processing on each channel of the given BerlineTower image (BerlinTower.bmp) using matlab as below.

image = imread('BerlinTower.bmp');

% (a), (b), (c) Show its R-channel data, G-channel data, B-channel data, that is, Image, histogram, and entropy.
channels = {'R', 'G', 'B'};
figure;
for i = 1:3
    channel_image = image(:,:,i);
    subplot(3, 3, (i-1)*3 + 1);
    imshow(channel_image);
    title([channels{i}, '-Channel Image']);
    subplot(3, 3, (i-1)*3 + 2);
    histogram(channel_image, 0:255);
    title([channels{i}, '-Channel Histogram']);
    xlabel('Pixel Intensity');
    ylabel('Frequency');
    channel_entropy = entropy(channel_image);
    subplot(3, 3, (i-1)*3 + 3);
    text(0.5, 0.5, ['Entropy: ', num2str(channel_entropy)], 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    axis off;
    title([channels{i}, '-Channel Entropy']);
end

ycbcr_image = rgb2ycbcr(image);

% (d), (e), (f) Show its Y-channel data, Cb-channel data, Cr-channel data, that is, Image, histogram, and entropy.
ycbcr_channels = {'Y', 'Cb', 'Cr'};
figure;
for i = 1:3
    channel_image = ycbcr_image(:,:,i);
    subplot(3, 3, (i-1)*3 + 1);
    imshow(channel_image);
    title([ycbcr_channels{i}, '-Channel Image']);
    subplot(3, 3, (i-1)*3 + 2);
    histogram(channel_image, 0:255);
    title([ycbcr_channels{i}, '-Channel Histogram']);
    xlabel('Pixel Intensity');
    ylabel('Frequency');
    channel_entropy = entropy(channel_image);
    subplot(3, 3, (i-1)*3 + 3);
    text(0.5, 0.5, ['Entropy: ', num2str(channel_entropy)], 'FontSize', 12, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    axis off;
    title([ycbcr_channels{i}, '-Channel Entropy']);
end

%%%%%
%%%%%
%%%%%
%%%%%
%%%%%

% Understanding Huffman coding
% 4. Perform Huffman coding of the given KoreanSoccer image (KoreanSoccer.bmp) using matlab as below.
% MATLAB Code for Huffman and Run-Length Coding on KoreanSoccer.bmp
% 하위 폴더를 MATLAB 경로에 추가
addpath('RLE___IRLE_coding');

% Load the image and convert to grayscale
image_rgb = imread('KoreanSoccer.bmp');
gray_image = rgb2gray(image_rgb);

% (a) Convert the KoreanSoccer image to a gray image. Show its gray scale image and histogram.
figure;
subplot(1, 2, 1);
imshow(gray_image);
title('Grayscale Image');
subplot(1, 2, 2);
histogram(gray_image, 0:255);
title('Histogram of Grayscale Image');
xlabel('Pixel Intensity');
ylabel('Frequency');

% (b) Reduce the bit depth of the gray image obtained in 4(a) to 4-bit depth. Show its image and histogram.
B = 8;
n = 4;
gray_4bit_image = floor(double(gray_image) / (2^(B-n))) * 2^(B-n);
% Display 4-bit grayscale image and histogram
figure;
subplot(1, 2, 1);
imshow(uint8(gray_4bit_image));
title('4-Bit Depth Grayscale Image');
subplot(1, 2, 2);
histogram(gray_4bit_image, 0:2^(n)-1);
title('Histogram of 4-Bit Depth Image');
xlabel('Pixel Intensity');
ylabel('Frequency');

% (c) Calculate the average code length if Huffman coding is applied to a 4-bit depth gray image.
data2encode = gray_4bit_image(:);
symbolList = unique(data2encode);
counts = histc(data2encode, symbolList);
probList = counts / sum(counts);
[dict, avgCodeLength] = huffmandict(symbolList, probList);
fprintf('\n\n');
fprintf('4. (c) Average code length with Huffman coding: %f\n', avgCodeLength);
fprintf('\n');

% (d) Calculate the compression ratio if Huffman coding is applied to a 4-bit depth gray image.
original_bits = numel(data2encode) * 4;

huff_encoded = huffmanenco(data2encode, dict);
if iscell(huff_encoded)
    huffman_bits = sum(cellfun('length', huff_encoded));
else
    huffman_bits = length(huff_encoded);
end

compression_ratio_huff = original_bits / huffman_bits;
fprintf('\n');
fprintf('4. (d) Compression ratio with Huffman coding: %f\n', compression_ratio_huff);
fprintf('\n');

% (e) Encode 4-bit depth gray scale image obtained in 4(b) with run-length coding.
data_run_length = gray_4bit_image(:);
Output = rle(data_run_length);
values = Output(1:2:end);       
run_lengths = Output(2:2:end);  

% (f) Calculate the average code length if Huffman coding is applied to run-length coded gray image.
symbolList_values = unique(values);
counts_values = histc(values, symbolList_values);
probList_values = counts_values / sum(counts_values);
[dict_values, avgCodeLength_values] = huffmandict(symbolList_values, probList_values);
symbolList_runs = unique(run_lengths);
counts_runs = histc(run_lengths, symbolList_runs);
probList_runs = counts_runs / sum(counts_runs);
[dict_runs, avgCodeLength_runs] = huffmandict(symbolList_runs, probList_runs);
total_symbols = numel(values) + numel(run_lengths);
avgCodeLength_total = (avgCodeLength_values * numel(values) + avgCodeLength_runs * numel(run_lengths)) / total_symbols;
fprintf('\n');
fprintf('4. (f) Average code length with Huffman coding on run-length coded image: %f\n', avgCodeLength_total);
fprintf('\n');

% (g) Calculate the compression ratio if Huffman coding is applied to run-length coded gray image.
max_value = max(values);
max_run_length = max(run_lengths);
bits_per_value = ceil(log2(double(max_value) + 1));
bits_per_run_length = ceil(log2(double(max_run_length) + 1));
original_bits_run = numel(values) * bits_per_value + numel(run_lengths) * bits_per_run_length;
code_lengths_values = cellfun('length', dict_values(:,2));
huffman_bits_values = sum(counts_values .* code_lengths_values');
code_lengths_runs = cellfun('length', dict_runs(:,2));
huffman_bits_runs = sum(counts_runs .* code_lengths_runs');
huffman_bits_run = huffman_bits_values + huffman_bits_runs;
compression_ratio_huff_run = original_bits_run / huffman_bits_run;
fprintf('4. (g) Compression ratio with Huffman coding on run-length coded image: %f\n', compression_ratio_huff_run);
