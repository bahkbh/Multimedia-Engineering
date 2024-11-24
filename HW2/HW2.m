clc
% 1. Capture your own image and import the data into computer.
img = imread('my_image.jpeg');
figure;
imshow(img);
title('Original Image');

% (a) Shows its horizontal and vertical sizes, file size, bit depth.
[height, width, channels] = size(img);
fileInfo = dir('my_image.jpeg');
fileSize = fileInfo.bytes;
bitDepth = imfinfo('my_image.jpeg');
fprintf('Horizontal size: %d\n', width);
fprintf('Vertical size: %d\n', height);
fprintf('File size: %d bytes\n', fileSize);
fprintf('Bit depth: %d bits per pixel\n', bitDepth.BitDepth);
% (b) Show a color image.
figure;
imshow(img);
title('Color Image (All RGB Channels)');
% (c) Show its R-channel image only.
figure;
imshow(img(:,:,1));  % Red channel
title('R-Channel Image');
% (d) Show its G-channel image only.
figure;
imshow(img(:,:,2));  % Green channel
title('G-Channel Image');
% (e) Show its B-channel image only.
figure;
imshow(img(:,:,3));  % Blue channel
title('B-Channel Image');


% 2. Process Gray Scale Image
% (a) Show its gray scale image and its histogram
gray_img = rgb2gray(img);
figure;
imshow(gray_img);
title('Grayscale Image');
figure;
histogram(gray_img, 256);
title('Histogram of Grayscale Image');
xlabel('Pixel Intensity');
ylabel('Frequency');
% (b) Show each image with lower bit depth along along with its histogram.
% Original Bit depth
B = 8;
% Reduce to 7 bits
n = 7;
gray_7bit = floor(double(gray_img) / 2^(B-n)) * 2^(B-n);
figure;
imshow(uint8(gray_7bit));
title('7-bit Grayscale Image');
figure;
histogram(uint8(gray_7bit), 128);  % 128 bins for 7-bit values
title('Histogram of 7-bit Grayscale Image');
xlabel('Pixel Intensity');
ylabel('Frequency');
% Reduce to 4 bits
n = 4;
gray_4bit = floor(double(gray_img) / 2^(B-n)) * 2^(B-n);
figure;
imshow(uint8(gray_4bit));
title('4-bit Grayscale Image');
figure;
histogram(uint8(gray_4bit), 16);  % 16 bins for 4-bit values
title('Histogram of 4-bit Grayscale Image');
xlabel('Pixel Intensity');
ylabel('Frequency');
% Reduce to 1 bit
n = 1;
gray_1bit = floor(double(gray_img) / 2^(B-n)) * 2^(B-n);
figure;
imshow(uint8(gray_1bit));
title('1-bit Grayscale Image');
figure;
histogram(uint8(gray_1bit), 2);  % 2 bins for 1-bit values
title('Histogram of 1-bit Grayscale Image');
xlabel('Pixel Intensity');
ylabel('Frequency');


% 3. Represent the RGB Image into indexed color image
% (a) Show the histogram, mean and standard deviation of each color channel.
R_channel = img(:,:,1);  
G_channel = img(:,:,2); 
B_channel = img(:,:,3);  
figure;
subplot(3,1,1), histogram(R_channel, 256), title('R-Channel Histogram');
subplot(3,1,2), histogram(G_channel, 256), title('G-Channel Histogram');
subplot(3,1,3), histogram(B_channel, 256), title('B-Channel Histogram');
mean_R = mean(R_channel(:));
std_R = std(double(R_channel(:)));
mean_G = mean(G_channel(:));
std_G = std(double(G_channel(:)));
mean_B = mean(B_channel(:));
std_B = std(double(B_channel(:)));
fprintf('R Channel Mean: %.2f, Standard Deivation: %.2f\n', mean_R, std_R);
fprintf('G Channel Mean: %.2f, Standard Deivation: %.2f\n', mean_G, std_G);
fprintf('B Channel Mean: %.2f, Standard Deivation: %.2f\n', mean_B, std_B);
% (b) Show the bit-depth reduced color channel image of R, G, and B channel
% image along with their histograms
% Original Bit depth
B = 8;
R_2bit = floor(double(R_channel) / 2^(B-2)) * 2^(B-2);
B_2bit = floor(double(B_channel) / 2^(B-2)) * 2^(B-2);
G_4bit = floor(double(G_channel) / 2^(B-4)) * 2^(B-4);
figure;
imshow(uint8(R_2bit));
title('R-Channel 2-bit Image');
figure;
imshow(uint8(G_4bit));
title('G-Channel 4-bit Image');
figure;
imshow(uint8(B_2bit));
title('B-Channel 2-bit Image');
figure;
subplot(3,1,1), histogram(uint8(R_2bit), 4), title('R-Channel 2-bit Histogram');
subplot(3,1,2), histogram(uint8(G_4bit), 16), title('G-Channel 4-bit Histogram');
subplot(3,1,3), histogram(uint8(B_2bit), 4), title('B-Channel 2-bit Histogram');
mean_R2 = mean(R_2bit(:));
std_R2 = std(double(R_2bit(:)));
mean_G4 = mean(G_4bit(:));
std_G4 = std(double(G_4bit(:)));
mean_B2 = mean(B_2bit(:));
std_B2 = std(double(B_2bit(:)));
fprintf('2 Bits R Channel Mean: %.2f, Standard Deivation: %.2f\n', mean_R2, std_R2);
fprintf('4 Bits G Channel Mean: %.2f, Standard Deivation: %.2f\n', mean_G4, std_G4);
fprintf('2 Bits B Channel Mean: %.2f, Standard Deivation: %.2f\n', mean_B2, std_B2);
% (c) Concatenate the three color image.
reduced_img = cat(3, uint8(R_2bit), uint8(G_4bit), uint8(B_2bit));
figure;
imshow(reduced_img);
title('RGB Image with Reduced Bit Depth');
imwrite(reduced_img, 'reduced_image.png');
fileInfo = dir('reduced_image.png');
fprintf('File size of reduced bit depth image: %d bytes\n', fileInfo.bytes);