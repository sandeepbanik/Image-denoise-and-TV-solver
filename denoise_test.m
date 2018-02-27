%% Denoise example
    clear;
    close all;
    clc;
    addpath('Denoise');
%% Data input
    img_clr = imread('test_img.jpg');
    img_gray = double(rgb2gray(img_clr));
    X = img_gray./255;
    b = X + 2e-1*rand(size(X));

%% Parameter
    lambda = 0.1;
    max_iter = 200;
    [X_opt,h] = denoise_algo(b, lambda, max_iter, -Inf, Inf);
    [FX_opt,g] = Fdenoise_algo(b, lambda, max_iter, -Inf, Inf);
    
%% Plot
    b_scale = 255*b;
    x_scale = X_opt*255;
    Fx_scale = FX_opt*255;
    imshow(uint8(x_scale));
    subplot(1,3,1);
    imshow(uint8(b_scale));
    title('Corrupted image');
    subplot(1,3,2);
    imshow(uint8(x_scale));
    title('FISTA image');
    subplot(1,3,3);
    imshow(uint8(Fx_scale));
    title('MFISTA image');
    
    