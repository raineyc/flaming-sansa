clear
clc

%load('/home/cameron/MATLAB/pothole data/CPGp83fab19110606.mat')

%im = imread('shapes.bmp');
%im = imread('Valve_original_(1).PNG');
%data = im2double(rgb2gray(im));


G_filt = fspecial('gaussian', 31, 5);
%G_filt = fspecial('unsharp');



data_filt = imfilter(data,G_filt,'symmetric');

[edges_bin, edges_angle, edges_angle_im] = find_canny_edges(data_filt);
figure; imshow(edges_angle_im)

