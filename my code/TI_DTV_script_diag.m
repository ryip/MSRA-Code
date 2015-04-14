% Improvement on original DTV solver by implementing vertical, horizontal,
% and diagonal encouragement norms.
% Original code done by C. Fernandez-Granda
% Modifications done by Richard Yip

% Script to super-resolve an image by applying transform-invariant 
% directional total-variation regularization. For more information see 
% "Super-resolution via Transform-invariant Group-sparse Regularization" by
% C. Fernandez-Granda and E. Candes
%
% Author: Carlos Fernandez-Granda
% Contact: cfgranda@stanford.edu

% addpath('C:\Users\v-riyip\Documents\MATLAB\TI_DTV_code\TILT_v1_04')
% addpath('C:\Users\v-riyip\Documents\MATLAB\TI_DTV_code\TFOCS')
% addpath('C:\Users\v-riyip\Documents\MATLAB\TI_DTV_code\code')
clear all
close all

% 
SRF =4 ;
sigma = SRF/2+1.5;
% lambda controls the penalization on the TI_DTV term in the cost function
lambda = 4; 
% Increase beta if bright spots appear in certain regions
beta=0.5;

image_number=6; % Change this number to super-resolve other images

img_name = ['C:\Users\Richard\Documents\ MSRA\TI_DTV_code\images\img' num2str(image_number) '.jpg'];
iml = imread(img_name);
iml_ycbcr = rgb2ycbcr(iml);
img_small = double(iml_ycbcr(:, :, 1));
[n1_small_aux, n2_small_aux]=size(img_small);
% The initial points determine the subimage. You can use the script
% select_subimage to save new initial points. 
if image_number == 2
    initial_point_string = ['C:\Users\Richard\Documents\ MSRA\TI_DTV_code\initial_points\initial_points_' num2str(image_number) '_1'];
else
    initial_point_string = ['C:\Users\Richard\Documents\ MSRA\TI_DTV_code\initial_points\initial_points_' num2str(image_number)];
end
load(initial_point_string);
iter_TFOCS=200;
tic
[image_estimate,left_bound, right_bound, top_bound, bottom_bound,x_1,x_2,y_1,y_2] = TI_DTV_diag( img_small,...
        initial_points, SRF,sigma,lambda,beta,iter_TFOCS );
toc
iml_cb = iml_ycbcr(:, :, 2);
iml_cr = iml_ycbcr(:, :, 3);
imh_cb_big = imresize(iml_cb, SRF , 'bicubic');
imh_cr_big = imresize(iml_cr, SRF, 'bicubic');
imh_cb=imh_cb_big(top_bound:bottom_bound, left_bound:right_bound);
imh_cr=imh_cr_big(top_bound:bottom_bound, left_bound:right_bound);
imh_ycbcr = cat(3, uint8(image_estimate), cat(3, imh_cb, imh_cr));
imh = ycbcr2rgb(uint8(imh_ycbcr));
x_1_lowres =(top_bound+x_1)./SRF;
x_2_lowres =(top_bound+x_2)./SRF;
y_1_lowres =(left_bound+y_1)./SRF;
y_2_lowres =(left_bound+y_2)./SRF;
figure;
imshow(iml(x_1_lowres:x_2_lowres,y_1_lowres:y_2_lowres,:),[])
title('Low resolution image')
figure
imshow(imh(x_1:x_2,y_1:y_2,:),[])
title('Super-resolved image')
