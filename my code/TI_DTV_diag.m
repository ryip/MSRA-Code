
% Improvement on original DTV solver by implementing vertical, horizontal,
% and diagonal encouragement norms.
% Original code done by C. Fernandez-Granda
% Modifications done by Richard Yip

% Function which applies transform-invariant directional total-variation 
% regularization to perform image super-resolution. For more information see 
% "Super-resolution via Transform-invariant Group-sparse Regularization" by
% C. Fernandez-Granda and E. Candes
%
% Author: Carlos Fernandez-Granda 
% Email: cfgranda@stanford.edu

function [image_estimate,left_bound, right_bound, top_bound, bottom_bound,x_1,x_2,y_1,y_2] = TI_DTV_diag( img_small,...
    initial_points_small, SRF,sigma,lambda,beta,iter_TFOCS)

[n1_small_aux n2_small_aux]=size(img_small);
n1_big = n1_small_aux * SRF;
n2_big = n2_small_aux * SRF;

% Upsampling using bicubic interpolation
img = imresize(img_small,[n1_big n2_big]);

% Choose an area of the image
initial_points=SRF.*(initial_points_small);
initial_points=floor(initial_points);
focus_size=[initial_points(2, 2)-initial_points(2, 1)+1 initial_points(1, 2)-initial_points(1, 1)+1];
img_center=floor(mean(initial_points, 2));

% Apply TILT to learn the transformation
tic
[Dotau, A, E,f,tfm_matrix, focus_size, error_sign, UData, VData, XData, YData, A_scale]=...
    TILT(uint8(img), 'HOMOGRAPHY', initial_points, 'SAVE_PATH', [], 'DISPLAY_INTER', 0);
toc

% Crop the image
expand_rate=1;
left_bound_aux=ceil(max(initial_points(1, 1)-expand_rate*(initial_points(1, 2)-initial_points(1, 1)), 1));
left_bound = SRF*ceil((left_bound_aux-1)/SRF)+1;
right_bound_aux = floor(min(initial_points(1, 2) + expand_rate*(initial_points(1, 2)-initial_points(1, 1)), size(img, 2)));
right_bound = left_bound + SRF*ceil((right_bound_aux-left_bound)/SRF)-1;
top_bound_aux=ceil(max(initial_points(2, 1)-expand_rate*(initial_points(2, 2)-initial_points(2, 1)), 1));
top_bound = SRF*ceil((top_bound_aux-1)/SRF)+1;
bottom_bound_aux =floor(min(initial_points(2, 2)+expand_rate*(initial_points(2, 2)-initial_points(2, 1)), size(img, 1)));
bottom_bound = top_bound + SRF*ceil((bottom_bound_aux-top_bound)/SRF)-1;
% Make sure that it can be downsampled exactly
n2_small = (right_bound-left_bound+1)/SRF;
n1_small =(bottom_bound-top_bound+1)/SRF;
input_image=img(top_bound:bottom_bound, left_bound:right_bound);
[n_1,n_2]=size(input_image);
center=img_center+[1-left_bound; 1-top_bound];
image_size=size(input_image);
image_center=floor(center);
focus_center=zeros(2, 1);
focus_center(1)=floor((1+focus_size(2))/2);
focus_center(2)=floor((1+focus_size(1))/2);
UData=[1-image_center(1) image_size(2)-image_center(1)];
VData=[1-image_center(2) image_size(1)-image_center(2)];
XData=[1-focus_center(1) focus_size(2)-focus_center(1)];
YData=[1-focus_center(2) focus_size(1)-focus_center(2)];
tfm=fliptform(maketform('projective', tfm_matrix'));

% Create sparse interpolation matrix
tic
interpolation_matrix_aux = create_transformation_matrix_fast(tfm,UData ,VData, XData, YData, n_1,n_2,focus_size);
toc
rect_im=reshape(interpolation_matrix_aux*input_image(:),focus_size(1),focus_size(2));
aux_ones = ones(size(rect_im));
aux_transf = interpolation_matrix_aux'*aux_ones(:);
transf_indices = find(aux_transf>0);
[n1,n2]=size(input_image);
[m1,m2]=size(rect_im);
interpolation_matrix = interpolation_matrix_aux(:,transf_indices);
L = length(transf_indices);

%THIS SECTION
% Create auxiliary operators
mode = 'r2r';
aux_op_f=@(x)aux_op_func(x,n_1,n_2,transf_indices);
aux_op_t=@(x)x(transf_indices);
aux_op = linop_handles( [n_1*n_2,L],aux_op_f, aux_op_t, mode);
TV_x_op=linop_TV_x( [m1,m2] );
TV_y_op=linop_TV_y( [m1,m2] );
TV_DR_op=linop_TV_DR([m1,m2]);
TV_DL_op=linop_TV_DL([m1,m2]);
op_x = @(x)TV_x_op(interpolation_matrix*x,1);
op_y = @(x)TV_y_op(interpolation_matrix*x,1);    
adj_op_x = @(x)interpolation_matrix'*TV_x_op(x,2);
adj_op_y = @(x)interpolation_matrix'*TV_y_op(x,2);
%ADDED following 4 lines of code
op_DR = @(x)TV_DR_op(interpolation_matrix*x,1);
op_DL = @(x)TV_DL_op(interpolation_matrix*x,1);
adj_op_DR = @(x)interpolation_matrix'*TV_DR_op(x,2);
adj_op_DL = @(x)interpolation_matrix'*TV_DL_op(x,2);

%this thing
trans_op_x = linop_handles( [m1*m2,L],op_x, adj_op_x, mode);
trans_op_y = linop_handles( [m1*m2,L],op_y, adj_op_y, mode);
trans_op_DR = linop_handles([m1*m2,L],op_DR,adj_op_DR,mode); %ADDED
trans_op_DL = linop_handles([m1*m2,L],op_DL,adj_op_DL,mode); %ADDED

% Parameters for downsampling filter
% kernel_length = 4*sigma+1;
half = 2*sigma;
[ind_x_ker ind_y_ker] = meshgrid(-half:half);
kernel= exp( -(((ind_x_ker.^2)+(ind_y_ker.^2)) ./ (2* sigma^2)) ); % formula for 2D gaussian
kernel= kernel./sum(kernel(:));
mode = 'r2r';
indices_down_x = ((1:n1_small)*SRF -floor(SRF/2)+1);
indices_down_y = ((1:n2_small)*SRF -floor(SRF/2)+1);

% Downsampling operator
Af  = @(x)A_downsamp_transf(x,kernel,indices_down_x,indices_down_y,transf_indices,n1_small,n2_small,n_1,n_2,0);
At  =  @(x)A_downsamp_transf(x,kernel,indices_down_x,indices_down_y,transf_indices,n1_small,n2_small,n_1,n_2,1);
A = linop_handles( [n1_small*n2_small,L], Af, At, mode);
aux_A = Af(ones(L,1)); 
transf_ind_small = aux_A>0;

% Mask to avoid border effects
half_m = 5;
[ind_m_x ind_m_y] = meshgrid(-half_m:half_m);
sigma_m = 1;
kernel_m= exp( -(((ind_m_x.^2)+(ind_m_y.^2)) ./ (2* sigma_m^2)) ); % formula for 2D gaussian
kernel_m= kernel_m./max(kernel_m(:));
outer_mask = zeros(n1_small,n2_small);
outer_mask(transf_ind_small)=1;
aux_mask = outer_mask;
for ind=1:half_m
    aux_mask(ind:end,half:end)=aux_mask(ind:end,half:end)+ outer_mask(1:(end-ind+1),1:(end-half+1));
    aux_mask(ind:end,1:(end-half+1))=aux_mask(ind:end,1:(end-half+1))+ outer_mask(1:(end-ind+1),half:end);
    aux_mask(1:(end-ind+1),half:end)=aux_mask(1:(end-ind+1),half:end)+ outer_mask(ind:end,1:(end-half+1));
    aux_mask(1:(end-ind+1),1:(end-half+1))=aux_mask(1:(end-ind+1),1:(end-half+1))+ outer_mask(ind:end,half:end);
    aux_mask(half:end,ind:end)=aux_mask(half:end,ind:end)+ outer_mask(1:(end-half+1),1:(end-ind+1));
    aux_mask(1:(end-half+1),ind:end)=aux_mask(1:(end-half+1),ind:end)+ outer_mask(half:end,1:(end-ind+1));
    aux_mask(half:end,1:(end-ind+1))=aux_mask(half:end,1:(end-ind+1))+ outer_mask(1:(end-half+1),ind:end);
    aux_mask(1:(end-half+1),1:(end-ind+1))=aux_mask(1:(end-half+1),1:(end-ind+1))+ outer_mask(half:end,ind:end);
end
inner_mask_small = double(aux_mask==max(aux_mask(:)));
transf_ind_small = find(outer_mask>0);
inner_mask_wind = conv2(inner_mask_small,kernel_m,'same');
inner_mask_wind = inner_mask_wind./max(inner_mask_wind(:));
[aux_ind_x,aux_ind_y]=find(inner_mask_wind>0);
x_1 = SRF*min(aux_ind_x);
x_2 = SRF*max(aux_ind_x);
y_1 = SRF*min(aux_ind_y);
y_2 = SRF*max(aux_ind_y);


% Masked downsampling operator
Af_wind  = @(x)A_downsamp_transf_mask(x,kernel,inner_mask_wind,indices_down_x,indices_down_y,transf_indices,n1_small,n2_small,n_1,n_2,0);
At_wind  = @(x)A_downsamp_transf_mask(x,kernel,inner_mask_wind,indices_down_x,indices_down_y,transf_indices,n1_small,n2_small,n_1,n_2,1);
A_wind = linop_handles( [n1_small*n2_small,L], Af_wind, At_wind, mode);

offset_left = (left_bound-1)/SRF;
offset_top = (top_bound-1)/SRF;
y_aux = img_small( (offset_top+1):(offset_top+n1_small), (offset_left+1):1:(offset_left+n2_small));
y_wind_aux = y_aux.*inner_mask_wind;
y_wind = zeros(size(y_aux));
y_wind(transf_ind_small)=y_wind_aux(transf_ind_small);
y_wind=y_wind(:);

% Setting TFOCS parameters
opts = [];
opts.maxIts     = iter_TFOCS;
opts.printEvery = 50;
opts.tol        = 1e-8;
z0  = [];   % we don't have a good guess for the dual
x0 = ones(L,1);
mu=5e-5;
opts.normA2 = linop_normest( A ).^2;
opts.nonneg = true;
   

if lambda ==0
    [estimate_v, out, optsOut ]=solver_transf_DTV_L2( A_wind, y_wind,trans_op_x,trans_op_y,lambda, mu, x0, z0,L,n_1,n_2,m1,m2, opts );
    image_estimate = zeros(n1,n2);
    image_estimate(transf_indices) = estimate_v;
else
    [estimate_v, out, optsOut ]=solver_transf_DTV_L2_TV_diag( A_wind, y_wind,trans_op_x,trans_op_y,trans_op_DR,trans_op_DL,aux_op,lambda, beta, mu, x0, z0,L,n_1,n_2,m1,m2, opts );
    image_estimate = zeros(n1,n2);
    image_estimate(transf_indices) = estimate_v;
end
