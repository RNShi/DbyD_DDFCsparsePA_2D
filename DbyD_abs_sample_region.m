% 2D Dimension-by-Dimension domain decomposition sparse PA method of chosen 
% sample region (add the same number of extra points to each boundaries or 
% the number of extra points depend on the average total variation)
% F1: the input image data (absolutely value)
% (a1:a2,b1:b2): the sample region
% res: position integer, change resolution of reconstruction
% n_more/m_more: set number of extra points added into the boundaries of each subdomain
% F3: the reconstruction of sample region(output)

close all;
clear all; 
format long

% p = parpool;     % open parallel pool using the default profile to define the number of workers
%% Input image data, F1(absolute)------------------------
sample = input('Please chosse the sample region (1-3): '); 
kind_of_extra = input('Please chosse the extra points = (1-length of sample region/10; 2-length of sample region/4;\n 3-length of sample region/2; 4-find by average total variation): '); 

% test image(shepp-logan image)
fn_r = importdata('shepp801_85_r.txt');
fn_i = importdata('shepp801_85_i.txt');
F1 = fn_r+fn_i*1i;
F1 = abs(F1);


[n1,m1]=size(F1); % size of the input data
res = 2;  % change resolution of reconstruction
maxff = max(max(F1));

% choose sample region:
if sample == 1
    a1 = 30; a2 = 50;
    b1 = 29; b2 = 49;
elseif sample == 2
    a1 = 30; a2 = 50;
    b1 = 53; b2 = 73;
elseif sample == 3
    a1 = 50; a2 = 70;
    b1 = 15; b2 = 35;
end

% size of sample region
n = a2-a1+1; m = b2-b1+1;

% set # of extra points to the boundary of each subdomain: 
if kind_of_extra == 1
    n_more = floor(n/10);
    m_more = floor(m/10);
    lambda = 0.0175;
elseif kind_of_extra == 2
    n_more = floor(n/4);
    m_more = floor(m/4);
    lambda = 0.015;
elseif kind_of_extra == 3
    n_more = floor(n/2);
    m_more = floor(m/2);
    lambda = 0.0125;
elseif kind_of_extra == 4
    lambda = 0.02; 
end

% the extension of the sample region: f_new
if kind_of_extra >= 1 && kind_of_extra <= 3
    F2 = zeros(n1+2*n_more,m1+2*m_more);
    F2(n_more+1:n_more+n1, m_more+1:m_more+m1) = F1;
    F2(          1:n_more, m_more+1:m_more+m1) = F1(end-n_more+1:end, :);
    F2(   n_more+n1+1:end, m_more+1:m_more+m1) = F1(1:n_more        , :);
    F2(n_more+1:n_more+n1, 1:m_more          ) = F1(:               , end-m_more+1:end);
    F2(n_more+1:n_more+n1, m_more+m1+1:end   ) = F1(:               , 1:m_more);
    F2(          1:n_more, 1:m_more          ) = F1(end-n_more+1:end, end-m_more+1:end);
    F2(          1:n_more, m_more+m1+1:end   ) = F1(end-n_more+1:end, 1:m_more);
    F2(   n_more+n1+1:end, 1:m_more          ) = F1(1:n_more        , end-m_more+1:end);
    F2(   n_more+n1+1:end, m_more+m1+1:end   ) = F1(1:n_more        , 1:m_more);

    aa1 = a1; aa2 = a2 + 2*n_more;
    bb1 = b1; bb2 = b2 + 2*m_more;
    f_new = F2(aa1:aa2,bb1:bb2);
elseif kind_of_extra == 4
    n_more_min = 2;
    if n/10 <= 10
        beta = min(max(10, floor(n/10)),floor(n/3));  
    else
        beta = max(10, floor(n/30));
    end
    n_more_max = beta;
    n_more_max = floor(n/2);
    n_new0 = n+n_more_max*2;

    F2 = zeros(n1+2*n_more_max,m1+2*n_more_max);
    F2(n_more_max+1:n_more_max+n1, n_more_max+1:n_more_max+m1) = F1;
    F2(          1:n_more_max, n_more_max+1:n_more_max+m1) = F1(end-n_more_max+1:end, :);
    F2(   n_more_max+n1+1:end, n_more_max+1:n_more_max+m1) = F1(1:n_more_max        , :);
    F2(n_more_max+1:n_more_max+n1, 1:n_more_max          ) = F1(:               , end-n_more_max+1:end);
    F2(n_more_max+1:n_more_max+n1, n_more_max+m1+1:end   ) = F1(:               , 1:n_more_max);
    F2(          1:n_more_max, 1:n_more_max          ) = F1(end-n_more_max+1:end, end-n_more_max+1:end);
    F2(          1:n_more_max, n_more_max+m1+1:end   ) = F1(end-n_more_max+1:end, 1:n_more_max);
    F2(   n_more_max+n1+1:end, 1:n_more_max          ) = F1(1:n_more_max        , end-n_more_max+1:end);
    F2(   n_more_max+n1+1:end, n_more_max+m1+1:end   ) = F1(1:n_more_max        , 1:n_more_max);
    
    aa1 = a1; aa2 = a2 + 2*n_more_max;
    bb1 = b1; bb2 = b2 + 2*n_more_max;
    f_new = F2(aa1:aa2,bb1:bb2);
end

% 1D domain decomposition Fourier continuation sparse PA reconstruction in y-direction
N = size(f_new,1); M = size(f_new,2);

if kind_of_extra >= 1 && kind_of_extra <= 3
    parfor ii = 1:N
        f_new1(:,ii) = function_1DDDFCSPA(f_new(ii,:).', res, lambda);
    end 
    f_new1 = f_new1.';
elseif kind_of_extra == 4
    n_more_aTV_y = zeros(n_new0,2);
    n_more_aTV_y = function_aTV_extra_point_y(F2, n_new0, n_more_min, n_more_max, beta, a1, a2, b1, b2);
    f_new1 = zeros(M*res-(res-1),N);
    for ii = 1:N
            f_new1((1+n_more_max-n_more_aTV_y(ii,1))*res-(res-1):...
                   (M-n_more_max+n_more_aTV_y(ii,2))*res-(res-1),ii) = function_1DDDFCSPA(...
                f_new(ii,1+n_more_max-n_more_aTV_y(ii,1):...
                         M-n_more_max+n_more_aTV_y(ii,2)).', res, lambda);
    end 
    f_new1 = f_new1.';
end

% 1D domain decomposition Fourier continuation sparse PA reconstruction in x-direction
N1 = size(f_new1,1); M1 = size(f_new1,2);

if kind_of_extra >= 1 && kind_of_extra <= 3
    parfor ii = 1:M1
        f_new2(:,ii) = function_1DDDFCSPA(f_new1(:,ii), res, lambda);
    end
elseif kind_of_extra == 4
    n_new1 = n_new0*res-(res-1);
    n_more_aTV_x = zeros(n_new1,2);
    n_more_aTV_x = function_aTV_extra_point_x(f_new1, n_new1, n_more_min, n_more_max, beta, n);
    f_new2 = zeros(N*res-(res-1),M1);
    for ii = (1+n_more_max-n_more_min)*res-(res-1):M1-(n_more_max-n_more_min)*res
        f_new2((1+n_more_max-n_more_aTV_x(ii,1))*res-(res-1):...
             (N1-n_more_max+n_more_aTV_x(ii,2))*res-(res-1),ii) = function_1DDDFCSPA(...
             f_new1(1+n_more_max-n_more_aTV_x(ii,1):...
                   N1-n_more_max+n_more_aTV_x(ii,2),ii), res, lambda);
    end  
end

% the final reconstruction
N2 = size(f_new2,1); M2 = size(f_new2,2);
if kind_of_extra >= 1 && kind_of_extra <= 3
    F3 = f_new2(n_more*res+1:n_more*res+(n-1)*res, ...
                m_more*res+1:m_more*res+(m-1)*res);
elseif kind_of_extra == 4
    F3 = f_new2(n_more_max*res+1:end-n_more_max*res, ...
                n_more_max*res+1:end-n_more_max*res);
end

% plot original image
figure, imagesc(F1(a1:a2,b1:b2),[0,maxff]), colormap(gray)
set(gcf, 'Position',  [100, 100, 320, 240])
set(gca,'xtick',[]),set(gca,'ytick',[])
title('original image (Gibbs oscillations)')
colorbar

% plot reconstructed image
figure, imagesc(F3,[0,maxff]),colormap(gray)
set(gcf, 'Position',  [100, 100, 320, 240])
set(gca,'xtick',[]),set(gca,'ytick',[])
title('reconstructed image')
colorbar



