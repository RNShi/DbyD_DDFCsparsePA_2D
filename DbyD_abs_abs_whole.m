% 2D Dimension-by-Dimension domain decomposition sparse PA method of whole
% image
% F1: the input image data(take absolute value if it is complex)
% part_x: the number of subdomains in x direction
% part_y: the number of subdomains in y direction
% res: position integer, change resolution of reconstruction
% n_more/m_more: set number of extra points added into the boundaries of each subdomain
% F3: the reconstructed image



close all;
clear all; 
format long

% p = parpool;     % open parallel pool using the default profile to define the number of workers
%% Input image data, F1(absolute)------------------------
% the # of subdomain in x and y-directions
% part_x = input('The number of subdomains in x direction: ');
% part_y = input('The number of subdomains in y direction: ');
% kind_of_extra = input('Please chosse the extra points = (1-length of sample region/10; 2-length of sample region/4;\n 3-length of sample region/2): '); 

part_x = 6;
part_y = 7;
kind_of_extra = 1;

% test image(shepp-logan image)
fn_r = importdata('shepp801_85_r.txt');
fn_i = importdata('shepp801_85_i.txt');
F1 = fn_r+fn_i*1i;
F1 = abs(F1);

[n1,m1]=size(F1); % size of the input data
res = 2;  % change resolution of reconstruction
maxff = max(max(F1));

%% domain decomposition Fourier continuation sparse PA method (total part_x*part_y subdomain)

% size of each subdomain
n = floor(n1/part_x)+1;
m = floor(m1/part_y)+1;

F3 = zeros(n1*2-1,m1*2-1); % output data(reconstruction)
for py = 6:7
for px = 5:6

% find endpoints for each subdomain
if px ~= part_x
    a1 = (px-1)*(n-1)+1; a2 = px*(n-1)+1;
else
    a1 = (px-1)*(n-1)+1; a2 = n1;
end       
if py ~= part_y
    b1 = (py-1)*(m-1)+1; b2 = py*(m-1)+1;
else
    b1 = (py-1)*(m-1)+1; b2 = m1;
end

% set # of extra points to the boundary of each subdomain: 
if kind_of_extra == 1
    n_more = ceil(n/10);
    m_more = ceil(m/10);
    lambda = 0.0175;
elseif kind_of_extra == 2
    n_more = ceil(n/4);
    m_more = ceil(m/4);
    lambda = 0.015;
elseif kind_of_extra == 3
    n_more = ceil(n/2);
    m_more = ceil(m/2);
    lambda = 0.0125;
end

% the extension of the sample region: f_new
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


N = size(f_new,1); 
M = size(f_new,2);

% 1D domain decomposition Fourier continuation sparse PA reconstruction in y-direction
N = size(f_new,1); M = size(f_new,2);
parfor ii = 1:N
    f_new1(:,ii) = function_1DDDFCSPA(f_new(ii,:).', res, lambda);
end 
f_new1 = f_new1.';

% 1D domain decomposition Fourier continuation sparse PA reconstruction in x-direction
N1 = size(f_new1,1); M1 = size(f_new1,2);
parfor ii = 1:M1
    f_new2(:,ii) = function_1DDDFCSPA(f_new1(:,ii), res, lambda);
end

% the final reconstruction
N2 = size(f_new2,1); M2 = size(f_new2,2);
if px ~= part_x && py ~= part_y
    F3((px-1)*(n-1)*res+1:px*(n-1)*res, (py-1)*(m-1)*res+1:py*(m-1)*res)...
        = f_new2(n_more*res+1:n_more*res+(n-1)*res,...
                 m_more*res+1:m_more*res+(m-1)*res);
elseif px == part_x && py ~= part_y
    F3((px-1)*(n-1)*res+1:end         , (py-1)*(m-1)*res+1:py*(m-1)*res)...
        = f_new2(n_more*res+1:n_more*res+(n1-((px-1)*(n-1)+1)+1)*res-1,...
                 m_more*res+1:m_more*res+(m-1)*res);
elseif px ~= part_x && py == part_y
    F3((px-1)*(n-1)*res+1:px*(n-1)*res, (py-1)*(m-1)*res+1:end         )...
        = f_new2(n_more*res+1:n_more*res+(n-1)*res,...
                 m_more*res+1:m_more*res+(m1-((py-1)*(m-1)+1)+1)*res-1);   
elseif px == part_x && py == part_y
    F3((px-1)*(n-1)*res+1:end         , (py-1)*(m-1)*res+1:end         )...
        = f_new2(n_more*res+1:n_more*res+(n1-((px-1)*(n-1)+1)+1)*res-1,...
                 m_more*res+1:m_more*res+(m1-((py-1)*(m-1)+1)+1)*res-1);
end

px,py
clear f_new1 f_new2
end
end

% plot original image
figure, imagesc(F1,[0,maxff]), colormap(gray)
set(gca,'xtick',[]),set(gca,'ytick',[])
title('original image')
colorbar

% plot reconstructed image
figure, imagesc(F3,[0,maxff]), colormap(gray)
set(gca,'xtick',[]),set(gca,'ytick',[])
title('reconstructed image')
colorbar



