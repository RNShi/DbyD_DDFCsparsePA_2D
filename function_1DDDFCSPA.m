% function: 1D domain decomposition Fourier continuation sparse PA
% reconstruction

function frecon = function_1DDDFCSPA(f_new, res, lambda)
N = size(f_new,1); % size of input data
x0 = 0; x1 = 1;
N1 = (N-1)*res+1; % size of output data

% ----- some parameters ------     
% beta = the size of the data near the boundaries used for the FC 
if N/10 <= 10
    beta = min(max(10, floor(N/10)),floor(N/3));  
else
    beta = max(10, floor(N/30));
end
% gamma = the size of the extension data for the FC 
gamma = 2*beta+1; 
% m=M<beta
if beta>5
    m = 5;
else
    m = beta-1;
end  
Q     = 100;     
mm    = 30 ;     % mm = g
dx    = (x1-x0)/(N-1);  % dx=dx
delta = beta*dx;     % delta=beta*dx
d     = (gamma-1)*dx;   % d=gamma*dx
d_Q   = delta/Q;    % d_Q=delta/Q
x = zeros(N,1);
for ii = 1:N
    x(ii) = x0 + (ii-1)*dx;    
end
sN = N+gamma-1; % size of the periodic function reconstructed by FC
bx = (x1-x0)/(N-1)*(sN-1); 
sN1 = (sN-1)*res+1;

%------------- The points near the boundary -------------------
x_left (1:beta+1) = x(N-beta:N);
x_right(1:beta+1) = x(1:beta+1);

%------------ recurrence Gram polynomial ----------------------
a0 = x1 - delta;  a1 = x1;
x_l = 2*(x_left-a0)/(a1-a0) - 1;    %x_l in [-1,1] beta+1 points
a0 = x0;  a1 = x0+delta;
x_r = 2*(x_right-a0)/(a1-a0) - 1;   %x_r in [-1,1] beta+1 points

alpha = zeros(m,1);
for n = 1: m
    alpha(n) = beta/n*sqrt((4*n^2-1)/((beta+1)^2-n^2));
end
P_l = zeros(beta+1, m+1);   %(beta+1,M+1)    
P_r = zeros(beta+1, m+1);  
P_l(:,1) = 1/sqrt(beta+1);  P_l(:,2) = alpha(1)*x_l'.*P_l(:,1);
P_r(:,1) = 1/sqrt(beta+1);  P_r(:,2) = alpha(1)*x_r'.*P_r(:,1);
for ii = 3: m+1
    P_l(:,ii) = alpha(ii-1)*x_l'.*P_l(:,ii-1) - alpha(ii-1)/alpha(ii-2)*P_l(:,ii-2);
    P_r(:,ii) = alpha(ii-1)*x_r'.*P_r(:,ii-1) - alpha(ii-1)/alpha(ii-2)*P_r(:,ii-2);
end

%------ function values in [x0,x0+delta] and [x1+d,x1+d+delta] ------------
f_left  = f_new(N-beta:N,1);
f_right = f_new(1:beta+1,1);
a_left = zeros(m+1,1);
a_right = zeros(m+1,1);
for ii = 1:m+1
    a_left(ii,1)  = sum(f_left(:,1) .*(P_l(:,ii))); %a_left=b_l
    a_right(ii,1) = sum(f_right(:,1).*(P_r(:,ii)));
end

% ---------- Compute f_even and f_odd ---------------------------
x_d = zeros(2*(beta+gamma-1),1);
for ii = 1:2*(beta+gamma-1)
    x_d(ii) = (x1-delta) + (ii-1)*dx;   % x in [b-delta,b+2*d+delta]  2*(beta+gamma) points
end
x_d_left = zeros(Q+1,1);
for ii = 1:Q+1
    x_d_left(ii) = (x1-delta) + (ii-1)*d_Q;    % x in [b-delta,b] Q+1  points
end

% ------------ recurrence Gram polynomial -------------------------
a0 = x1 - delta;  a1 = x1;
x_d_l = 2*(x_d_left-a0)/(a1-a0) - 1;   % x_d_l in [-1,1] Q+1 points

P_d_l = zeros(Q+1, m+1);
P_d_l(:,1) = 1/sqrt(beta+1); 
P_d_l(:,2) = alpha(1)*x_d_l.*P_d_l(:,1);
for ii = 3: m+1
    P_d_l(:,ii) = alpha(ii-1)*x_d_l.*P_d_l(:,ii-1) - alpha(ii-1)/alpha(ii-2)*P_d_l(:,ii-2);
end

%-------- Get the even and odd k ----------------------------
if mod(mm, 2) == 0
    t_g = -mm/2+1:mm/2;   %t_g=g_mm
else
    t_g = -(mm-1)/2:(mm-1)/2;
end

if mod(t_g(1), 2) == 0
    t_g_even = t_g(1:2:mm);
    t_g_odd  = t_g(2:2:mm);
else
    t_g_even = t_g(2:2:mm);
    t_g_odd  = t_g(1:2:mm);
end

%------------ matrix A (for calculating a_hat) -------------------------
for ii = 1:Q+1
    angle = pi/(d+delta)*x_d_left(ii);
%     angle = pi/(d_len+delta)*x_d_l(ii);
    for jj = 1:length(t_g_even)
        A(ii,jj) = exp(1i*angle*t_g_even(jj));  
    end
    for jj = 1:length(t_g_odd)
        B(ii,jj) = exp(1i*angle*t_g_odd(jj));
    end
end

[U_even, S_even, V_even] = svd(A,0);
[U_odd,  S_odd,  V_odd ] = svd(B,0);
S_even_inv = inv(S_even);
S_odd_inv = inv(S_odd);
svd_tol = 1e-11; 
for ii = 1:length(t_g_even)
    if(abs(S_even_inv(ii, ii)) > 1/svd_tol)
        S_even_inv(ii, ii) = 0;
    end
end
for ii = 1:length(t_g_odd)
    if(abs(S_odd_inv(ii, ii)) > 1/svd_tol)
        S_odd_inv(ii, ii) = 0;
    end
end
a_hat = zeros(length(t_g_odd), m+1);
b_hat = zeros(length(t_g_odd), m+1);
for ii=1:m+1
    a_hat(:,ii) = V_even*S_even_inv*(U_even')*P_d_l(:,ii);
    b_hat(:,ii) = V_odd *S_odd_inv *(U_odd' )*P_d_l(:,ii);
end
a_hat = conj(a_hat);
b_hat = conj(b_hat);

f_even = zeros(gamma-2, m+1);
f_odd  = zeros(gamma-2, m+1);
for ii = 1:m+1
    for jj = 1:2*(beta+gamma-1)
        angle = pi/(d+delta)*x_d(jj);
        f_even(jj,ii) = sum(a_hat(:,ii)'.*exp(1i*angle*t_g_even));
        f_odd (jj,ii) = sum(b_hat(:,ii)'.*exp(1i*angle*t_g_odd ));
    end
end

f_match = zeros(2*(beta+gamma-1),1);
for ii = 1:2*(beta+gamma-1)
    f_match(ii,1) = sum((a_left(:,1)+a_right(:,1))'.*f_even(ii,:) +...
                         (a_left(:,1)-a_right(:,1))'.*f_odd(ii,:))/2;
end

%--------- periodic function by the FC, f_c---------------------------------------------------
M = N+gamma-1-1;  %M = sN-1?
f_fc = zeros(M+1,1);
f_fc(1:N,1) = f_new(1:N,1);  
f_fc(N+1:M,1) = real(f_match(beta+2:beta+gamma-1,1));
f_fc(M+1,1) = f_new(1,1); 

%% Fourier Coefficients 
if rem(sN1,2) == 0  %even expansion:M=evem number M+1=sN
    j  = 0:sN1;
    z1 = 2*pi*dx/res*j/bx;   %2*pi*j/M;
    j  = 0:M;
    z2 = 2*pi*dx*j/bx;    %2*pi*j/M;
    k1 = -M/2:M/2;k1=k1';
    A1 = (1/sN1)*exp(-1i*k1*z1);
    A1(1,:)   = A1(1,:)/2;
    A1(M+1,:) = A1(M+1,:)/2;
    f_hat = zeros(M+1,1);
    for ii = 1:M+1
        for jj = 1:M
            f_hat(ii,1) = f_hat(ii,1)+f_fc(jj,1)*exp(-1i*k1(ii)*z2(jj));
        end
    end
    f_hat(1,1)   = f_hat(1,1)/(2*M);
    f_hat(M+1,1) = f_hat(M+1,1)/(2*M);
    f_hat(2:M,1) = f_hat(2:M,1)/M;
else        %odd expansion:new N1=odd number
    j  = 0:sN1-1;
    z1 = 2*pi*dx/res*j/bx;
    j  = 0:M-1;
    z2 = 2*pi*dx*j/bx;
    k1 = -(M-1)/2:(M-1)/2;k1=k1';
    A1 = (1/sN1)*exp(-1i*k1*z1);
    f_hat = zeros(M,1);    %zeros(M+1,1);
    for ii = 1:M   %M+1
        for jj = 1:M
            f_hat(ii,1) = f_hat(ii,1)+f_fc(jj,1)*exp(-1i*k1(ii)*z2(jj));
        end
    end
    f_hat(:,1) = f_hat(:,1)/M;
end

%% Jump Function Approximation, assuming uniform points
do_periodic_bd = 0; % Does Local Edge detection using periodic on non-perodic conditions
m_1=2;
m_1=1;

C_sub = ones(m_1+1,1);
for j = 1:m_1+1
    for js =1:m_1+1
        if js ~= j
            C_sub(j) = C_sub(j)/(j-js);
        end
    end
end
C_sub = C_sub/sum(C_sub(1:floor(m_1/2)+1)); % Normalize Coefficients

% C is the edge detection matrix    
if rem(sN1,2)==0
    C = zeros(sN1,sN1+1);     
    for j=1:sN1
        inds = j-floor(m_1/2):j+m_1-floor(m_1/2);
        if min(inds)<1
            gap=1-min(inds);
            inds = inds+gap;  
        elseif max(inds)>sN1+1
            gap=max(inds)-(sN1+1);
            inds = inds-gap;
        end
        C(j,inds) = C_sub;
    end
else
    C = zeros(sN1,sN1);% C = zeros(n+gamma+1); %C will be our first edge detection matrix        
    for j=1:sN1
        inds = j-floor(m_1/2):j+m_1-floor(m_1/2);
        if min(inds)<1
            gap=1-min(inds);
            inds = inds+gap;  
        elseif max(inds)>sN1
            gap=max(inds)-sN1;
            inds = inds-gap;
        end
        C(j,inds) = C_sub;
    end
end

%------------L1 regularization(CVX)---------------
if rem(sN1,2)==0
    frecon = zeros(N1,1);
    cvx_begin quiet
    %         cvx_precision(0.01)
        variable f(sN1+1) 
        minimize(norm(A1*f-f_hat(:,1),2)+ lambda*norm(C*f,1));
    cvx_end
    frecon(1:N1,1) = f(1:N1);
else
    frecon = zeros(N1,1);
    cvx_begin quiet
        variable f(sN1) 
        minimize(norm(A1*f-f_hat(:,1),2)+ lambda*norm(C*f,1));
    cvx_end
    frecon(1:N1,1) = f(1:N1);  
end

return
