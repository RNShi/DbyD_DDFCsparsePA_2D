

function n_more_aTV_x = function_aTV_extra_point_x(f_new1, n_new1, n_more_min, n_more_max, beta, n)

%% Test # of extended points should be used (TV)
K=1.5;
n_more_TV_x = zeros(n_new1,2);
TV = zeros(n_new1,2);
%top
for i=1:n_new1
    n_more_TV_x(i,1) = n_more_min;
    for j=1:n_more_min
        TV(i,1) = TV(i,1)+abs(f_new1(1+n_more_max-j+1,i)-f_new1(1+n_more_max-j,i));
    end
    for j=3:n_more_max
        TV(i,1)=TV(i,1)+abs(f_new1(1+n_more_max-j+1,i)-f_new1(1+n_more_max-j,i));
        if abs(f_new1(1+n_more_max-j+1,i)-f_new1(1+n_more_max-j,i))>TV(i,1)/(n_more_TV_x(i,1)+1)*K
            break;
        end
        n_more_TV_x(i,1)=n_more_TV_x(i,1)+1;
    end
end

%bottom
for i=1:n_new1
    n_more_TV_x(i,2) = n_more_min;
    for j=1:n_more_min
        TV(i,2) = TV(i,2)+abs(f_new1(n_more_max+n+j,i)-f_new1(n_more_max+n+j-1,i));
    end
    for j=3:n_more_max
        TV(i,2)=TV(i,2)+abs(f_new1(n_more_max+n+j,i)-f_new1(n_more_max+n+j-1,i));
        if abs(f_new1(n_more_max+n+j,i)-f_new1(n_more_max+n+j-1,i))>TV(i,2)/(n_more_TV_x(i,2)+1)*K
            break;
        end
        n_more_TV_x(i,2)=n_more_TV_x(i,2)+1;
    end
end

%% Test # of extended points should be used (aTV)
n_more_aTV_x=zeros(n_new1,2);
aTV_min=zeros(n_new1,2);
aTV_temp=zeros(n_new1,2);
%top
for i=1:n_new1
    n_more_aTV_x(i,1)=n_more_min;
    for j=1:beta
        aTV_min(i,1)= aTV_min(i,1)+abs(f_new1(1+n_more_max-n_more_min+j, i)-...
                                       f_new1(1+n_more_max-n_more_min+j-1, i));
    end
    if n_more_TV_x(i,1) == n_more_min
        continue;
    else
        for j=1:n_more_TV_x(i,1)-n_more_min
            aTV_temp(i,1)=0;
            for k=1:beta
                aTV_temp(i,1)= aTV_temp(i,1)+abs(f_new1(1+n_more_max-n_more_min+k-j, i)-...
                                                 f_new1(1+n_more_max-n_more_min+k-1-j, i));
            end
            if aTV_temp(i,1)<=aTV_min(i,1)
                aTV_min(i,1)=aTV_temp(i,1);
                n_more_aTV_x(i,1)=j+n_more_min;
            end
        end
    end
end

%bottom
for i=1:n_new1
    n_more_aTV_x(i,2)=n_more_min;
    for j=1:beta
        aTV_min(i,2)= aTV_min(i,2)+abs(f_new1(n_more_max+n-beta+n_more_min+j, i)-...
                                       f_new1(n_more_max+n-beta+n_more_min+j-1, i));
    end
    if n_more_TV_x(i,2) == n_more_min
        continue;
    else
        for j=1:n_more_TV_x(i,2)-n_more_min
            aTV_temp(i,2)=0;
            for k=1:beta
                aTV_temp(i,2)= aTV_temp(i,2)+abs(f_new1(n_more_max+n-beta+n_more_min+k+j, i)-...
                                                 f_new1(n_more_max+n-beta+n_more_min+k-1+j, i));
            end
            if aTV_temp(i,2)<=aTV_min(i,2)
                aTV_min(i,2)=aTV_temp(i,2);
                n_more_aTV_x(i,2)=j+n_more_min;
            end
        end
    end
end

return