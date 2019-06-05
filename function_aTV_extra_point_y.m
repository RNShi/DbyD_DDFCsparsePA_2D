

function n_more_aTV_y = function_aTV_extra_point_y(F2, n_new1, n_more_min, n_more_max, beta, a1, a2, b1, b2)

%% Test # of extended points should be used (TV)
K=1.5;
n_more_TV_y = zeros(n_new1,2);
TV = zeros(n_new1,2);
%left
for i=1:n_new1
    n_more_TV_y(i,1) = n_more_min;
    for j=1:n_more_min
        TV(i,1) = TV(i,1)+abs(F2(a1-1+i,b1-j+1)-F2(a1-1+i,b1-j));
    end
    for j=2:n_more_max-1
        TV(i,1)=TV(i,1)+abs(F2(a1-1+i,b1-j)-F2(a1-1+i,b1-j-1));
        if abs(F2(a1-1+i,b1-j)-F2(a1-1+i,b1-j-1))>TV(i,1)/(n_more_TV_y(i,1)+1)*K
            break;
        end
        n_more_TV_y(i,1)=n_more_TV_y(i,1)+1;
    end
end

%right
for i=1:n_new1
    n_more_TV_y(i,2) = n_more_min;
    for j=1:n_more_min
        TV(i,2) = TV(i,2)+abs(F2(a1-1+i,b2+j)-F2(a1-1+i,b2+j-1));
    end
    for j=2:n_more_max-1
        TV(i,2)=TV(i,2)+abs(F2(a1-1+i,b2+j+1)-F2(a1-1+i,b2+j));
        if abs(F2(a1-1+i,b2+j+1)-F2(a1-1+i,b2+j))>TV(i,2)/(n_more_TV_y(i,2)+1)*K
            break;
        end
        n_more_TV_y(i,2)=n_more_TV_y(i,2)+1;
    end
end

%% Test # of extended points should be used (aTV)
n_more_aTV_y=zeros(n_new1,2);
aTV_min=zeros(n_new1,2);
aTV_temp=zeros(n_new1,2);
%left
for i=1:n_new1
    n_more_aTV_y(i,1)=n_more_min;
    for j=1:beta
        aTV_min(i,1)= aTV_min(i,1)+abs(F2(a1-1+i, b1+n_more_max-1+j+1-n_more_min)-...
                                       F2(a1-1+i, b1+n_more_max-1+j-n_more_min));
    end
    if n_more_TV_y(i,1) == n_more_min
        continue;
    else
        for j=1:n_more_TV_y(i,1)-n_more_min
            aTV_temp(i,1)=0;
            for k=1:beta
                aTV_temp(i,1)= aTV_temp(i,1)+abs(F2(a1-1+i, b1+n_more_max-1+k-j+1-n_more_min)-...
                                                 F2(a1-1+i, b1+n_more_max-1+k-j-n_more_min));
            end
            if aTV_temp(i,1)<=aTV_min(i,1)
                aTV_min(i,1)=aTV_temp(i,1);
                n_more_aTV_y(i,1)=j+n_more_min;
            end
        end
    end
end
%right
for i=1:n_new1
    n_more_aTV_y(i,2)=n_more_min;
    for j=1:beta
        aTV_min(i,2)= aTV_min(i,2)+abs(F2(a1-1+i, b2+n_more_max-beta+j+n_more_min)-...
                                       F2(a1-1+i, b2+n_more_max-beta+j+n_more_min-1));
    end
    if n_more_TV_y(i,2) == n_more_min
        continue;
    else
        for j=1:n_more_TV_y(i,2)-n_more_min
            aTV_temp(i,2)=0;
            for k=1:beta
                aTV_temp(i,2)= aTV_temp(i,2)+abs(F2(a1-1+i, b2+n_more_max-beta+k+j+n_more_min)-...
                                                 F2(a1-1+i, b2+n_more_max-beta+k+j+n_more_min-1));
            end
            if aTV_temp(i,2)<=aTV_min(i,2)
                aTV_min(i,2)=aTV_temp(i,2);
                n_more_aTV_y(i,2)=j+n_more_min;
            end
        end
    end
end

return