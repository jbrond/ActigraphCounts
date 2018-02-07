function rs = runsum(data,len,threshold)
%
%
% Running sum
%
% return data on running sum
%
% Auther : Jan Brønd
%
N = size(data);
cnt = ceil (N(1) / len);

rs = zeros(cnt,N(2));

for k=1:N(2)
    for n=1:cnt
        rs(n,k) = 0;
        for p=1+len*(n-1):len*n

            if (p<=N(1) && data(p,k)>=threshold)
                rs(n,k) = rs(n,k) + data(p,k) - threshold;
            end; 
        end
    end
end