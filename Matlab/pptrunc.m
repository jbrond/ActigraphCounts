function [ outd ] = pptrunc( data, max_value)
%peak to peak truncation
outd = data;

N = length(data(1,:));

for n=1:N
    I = find(data(:,n) > max_value);
    outd(I,n) = max_value;
    I = find(data(:,n) < -max_value);
    outd(I,n) = -max_value;
end

