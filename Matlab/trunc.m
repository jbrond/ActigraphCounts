function [ outd ] = trunc( data, min_value )

outd = data;
    
I = find(data(:,1)<min_value);
        
outd(I,1) = 0;


    
    