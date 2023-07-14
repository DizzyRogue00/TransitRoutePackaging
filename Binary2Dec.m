function [delta]=Binary2Dec(l,u,binary_delta)
bits=max(u-l);
bits=ceil(log2(bits+1));
row=size(binary_delta,1);
binary_delta=reshape(binary_delta.',bits,[]);
binary_delta=binary_delta.';
delta=num2str(binary_delta);
delta=bin2dec(delta);
delta=reshape(delta,[],row).';
delta=l+delta;
end