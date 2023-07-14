function [binary_delta]=Dec2Binary(l,u,delta)
bits=max(u-l);
bits=ceil(log2(bits+1));
delta=delta-l;
delta_row=size(delta,1);
binary_delta=dec2bin(delta.',bits);
binary_delta=boolean(binary_delta-'0');
binary_delta=double(binary_delta);
binary_delta=reshape(binary_delta.',[],delta_row);
binary_delta=binary_delta.';
end