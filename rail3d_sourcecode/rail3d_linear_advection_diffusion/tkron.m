function A = tkron(B,C)
% TKRON Create the tensor Kronecker product.
%
% A = tkron(B,C) Forms the tensor Kronecker product of the two order-d
% tensors B and C.
% 
% S. Ragnarsson. Structured tensor computations: blocking, symmetries and
% Kronecker factorizations. PhD thesis, Cornell University, 2012.
%
% Modified the creation of A since the MATLAB toolbox extension mentioned
% in the thesis is no longer available (based on my investigating).
% 
% 

% if isnumeric(B), B = tensor(B); end
% if isnumeric(C), C = tensor(C); end
[n1,n2,n3] = size(B); [m1,m2,m3] = size(C); %needed in case B,C order-2.
n = [n1,n2,n3]; m = [m1,m2,m3]; N = n.*m;
d = length(n);
if length(m)~=d
    error('Tensor Kronecker products of different order tensors not supported.')
end
A = zeros(N);

for i1 = 1:n(1)
    for i2 = 1:n(2)
        for i3 = 1:n(3)
            A(1+(i1-1)*m(1):i1*m(1),1+(i2-1)*m(2):i2*m(2),1+(i3-1)*m(3):i3*m(3)) = B(i1,i2,i3)*C;
        end
    end
end
end


% function A = tkron(B,C)
% % TKRON Create the tensor Kronecker product.
% %
% % A = tkron(B,C) Forms the tensor Kronecker product of the two order-d
% % tensors B and C.
% if isnumeric(B), B = tensor(B); end
% if isnumeric(C), C = tensor(C); end
% n = size(B); m = size(C); N = n.*m;
% d = length(n);
% if length(m)~=d
%     error('Tensor Kronecker products of different order tensors not supported.')
% end
% 
% 
% %A = bltensor(zeros(N),m);
% 
% 
% 
% 
% 
% for i = 1:prod(n)
%     A(i) = B(i)*C;
% end
