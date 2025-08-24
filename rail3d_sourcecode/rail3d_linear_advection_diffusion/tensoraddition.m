function [U,G] = tensoraddition(U1,G1,U2,G2)
% Direct sum of two third-order tensors.

U = cell(1,3);
U{1} = [U1{1},U2{1}];
U{2} = [U1{2},U2{2}];
U{3} = [U1{3},U2{3}];

sz1 = size(G1);
sz2 = size(G2);
if numel(sz1)==2 %i.e., if G1 is size r1 x r2 x 1 (not an issue if 1 is in the first or second dimension)
    sz1 = [sz1(1),sz1(2),1];
end
if numel(sz2)==2 %i.e., if G2 is size r1 x r2 x 1 (not an issue if 1 is in the first or second dimension)
    sz2 = [sz2(1),sz2(2),1];
end

G = zeros(sz1(1)+sz2(1),sz1(2)+sz2(2),sz1(3)+sz2(3));
G(1:sz1(1),1:sz1(2),1:sz1(3)) = G1;
G(sz1(1)+1:end,sz1(2)+1:end,sz1(3)+1:end) = G2;
end