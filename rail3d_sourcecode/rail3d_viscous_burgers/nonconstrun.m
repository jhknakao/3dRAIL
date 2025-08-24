function [Uout,Gout,rank] = nonconstrun(Uin,Gin,tol)
[SU,SG,~] = mlsvd(Gin,tol); %truncate wrt tol
if numel(SU)==2 %if resulting hosvd/mlsvd is, e.g., size 1x1x1 (output as matrix)
    SU = {SU{1},SU{2},[1]};
end
sz = size(SG);
if numel(sz)==2
    rank = [sz(1),sz(2),1];
else
    rank = [sz(1),sz(2),sz(3)];
end
Uout = {Uin{1}*SU{1},Uin{2}*SU{2},Uin{3}*SU{3}};
Gout = SG;
[r1_out,r2_out,r3_out] = size(Gout);
rank = [r1_out,r2_out,r3_out];
end