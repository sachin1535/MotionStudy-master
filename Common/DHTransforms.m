function [H]=DHTransforms(thetas,alphas,disps,offsets)
%
%   Return the homogeneous transforms for the DH parameters 
%   
%   Input variables:
%       thetas(i)   -   rotation of link i for i=1...ndof
%       alphas(i)   -   twist of link i for i=1...ndof
%       disps(i)    -   displacement of link i for i=1...ndofs
%       offsets(i)  -   offset of link i for i=1...ndofs
%
%   Output variables
%       As(m,n,i)       -   is the 4x4xndof array of homogeneous transforms
%
[ndof, ~]=size(thetas);
As=zeros(4,4);
H = eye(4,4);
if isa(thetas,'sym')==1
    As=sym(As);
end
%
for i=1:ndof
    ct=cos(thetas(i));
    st=sin(thetas(i));
    ca=cos(alphas(i));
    sa=sin(alphas(i));
    a=offsets(i);
    d=disps(i);
    %
    %   rotation
    %
    As(1,1)=ct;
    As(1,2)=-st*ca;
    As(1,3)=st*sa;
    As(2,1)=st;
    As(2,2)=ct*ca;
    As(2,3)=-ct*sa;
    As(3,1)=0;
    As(3,2)=sa;
    As(3,3)=ca;
    %
    %   displacement
    %
    As(1,4)=a*ct;
    As(2,4)=a*st;
    As(3,4)=d;
    As(4,4)=1;
    
    H = H*As;
end
    