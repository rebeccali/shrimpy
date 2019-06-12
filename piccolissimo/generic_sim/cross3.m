function [ c ] = cross3( a,b )
%CROSS3(a,b) Fact vector cross product.
%   A and B.  That is, C = A x B.  A and B must be 3 element
%   vectors.
%
%   C = CROSS3(A,B) returns the cross product of A and B along the
%   first dimension of length 3.  A and B must be size 3x1.

c = [a(2,:).*b(3,:)-a(3,:).*b(2,:) 
    a(3,:).*b(1,:)-a(1,:).*b(3,:)
    a(1,:).*b(2,:)-a(2,:).*b(1,:)];

end

