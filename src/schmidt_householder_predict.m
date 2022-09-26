function [A,B] =  schmidt_householder_predict(A, B)
% function [A,B] =  schmidt_householder_predict(A, B)
% Inputs:  A = Phi C(+)
%          B = G Cq
% Outputs: A = C(-)
%          B = Zeroed
% Ref.: Grewal, M. S. and Andrews, A. P. (2001), 
%       Kalman Filtering: Theory and Practice, 
%       John Wiley & Sons, Inc., p. 244.
n = size(A, 2);
r = size(B, 2);

v = zeros(n,1);
w = zeros(r,1);

for k=n:-1:1,
    sigma = 0;
    for j=1:r,
        sigma = sigma + B(k,j)^2;
    end
    for j=1:k,
        sigma = sigma + A(k,j)^2;
    end
    alpha = sqrt(sigma);
    sigma = 0;
    for j=1:r,
        w(j) = B(k,j);
        sigma = sigma + w(j)^2;
    end
    for j=1:k,
        if j==k
            v(j)=A(k,j) - alpha;
        else
            v(j) = A(k,j);
        end
        sigma = sigma + v(j)^2;
    end
    alpha = 2/sigma;
	
    for i=1:k,
        sigma = 0;
        for j=1:r,
            sigma = sigma + B(i,j)*w(j);
        end
        for j=1:k,
            sigma = sigma + A(i,j)*v(j);
        end
        beta = alpha * sigma;
        for j=1:r,
            B(i,j) = B(i,j) - beta*w(j);
        end
        for j=1:k,
            A(i,j) = A(i,j) - beta*v(j);
        end
    end
    
end