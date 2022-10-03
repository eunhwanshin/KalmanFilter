function [Um,Dm] = ud_cov_predict(PhiU, D, G, DQ)
% Function [U,D] = ud_cov_predict(PhiU, D, G, DQ) implements UD-factorized  
% Kalman covariance prediction; see Grewal and Andrews (2001, p 251).
%
n = length(D);
p = length(DQ);

Um = zeros(n);
Dm = zeros(n,1);

for i=n:-1:1,
    % compute weighted norm-squared of i-th row of [PhiU, G]
    sigma = 0;
    for j=1:n,
        sigma = sigma + PhiU(i,j)^2 * D(j);
    end
    for j=1:p,
        sigma = sigma + G(i,j)^2 * DQ(j);
    end
    Dm(i) = sigma;
    Um(i,i) = 1;
    for j=1:i-1,
	    % compute weighted inner product between j-th row and i-th row
        sigma = 0;
        for k=1:n,
            sigma = sigma + PhiU(i,k)*D(k)*PhiU(j,k);
        end
        for k=1:p,
            sigma = sigma + G(i,k)*DQ(k)*G(j,k);
        end
		
		% Compute predicted U
		
        Um(j,i) = sigma / Dm(i);
		
		% orthogonalize j-th row w.r.t. i-th row 
        for k=1:n,
            PhiU(j,k) = PhiU(j,k)-Um(j,i)*PhiU(i,k);
        end
        for k=1:p,
            G(j,k) = G(j,k) - Um(j,i)*G(i,k);
        end
    end
        
end