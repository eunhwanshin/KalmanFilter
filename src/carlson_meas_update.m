function [x,C] = carlson_meas_update(x,C,z,H,R)
% function [x,C] = carlson_meas_update(x,C,z,H,R)
% Ref: Table 6.11, Grewal and Andrews (2001), p. 241 
n = length(x);

alpha = R;
inno  = z;
w = zeros(n,1);
for j=1:n,
    inno = inno - H(j) * x(j);
	
	% sigma stores v_j
    sigma = 0;
    for i=1:j,
        sigma = sigma + C(i,j) * H(i);
    end
    
	% \beta = R + \sum_{k=1} ^{j-1} v_k^2
	beta = alpha;
	
	% \alpha = R + \sum_{k=1} ^{j} v_k^2
    alpha = alpha + sigma^2;
	
	% \gamma = \sqrt{R + \sum_{k=1} ^{j-1} v_k^2} \sqrt{R + \sum_{k=1} ^{j} v_k^2}
    gamma = sqrt(alpha*beta);
	
	% \eta = \sqrt{R + \sum_{k=1} ^{j-1} v_k^2}/\sqrt{R + \sum_{k=1} ^{j} v_k^2}
    eta = beta/gamma;
	
	% \zeta  =  \frac{v_j }{ \sqrt{R + \sum_{k=1} ^{j-1} v_k^2} \sqrt{R + \sum_{k=1} ^{j} v_k^2}}
    zeta = sigma/gamma;

    w(j) = 0;
    for i=1:j,
        tau = C(i,j);
        C(i,j) = eta*C(i,j) - zeta*w(i);
        w(i) = w(i)+tau*sigma;
    end
    
end

epsilon = inno/alpha;
for i=1:n,
    x(i) = x(i) + w(i)*epsilon;
end

