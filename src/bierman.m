function [x,U,D,K] = bierman(x,U,D,z,H,R)
% function [x,U,D] = bierman(x,U,D,z,H,R) applies UD measurement update.
n = length(x);
v = zeros(n,1);
w = zeros(n,1);
delta = z;
for j=1:n,
    delta = delta - H(j)*x(j);
	v(j) = H(j);
	for i=1:j-1,
	    v(j) = v(j) + U(i,j)*H(i);
	end
end
sigma = R;
for j=1:n,
    nu=v(j);
	v(j) = v(j) * D(j,j);
	w(j) = v(j);
	for i=1:j-1,
	    tau = U(i,j) * v(j);
		U(i,j) = U(i,j) -nu *w(i)/sigma;
		w(i) = w(i) + tau;
	end
	D(j,j) = D(j,j)*sigma;
	sigma = sigma + nu*v(j);
	D(j,j) = D(j,j)/sigma;
end
%epsilon = delta/sigma;
K = w / sigma
for i=1:n,
    x(i) = x(i)+K(i)*delta;
end