function K = formK(x1, x2, param)

% form a kernel given parameters and input x

% 1. unpack parameters
theta0 = param(1);
theta1 = param(2);

% 2. compute the l-2 norm
x_dim = length(x1(1,:));
nsamps_x1 = size(x1, 1);
nsamps_x2 = size(x2, 1);

xn_i = repmat(x1, 1, nsamps_x2);
xn = reshape(xn_i', x_dim, []);
xn = xn';
xm = repmat(x2, nsamps_x1, 1);

norm_vec = sum(abs(xn-xm).^2, 2);
norm_mat_tr = reshape(norm_vec, nsamps_x2, nsamps_x1);
norm_mat = norm_mat_tr';
K =  theta0*exp(-.5/theta1.*norm_mat);



