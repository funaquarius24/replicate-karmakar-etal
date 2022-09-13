H = 1;
L = 1;
U = 1;

Da = 3; % Darcy number
Ar = H/L; % Aspect ratio
ls = 0.1; % Non dimensional slip length
M = 10;   % Viscosity ratio 
phi = pi / 5;

y = linspace(0, 1, 10);
u = zeros(size(y));

u(end) = U;

term1 = y.^2 / 2 * M * -1;
term2 = ((2 * M * U + 1) * y) / ( 2 * M * (ls + 1));
term3 = ((2 * M * U + 1) * ls) / ( 2 * M * (ls + 1));


u = term1 + term2 + term3;

