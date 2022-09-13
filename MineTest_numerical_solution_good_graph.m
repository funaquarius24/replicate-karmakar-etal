H = 1;
L = 1;
U = 1;



Da = 0.01; % Darcy number
Ar = H/L; % Aspect ratio
ls = 0.1; % Non dimensional slip length
M = 10;   % Viscosity ratio 
K = 1;    % Permeability ratio
F = 10;   % Non-dimensional inertial coeffcient
phi = pi / 5;

small_e = 1 / sqrt(Da);


a = sin(phi).^2 + K*cos(phi).^2;
b = sin(phi).^2 + sqrt(K)*cos(phi).^2;
c = (K-1)*sin(phi)*cos(phi);
d = sqrt((K-1))*sin(phi)*cos(phi);

y = linspace(0, 1, 10);
u = zeros(size(y));

u(end) = U;

P = F*b*small_e*u;

term0_1 = sqrt((a*small_e^2 + P) / M) .* y;
term0_2 = small_e * sqrt(a / M) .* y;

D_term_0 = (small_e*sqrt(a))/sqrt(M);
D_term_1 = (small_e*sqrt(a))/M;
D_term_2 = (small_e*a)/M;

D1_numerator = a^(5/2)*U*ls*small_e^3-ls*a^(3/2)*small_e-sqrt(M)*sinh(D_term_0);
D1_denominator = a*small_e^2*(a^(3/2)*cosh(D_term_0)*ls*small_e + sqrt(M)*sinh(D_term_1)*a);
D1 = D1_numerator / D1_denominator;

D2_numerator = sqrt(M)*(U*a*small_e^2 + cosh(D_term_2) -1);
D2_denominator = a*small_e^2*(a^(3/2)*cosh(D_term_0)*ls*small_e + sqrt(M)*sinh(D_term_1)*a);
D2 = D2_numerator / D2_denominator;

C_term0 = sqrt((a*small_e^2+P)/M);
C_term1 = (a*small_e^2+P);

C1_numerator = -sinh(C_term0).*sqrt(M) + sqrt(C_term1).*ls.*(U*a*small_e^2 + P*U - 1);
C1_denominator = C_term1.*(cosh(C_term0).*sqrt(C_term1).*ls + sqrt(M).*sinh(C_term0));
C1 = C1_numerator ./ C1_denominator;

C2_numerator = sqrt(M).*(cosh(C_term0) - 1 + C_term1*U);
C2_denominator = C_term1.*(cosh(C_term0).*sqrt(C_term1)*ls + sqrt(M)*sinh(C_term0));

C2 = C2_numerator ./ C2_denominator;



u_of_y = C1.*cosh(term0_1) + C2.*sinh(term0_1) + 1./(a*small_e^2 + P);

u_0_of_y = D1.*cosh(term0_2) + D2.*sinh(term0_2) + 1./(a*small_e^2);

u_of_y
u_0_of_y
term3 = ((2 * M * U + 1) * ls) / ( 2 * M * (ls + 1));


% u = term1 + term2 + term3;

