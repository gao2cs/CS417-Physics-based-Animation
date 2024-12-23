clear
clc
close all

% Matrices and Vectors:
I = eye(3); % 3x3 Identity matrix
B = [-I, I]; % Define B explicitly as a 3x6 matrix
q = sym('q', [6, 1], 'real'); % Declare q as a real 6x1 vector

% Scalar Variables:
syms lres real % Declare lres as a real scalar
syms k real    % Declare k as a real scalar

% Expression for V:
V = (1/2) * k * (sqrt(q' * B' * B * q) - lres)^2;

% Gradient of V = [dV/dq1, dV/dq2, ..., dV/dq6]^T
dV_dq = gradient(V, q);

% Hessian of V (Second derivatives w.r.t q):
dV2_dq2 = hessian(V, q);

% Display results
disp('Gradient of V:');
disp(dV_dq);

disp('Hessian of V:');
disp(dV2_dq2);

% Convert to C code style for use in C/C++:
%ccode(dV_dq)    % Generates C code for the gradient
ccode(dV2_dq2)  % Generates C code for the Hessian
