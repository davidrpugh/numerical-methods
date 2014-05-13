% Usage:
%       out = rbc_benchmark_answers_f(params, y)
%   where
%       out    is a (9,1) column vector of the residuals
%              of the static system
%       params is a (11,1) vector of parameter values
%              in the ordering as declared
%       y      is a (9,1) vector of endogenous variables
%              in the ordering as declared
%
% Created by Dynare++ v. 4.3.3

% params ordering
% =====================
% b
% alpha
% beta
% delta
% sigma
% theta
% omega
% rho_z
% KSS
% LSS
% YSS
%
% y ordering
% =====================
% K
% C
% Z
% Y
% I
% W
% r
% L
% zero_profit

function out = rbc_benchmark_answers_f(params, y)
if size(y) ~= [9,1]
	error('Wrong size of y, must be [9,1]');
end
if size(params) ~= [11,1]
	error('Wrong size of params, must be [11,1]');
end

% hardwired constants
a0 =            0;
a1 =            1;
a2 = NaN;
a3 =    1.1283792;
% numerical constants
a10 =            1;
% parameter values
a73 = params(1); % b
a7 = params(2); % alpha
a64 = params(3); % beta
a56 = params(4); % delta
a9 = params(5); % sigma
a61 = params(6); % theta
a75 = params(7); % omega
a24 = params(8); % rho_z
% KSS not used in the model
% LSS not used in the model
% YSS not used in the model
% exogenous variables to zeros
a27 = 0.0; % eps_z
% endogenous variables to y
a8 = y(1); % K
a55 = y(1); % K
a50 = y(2); % C
a65 = y(2); % C
a25 = y(3); % Z
a5 = y(3); % Z
a4 = y(4); % Y
a51 = y(5); % I
a30 = y(6); % W
a38 = y(7); % r
a68 = y(7); % r
a16 = y(8); % L
a44 = y(9); % zero_profit

t6 = exp(a5);
t11 = a9 - a10;
t12 = t11 / a9;
t13 = a8 ^ t12;
t14 = a7 * t13;
t15 = a10 - a7;
t17 = a16 ^ t12;
t18 = t15 * t17;
t19 = t14 + t18;
t20 = a9 / t11;
t21 = t19 ^ t20;
t22 = t6 * t21;
t23 = a4 - t22;
t26 = a24 * a25;
t28 = t26 + a27;
t29 = a5 - t28;
t31 = -(a10);
t32 = t31 / a9;
t33 = a16 ^ t32;
t34 = t15 * t33;
t35 = t34 / t19;
t36 = a4 * t35;
t37 = a30 - t36;
t39 = a8 ^ t32;
t40 = a7 * t39;
t41 = t40 / t19;
t42 = a4 * t41;
t43 = a38 - t42;
t45 = a16 * a30;
t46 = a4 - t45;
t47 = a8 * a38;
t48 = t46 - t47;
t49 = a44 - t48;
t52 = a50 + a51;
t53 = t45 + t47;
t54 = t52 - t53;
t57 = a10 - a56;
t58 = a8 * t57;
t59 = a51 + t58;
t60 = a55 - t59;
t62 = -(a61);
t63 = a50 ^ t62;
t66 = a65 ^ t62;
t67 = a64 * t66;
t69 = a10 + a68;
t70 = t69 - a56;
t71 = t67 * t70;
t72 = t63 - t71;
t74 = a10 - a16;
t76 = -(a75);
t77 = t74 ^ t76;
t78 = a73 * t77;
t79 = a30 * t63;
t80 = t78 - t79;
% setting the output variable
out = zeros(9, 1);
out(1) = t23;
out(2) = t29;
out(3) = t37;
out(4) = t43;
out(5) = t49;
out(6) = t54;
out(7) = t60;
out(8) = t72;
out(9) = t80;
