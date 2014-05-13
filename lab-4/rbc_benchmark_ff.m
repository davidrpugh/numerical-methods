% Usage:
%       out = rbc_benchmark_ff(params, y)
%   where
%       out    is a (9,9) matrix of the first order
%              derivatives of the static system residuals
%              columns correspond to endo variables in
%              the ordering as declared
%       params is a (12,1) vector of parameter values
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
% sigma_z
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

function out = rbc_benchmark_ff(params, y)
if size(y) ~= [9,1]
	error('Wrong size of y, must be [9,1]');
end
if size(params) ~= [12,1]
	error('Wrong size of params, must be [12,1]');
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
% sigma_z not used in the model
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

t1 = a1;
t6 = exp(a5);
t15 = a10 - a7;
t11 = a9 - a10;
t12 = t11 / a9;
t83 = t12 - a1;
t84 = a16 ^ t83;
t85 = t12 * t84;
t86 = t15 * t85;
t20 = a9 / t11;
t13 = a8 ^ t12;
t14 = a7 * t13;
t17 = a16 ^ t12;
t18 = t15 * t17;
t19 = t14 + t18;
t89 = t20 - a1;
t90 = t19 ^ t89;
t91 = t20 * t90;
t92 = t86 * t91;
t93 = t6 * t92;
t94 = -(t93);
t21 = t19 ^ t20;
t22 = t6 * t21;
t95 = -(t22);
t98 = a8 ^ t83;
t99 = t12 * t98;
t100 = a7 * t99;
t101 = t91 * t100;
t102 = t6 * t101;
t103 = -(t102);
t104 = -(a24);
t31 = -(a10);
t32 = t31 / a9;
t33 = a16 ^ t32;
t34 = t15 * t33;
t35 = t34 / t19;
t106 = -(t35);
t108 = t32 - a1;
t109 = a16 ^ t108;
t110 = t32 * t109;
t111 = t15 * t110;
t113 = t19 * t111;
t112 = t34 * t86;
t114 = t113 - t112;
t115 = t19 * t19;
t116 = t114 / t115;
t117 = a4 * t116;
t118 = -(t117);
t119 = t34 * t100;
t120 = -(t119);
t121 = t120 / t115;
t122 = a4 * t121;
t123 = -(t122);
t39 = a8 ^ t32;
t40 = a7 * t39;
t41 = t40 / t19;
t124 = -(t41);
t125 = t40 * t86;
t126 = -(t125);
t127 = t126 / t115;
t128 = a4 * t127;
t129 = -(t128);
t131 = a8 ^ t108;
t132 = t32 * t131;
t133 = a7 * t132;
t135 = t19 * t133;
t134 = t40 * t100;
t136 = t135 - t134;
t137 = t136 / t115;
t138 = a4 * t137;
t139 = -(t138);
t105 = -(a1);
t140 = -(a16);
t141 = -(t140);
t142 = -(a30);
t143 = -(t142);
t144 = -(a8);
t145 = -(t144);
t146 = -(a38);
t147 = -(t146);
t57 = a10 - a56;
t148 = -(t57);
t69 = a10 + a68;
t70 = t69 - a56;
t62 = -(a61);
t151 = t62 - a1;
t152 = a65 ^ t151;
t153 = t62 * t152;
t154 = a64 * t153;
t155 = t70 * t154;
t156 = -(t155);
t66 = a65 ^ t62;
t67 = a64 * t66;
t157 = -(t67);
t160 = a50 ^ t151;
t161 = t62 * t160;
t63 = a50 ^ t62;
t162 = -(t63);
t76 = -(a75);
t74 = a10 - a16;
t165 = t76 - a1;
t166 = t74 ^ t165;
t167 = t76 * t166;
t168 = t105 * t167;
t169 = a73 * t168;
t170 = a30 * t161;
t171 = -(t170);
% setting the output variable
out = zeros(9, 9);
out(1,4) = out(1,4) + t1; % Y(0)
out(1,8) = out(1,8) + t94; % L(0)
out(1,3) = out(1,3) + t95; % Z(0)
out(1,1) = out(1,1) + t103; % K(-1)
out(2,3) = out(2,3) + t1; % Z(0)
out(2,3) = out(2,3) + t104; % Z(-1)
out(3,4) = out(3,4) + t106; % Y(0)
out(3,6) = out(3,6) + t1; % W(0)
out(3,8) = out(3,8) + t118; % L(0)
out(3,1) = out(3,1) + t123; % K(-1)
out(4,4) = out(4,4) + t124; % Y(0)
out(4,8) = out(4,8) + t129; % L(0)
out(4,7) = out(4,7) + t1; % r(0)
out(4,1) = out(4,1) + t139; % K(-1)
out(5,4) = out(5,4) + t105; % Y(0)
out(5,6) = out(5,6) + t141; % W(0)
out(5,8) = out(5,8) + t143; % L(0)
out(5,9) = out(5,9) + t1; % zero_profit(0)
out(5,7) = out(5,7) + t145; % r(0)
out(5,1) = out(5,1) + t147; % K(-1)
out(6,5) = out(6,5) + t1; % I(0)
out(6,6) = out(6,6) + t140; % W(0)
out(6,8) = out(6,8) + t142; % L(0)
out(6,2) = out(6,2) + t1; % C(0)
out(6,7) = out(6,7) + t144; % r(0)
out(6,1) = out(6,1) + t146; % K(-1)
out(7,5) = out(7,5) + t105; % I(0)
out(7,1) = out(7,1) + t1; % K(0)
out(7,1) = out(7,1) + t148; % K(-1)
out(8,2) = out(8,2) + t156; % C(1)
out(8,7) = out(8,7) + t157; % r(1)
out(8,2) = out(8,2) + t161; % C(0)
out(9,6) = out(9,6) + t162; % W(0)
out(9,8) = out(9,8) + t169; % L(0)
out(9,2) = out(9,2) + t171; % C(0)
