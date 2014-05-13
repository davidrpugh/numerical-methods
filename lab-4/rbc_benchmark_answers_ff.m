% Usage:
%       out = rbc_benchmark_answers_ff(params, y)
%   where
%       out    is a (9,9) matrix of the first order
%              derivatives of the static system residuals
%              columns correspond to endo variables in
%              the ordering as declared
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

function out = rbc_benchmark_answers_ff(params, y)
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
t436 = -(a24);
t31 = -(a10);
t32 = t31 / a9;
t33 = a16 ^ t32;
t34 = t15 * t33;
t35 = t34 / t19;
t438 = -(t35);
t440 = t32 - a1;
t441 = a16 ^ t440;
t442 = t32 * t441;
t443 = t15 * t442;
t445 = t19 * t443;
t444 = t34 * t86;
t446 = t445 - t444;
t447 = t19 * t19;
t448 = t446 / t447;
t449 = a4 * t448;
t450 = -(t449);
t451 = t34 * t100;
t452 = -(t451);
t453 = t452 / t447;
t454 = a4 * t453;
t455 = -(t454);
t39 = a8 ^ t32;
t40 = a7 * t39;
t41 = t40 / t19;
t1116 = -(t41);
t1117 = t40 * t86;
t1118 = -(t1117);
t1119 = t1118 / t447;
t1120 = a4 * t1119;
t1121 = -(t1120);
t1123 = a8 ^ t440;
t1124 = t32 * t1123;
t1125 = a7 * t1124;
t1127 = t19 * t1125;
t1126 = t40 * t100;
t1128 = t1127 - t1126;
t1129 = t1128 / t447;
t1130 = a4 * t1129;
t1131 = -(t1130);
t437 = -(a1);
t1678 = -(a16);
t1679 = -(t1678);
t1680 = -(a30);
t1681 = -(t1680);
t1682 = -(a8);
t1683 = -(t1682);
t1684 = -(a38);
t1685 = -(t1684);
t57 = a10 - a56;
t1687 = -(t57);
t69 = a10 + a68;
t70 = t69 - a56;
t62 = -(a61);
t1690 = t62 - a1;
t1691 = a65 ^ t1690;
t1692 = t62 * t1691;
t1693 = a64 * t1692;
t1694 = t70 * t1693;
t1695 = -(t1694);
t66 = a65 ^ t62;
t67 = a64 * t66;
t1696 = -(t67);
t1699 = a50 ^ t1690;
t1700 = t62 * t1699;
t63 = a50 ^ t62;
t1765 = -(t63);
t76 = -(a75);
t74 = a10 - a16;
t1768 = t76 - a1;
t1769 = t74 ^ t1768;
t1770 = t76 * t1769;
t1771 = t437 * t1770;
t1772 = a73 * t1771;
t1773 = a30 * t1700;
t1774 = -(t1773);
% setting the output variable
out = zeros(9, 9);
out(1,4) = out(1,4) + t1; % Y(0)
out(1,8) = out(1,8) + t94; % L(0)
out(1,3) = out(1,3) + t95; % Z(0)
out(1,1) = out(1,1) + t103; % K(-1)
out(2,3) = out(2,3) + t1; % Z(0)
out(2,3) = out(2,3) + t436; % Z(-1)
out(3,4) = out(3,4) + t438; % Y(0)
out(3,6) = out(3,6) + t1; % W(0)
out(3,8) = out(3,8) + t450; % L(0)
out(3,1) = out(3,1) + t455; % K(-1)
out(4,4) = out(4,4) + t1116; % Y(0)
out(4,8) = out(4,8) + t1121; % L(0)
out(4,7) = out(4,7) + t1; % r(0)
out(4,1) = out(4,1) + t1131; % K(-1)
out(5,4) = out(5,4) + t437; % Y(0)
out(5,6) = out(5,6) + t1679; % W(0)
out(5,8) = out(5,8) + t1681; % L(0)
out(5,9) = out(5,9) + t1; % zero_profit(0)
out(5,7) = out(5,7) + t1683; % r(0)
out(5,1) = out(5,1) + t1685; % K(-1)
out(6,5) = out(6,5) + t1; % I(0)
out(6,6) = out(6,6) + t1678; % W(0)
out(6,8) = out(6,8) + t1680; % L(0)
out(6,2) = out(6,2) + t1; % C(0)
out(6,7) = out(6,7) + t1682; % r(0)
out(6,1) = out(6,1) + t1684; % K(-1)
out(7,5) = out(7,5) + t437; % I(0)
out(7,1) = out(7,1) + t1; % K(0)
out(7,1) = out(7,1) + t1687; % K(-1)
out(8,2) = out(8,2) + t1695; % C(1)
out(8,7) = out(8,7) + t1696; % r(1)
out(8,2) = out(8,2) + t1700; % C(0)
out(9,6) = out(9,6) + t1765; % W(0)
out(9,8) = out(9,8) + t1772; % L(0)
out(9,2) = out(9,2) + t1774; % C(0)
