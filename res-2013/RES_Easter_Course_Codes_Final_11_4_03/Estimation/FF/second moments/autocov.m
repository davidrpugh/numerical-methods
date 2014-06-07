function a = autocov(tseries)
% AUTOCOV(TSERIES) sample autocovariance function
% TSERIES must be a column vector
% returns a column vector

n = length(tseries);
na = ceil(10*log10(n));
tscent = tseries - mean(tseries);
temp = zeros(n,na);
temp(:,1) = tscent.^2;
for j = 2:na
   temp(1:(n-j+1),j) = tscent(j:n);
   temp(:,j) = temp(:,j).*tscent;
end
a = (sum(temp))'/n;
                   