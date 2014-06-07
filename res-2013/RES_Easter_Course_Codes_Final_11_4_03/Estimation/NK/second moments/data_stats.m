%data statistics

clc
close all
clear all

load us_data
%std
fprintf('STD Output ')
std(dy(5:end))
fprintf('STD Inflation ')
std(pinfobs(5:end))
fprintf('STD Interest rate ')
std(robs(5:end))

fprintf('Cross-Correlation Inflation-Output ')
corrcoef(pinfobs(5:end),dy(5:end))

fprintf('Cross-Correlation Interest Rate-Output ')
corrcoef(robs(5:end),dy(5:end))

fprintf('Autocorrelation Output ')
acfcomp(dy(5:end))

fprintf('Autocorrelation Inflation ')
acfcomp(pinfobs(5:end))

fprintf('Autocorrelation Interest Rate ')
acfcomp(robs(5:end))