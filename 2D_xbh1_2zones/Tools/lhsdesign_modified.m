function [X_scaled,X_normalized]=lhsdesign_modified(n,min_ranges_p,max_ranges_p)
%lhsdesign_modified is a modification of the Matlab Statistics function lhsdesign.
%It might be a good idea to jump straight to the example to see what does
%this function do.
%The following is the description of lhsdesign from Mathworks documentation
% X = lhsdesign(n,p) returns an n-by-p matrix, X, containing a latin hypercube sample of n values on each of p variables.
%For each column of X, the n values are randomly distributed with one from each interval (0,1/n), (1/n,2/n), ..., (1-1/n,1), and they are randomly permuted.

%lhsdesign_modified provides a latin hypercube sample of n values of
%each of p variables but unlike lhsdesign, the variables can range between
%any minimum and maximum number specified by the user, where as lhsdesign
%only provide data between 0 and 1 which might not be very helpful in many
%practical problems where the range is not bound to 0 and 1
%
%Inputs: 
%       n: number of radomly generated data points
%       min_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the minimum value of each variable
%       max_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the maximum value of each variable
%Outputs
%       X_scaled: [nxp] matrix of randomly generated variables within the
%       min/max range that the user specified
%       X_normalized: [nxp] matrix of randomly generated variables within the
%       0/1 range 
%
%Example Usage: 
%       [X_scaled,X_normalized]=lhsdesign_modified(100,[-50 100 ],[20  300]);
%       figure
%       subplot(2,1,1),plot(X_scaled(:,1),X_scaled(:,2),'*')
%       title('Random Variables')
%       xlabel('X1')
%       ylabel('X2')
%       grid on
%       subplot(2,1,2),plot(X_normalized(:,1),X_normalized(:,2),'r*')
%       title('Normalized Random Variables')
%       xlabel('Normalized X1')
%       ylabel('Normalized X2')
%       grid on


p=length(min_ranges_p);
[M,N]=size(min_ranges_p);
if M<N
    min_ranges_p=min_ranges_p';
end
    
[M,N]=size(max_ranges_p);
if M<N
    max_ranges_p=max_ranges_p';
end

slope=max_ranges_p-min_ranges_p;
offset=min_ranges_p;

SLOPE=ones(n,p);
OFFSET=ones(n,p);

for i=1:p
    SLOPE(:,i)=ones(n,1).*slope(i);
    OFFSET(:,i)=ones(n,1).*offset(i);
end
X_normalized = lhsdesign(n,p);

X_scaled=SLOPE.*X_normalized+OFFSET;

