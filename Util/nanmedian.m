function y=nanmedian(x)
%NANMEDIAN  Median ignoring NaNs 
%
%Syntax: y=nanmedian(x)
% x is a matrix and y its median ignoring NaN

%Author: Caroline Lafleur, physical oceanography
%Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%email: lafleurc@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%June 1999; Last revision: 18-Jun-1999 CL

[m,n]=size(x);
if m==1, x=x(:); m=n; n=1; end
y=repmat(nan,1,n);
for i=1:n
  I=find(~isnan(x(:,i)));
  if ~isempty(I)
     y(i)=median(x(I,i));
  end
end

