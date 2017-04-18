function y=nanmin(x,dim)
%NANMIN  Minimum ignoring NaNs 
%
%Syntax: y=nanmin(x)
% x is a matrix and y its minimum ignoring NaN

%Author: Caroline Lafleur, physical oceanography
%Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%email: lafleurc@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%June 1999; Last revision: 18-Jun-1999 CL

%[m,n]=size(x);
%if m==1, x=x(:); m=n; n=1; end
%y=repmat(nan,1,n);
%for i=1:n
%  I=find(~isnan(x(:,i)));
%  if ~isempty(I)
     %y(i)=min(x(I,i));
     %end
     %end

[m,n]=size(x);
if nargin==1, 
  % Determine which dimension MEAN will use
  if m==1, x=x(:);
     dim=1; m=n; n=1;
  else
     dim = min(find(size(x)~=1));
  end
elseif nargin==2
    if dim==2 & m==1
    x=x(:);
    dim=1; m=n; n=1;
    end
end
switch dim
case 1
   y=repmat(nan,1,n);
   for i=1:n
     I=find(~isnan(x(:,i)));
     if ~isempty(I)
        y(i)=min(x(I,i),[],dim);
     end
   end
case 2
   y=repmat(nan,m,1);
   for i=1:m
     I=find(~isnan(x(i,:)));
     if ~isempty(I)
        y(i)=min(x(i,I),[],dim);
     end
   end    
end
