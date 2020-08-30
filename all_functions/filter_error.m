% filter_error.m	26 june 1999

% this function performs proper two-sided filtering of error to aleviate edge effects
% true filter size is twice the size of one-sided filter!!!

function fdata=filter_error(data,filt_onesided,ind)

% data is just one column
if length(data(:,1))==1
   data=data';
end
ld=length(data);
one=ones(ld,1);
if ind==0
   scale1=filter(filt_onesided,1,one);
   scale2=flipud(filter(filt_onesided,1,flipud(one)));
   scale=(scale1+scale2)/2;
else
   % in this case we discard the influence of error in part of signal we are not interested at the moment
   q=find(data==0);	
   one(q)=0;
   scale1=filter(filt_onesided,1,one);
   scale2=flipud(filter(filt_onesided,1,flipud(one)));
   scale=(scale1+scale2)/2;   
end
fdata1=filter(filt_onesided,1,data);
fdata2=flipud(filter(filt_onesided,1,flipud(data)));
fdata=(fdata1+fdata2)/2;
fdata=fdata./scale;
return;