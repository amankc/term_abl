function [ind,indMin]=directionaldist(p,pk,r1,r2,dirAngle,margin,dirType)
% 'directionaldist' extract points in 'pk' that are in a direction relative
% to reference point 'p' (Output: 'ind'). You can define distance interval 
% with 'r1' and 'r2', and also angle interval by 'dirAngle' and 'margin'. 
% It also extracts the nearest neighbor in the specified direction and
% distance interval (Output: 'indMin').
% Inputs 
%
% Data
% 'p': Reference point location 2x1
% 'pk': location of N points 2xN
%
% Parameters
% 'r1': max distance to 'p' (Euclidean distance is used)
% 'r2': min distance to 'p' (Euclidean distance is used)
% 'dirAngle': reltaive direction angle from point 'p' to 'pk' (in degrees)
% 'margin': angle margin, will be used as defining the angle interval
% [dirAngle+margin, dirAngle-margin]. If margin is zero, it is possible 
% that you will not find any points that are exactly at dirAngle direction.
% So it is good to give some margin. Or it can be used for other purposes.
% 'dirType': direction type 'none' or 'symmetric', just try it out!
%
% Outputs
%
% 'ind': indices of points of 'pk', which are in the direction of 
% [dirAngle+margin, dirAngle-margin] relative to reference point 'p',
% where max distance to 'p' is 'r1', and min distance to 'p' is 'r2'.
% 'indMin': index of nearest neighbor point in 'pk' to 'p' in the specified
% direction angle interval
% Note: This function uses a function called 'nearestneighbour' 
% which is downloaded from fileexchange linke below:
% (http://www.mathworks.com/matlabcentral/fileexchange/12574-nearestneighbour-m) 
% 
% Copyright 2016 Abdullah H. Ozcan
%--------------------------------------------------------------------------
if r2==0
    r2=eps;
end
if r2>r1
    error('r1 should be bigger than r2...')
end


idx1 = nearestneighbour(p, pk, 'r', r1);
idx2 = nearestneighbour(p, pk, 'r', r2);
% idx = rangesearch(pk',p',r1); idx1=cell2mat(idx);
% idx = rangesearch(pk',p',r2); idx2=cell2mat(idx);
idx=setdiff(idx1,idx2);
[THETA,~] = cart2pol(p(1)-pk(1,idx),p(2)-pk(2,idx));
th=ceil(THETA/pi*180);
th=mod(th+180,360);
ind1=find(abs(dirAngle-th)<margin);
ind2=find(abs(dirAngle-th)>(360-margin));
ind=union(ind1,ind2);
if strcmp(dirType,'symmetric')
    a2=mod(180+dirAngle,360);
    ind1=find(abs(a2-th)<margin);
    ind2=find(abs(a2-th)>(360-margin));
    ind3=union(ind1,ind2);
    ind=union(ind,ind3);
end
ind=idx(ind);
[~,indMin]=min((p(1)-pk(1,ind)).^2+(p(2)-pk(2,ind)).^2);
indMin=ind(indMin);
