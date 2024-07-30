% Function determines if Start/Endpoint of terminus is above, below or on
% the Fjord wall boundary
%
% arguments: (input)
% Terminus Position - Variable that contains all terminus positions for glacier.
% Polygon - Fjord walls converted to polygon
% StartP - set to 1 for lower fjord walls, otherwise set to 0
% Endp - set to 1 for upper fjord walls, otherwise set to 0
%
%
%
% arguments: (output)
% InPolygon - Binary value, 0 = outside Polygon (Above Fjord wall), 1 = inside
% polygon (Below Fjord wall)
% OnPolygon - Bionary value, 0 = not on Fjord wall, 1 = on Fjord wall.      



function [InPolygon,OnPolygon,Number] = InPolygonCheck(TerminusPosition,Polygon, StartP,EndP)

TermPos = [TerminusPosition.X' TerminusPosition.Y'];

if StartP==1
Point1 = [TermPos(1,1) TermPos(1,2)];
end
if EndP==1
Point1 = [TermPos(end,1) TermPos(end,2)];
end
[InPolygon,OnPolygon] = inpolygon(Point1(:,1),Point1(:,2),Polygon.Vertices(:,1),Polygon.Vertices(:,2));
NoInPoly = inpolygon(TermPos(:,1),TermPos(:,2),Polygon.Vertices(:,1),Polygon.Vertices(:,2));
Number = numel(TermPos(NoInPoly));
end
