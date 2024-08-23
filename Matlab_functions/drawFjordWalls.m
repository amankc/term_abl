
function [UpperFW,LowerFW,UpperPolygon,LowerPolygon]=drawFjordWalls(Landsat8image, ReferenceMatrix,GlacierName,outFolder)
Red = Landsat8image(:,:,1);
Blue = Landsat8image(:,:,2);
Green = Landsat8image(:,:,3);
NormRed = (Red - min(Red(:)))/(max(Red(:)) - min(Red(:)));
NormBlue = (Blue - min(Blue(:)))/(max(Blue(:)) - min(Blue(:)));
NormGreen = (Green - min(Green(:)))/(max(Green(:)) - min(Green(:)));
NormRGB = double(cat(3,NormGreen,NormBlue,NormRed));
%%
f1=figure('units','normalized','outerposition',[0 0 1 1]);
mapshow(NormRGB,ReferenceMatrix)
waitfor(msgbox("Draw upper Fjordwall boundaries"));
DrawUpper =drawpolyline();
UpperFjordWall = DrawUpper.Position;
waitfor(msgbox("Draw lower Fjordwall boundaries"));
DrawLower = drawpolyline();
LowerFjordWall = DrawLower.Position;
close(f1)
%%
UfwX=[];
UfwY=[];
LfwX=[];
LfwY=[];

for i=1:length(UpperFjordWall(:,1))-1
    UfwX = [UfwX; linspace(UpperFjordWall(i,1),UpperFjordWall(i+1,1),100)'];
    UfwY = [UfwY; linspace(UpperFjordWall(i,2),UpperFjordWall(i+1,2),100)'];
end

for i=1:length(LowerFjordWall(:,1))-1
    LfwX = [LfwX; linspace(LowerFjordWall(i,1),LowerFjordWall(i+1,1),100)'];
    LfwY = [LfwY; linspace(LowerFjordWall(i,2),LowerFjordWall(i+1,2),100)'];
end


UpFW(:,1)=UfwX;
UpFW(:,2)=UfwY;
LoFW(:,1)=LfwX;
LoFW(:,2)=LfwY;

UpperFW = UpFW;
LowerFW = LoFW;

xmin = min(LowerFW(:,1));
ymin = min(LowerFW(:,2));

LoLine=[];
LoLine(1,:) = [xmin ymin];
LoLine= [LoLine;LowerFW];
LoLine(end,:) = [xmin ymin];
LowerPolygon =polyshape(LoLine);


% Upper boundary polygon
xmax = max(UpperFW(:,1));
ymax = max(UpperFW(:,2));

UpLine=[];
UpLine(1,:) = [xmax ymax];
UpLine= [UpLine;UpperFW];
UpLine(end,:) = [xmax ymax];
UpperPolygon =polyshape(UpLine);
plot(UpperPolygon)

BoundaryFolder=strcat(outFolder,'FjordBoundaries/');
if ~exist(BoundaryFolder, 'dir')
   mkdir(BoundaryFolder)
end


writematrix(UpperFW,strcat(BoundaryFolder,'UpperFW_',GlacierName,'.csv'))
writematrix(LowerFW,strcat(BoundaryFolder,'LowerFW_',GlacierName,'.csv'))
end

