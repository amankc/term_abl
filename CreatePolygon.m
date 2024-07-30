function [UpperPoly,LowerPoly] = CreatePolygon(FW,FWl)
    xmin = min(FWl(:,1));
    ymin = min(FWl(:,2));
    
    LoLine=[];
    LoLine(1,:) = [xmin ymin];
    LoLine= [LoLine;FWl];
    LoLine(end,:) = [xmin ymin];
    LowerPoly =polyshape(LoLine);
    
    
    % Upper boundary polygon
    xmax = max(FW(:,1));
    ymax = max(FW(:,2));
    
    UpLine=[];
    UpLine(1,:) = [xmax ymax];
    UpLine= [UpLine;FW];
    UpLine(end,:) = [xmax ymax];

    UpperPoly =polyshape(UpLine);
end