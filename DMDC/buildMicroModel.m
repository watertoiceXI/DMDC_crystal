function [macroVars,microModel] = buildMicroModel(info,micromodel,junk)

    if (micromodel == 1)     %NO MICRO STRUCTURE
  
%     s    keyboard
        [macroVars,microModel] = microNoStructure(info.ims,[],info.pr); 
   
    elseif (micromodel == 2) %SPATIAL MICRO STRUCTURE
        
        microXCoords = meshgrid(1:size(info.ims,2),1:size(info.ims,3));
        microYCoords = microXCoords';
        microCoords = [microXCoords(:) microYCoords(:)];
        [macroVars,microModel] = microSpatialStructure(info.ims(1:end,:),microCoords,64,[],256);
        
    elseif (micromodel == 3) %ADAPTED MICRO STRUCTURE
        
        %%% COMING SOON
        
        
    else
        
        [macroVars,microModel] = microNoStructure(info.ims,[],512); 

    end


end

