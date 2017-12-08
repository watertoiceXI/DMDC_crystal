function [ macroVars, microModel ] = microAdaptedSpatialStructure( macroData, microCoords, microk, microModel, macrodim )
        % macroData is T x N1 x N2 x ... x Nk where T is time and Ni are
        % spatial dimensions and N = N1*N2*...*Nk is the number of spatial
        % nodes.
        
        % microCoords is N x L dimensional giving the L-dimensional
        % coordinates of each spatial node.  For example, in a 2D image,
        % the macroData is size T x W x H where width and height are the
        % number of horizontal and vertical pixels, N = W*H is the total
        % number of pixels, and microCoords is N x 2 and contains the
        % location of each pixel within the image.
        
        fprintf('Adapted Microscopic Structure...      ');    
        tic;
        
        if (isempty(microModel))
            spatialdim=size(macroData);
            spatialdim=prod(spatialdim(2:end));
            microModel.spatialdim = spatialdim;
            microModel.macrodim = macrodim;
           [microModel.y,microModel.q,microModel.b,microModel.d,microModel.sigma] = DiffusionMapping(microCoords,microk,macrodim,0,1,0); 
            
        end

        macroVars = double(reshape(macroData,[size(macroData,1),microModel.spatialdim]))*microModel.y;  

        toc;
end

