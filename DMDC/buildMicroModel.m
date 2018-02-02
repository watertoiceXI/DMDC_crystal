function [macroVars,microModel] = buildMicroModel(info,micromodel)

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
        
        
    elseif (micromodel == 4) %LDA MICRO STRUCTURE
        sz = size(info.ims);
        aux_sz = size(info.aux);
        spatialdim = prod(sz(2:end));
        density = info.sparse_average/spatialdim;
        microModel.spatialdim = spatialdim;
        microModel.projection = zeros(spatialdim,aux_sz(2));
        macroVars = double(info.ims(:,:));
        for q = 1:info.pr
            if density == 0
                rnd_proj_basis = randn(spatialdim,info.pr_prime);
                [rnd_proj_basis,~]=qr(rnd_proj_basis,0);
            else
                rnd_proj_basis = sprand(spatialdim,info.pr,density);
                for w = 1:info.pr
                    %verify that there's something in the basis. LDA hates it
                    %when there's not. 
                    while ~any(rnd_proj_basis(:,w))
                        rnd_proj_basis(:,w)=sprand(spatialdim,1,density);
                    end
                end
            end
            rnd_proj = info.ims(:,:)*rnd_proj_basis;
            W = LDA(rnd_proj,info.aux(:,q)>0);
            microModel.projection(:,q) = W(1,1)+rnd_proj_basis*W(1,2:end)';
        end
    else
        
        [macroVars,microModel] = microNoStructure(info.ims,[],512); 

    end


end

