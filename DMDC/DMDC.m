function [ components ] = DMDC(info,showComponents,createVideo,videoTitle)

validFirstInds = [0 info.segments(1:end-1)] + max(info.delays) + 1;
validFinalInds = info.segments;
info.y=info.q;

validInds = validFirstInds(1):validFinalInds(1);
components = double(info.ims(validFirstInds(1):validFinalInds(1),:))'*info.y(1:info.segments(1)-max(info.delays),:);

for i = 2:length(validFirstInds)
    validInds = [validInds validFirstInds(i):validFinalInds(i)];
    components = components+double(info.ims(validFirstInds(i):validFinalInds(i),:))'*info.y(1+info.segments(i-1)-max(info.delays)*(i-1):info.segments(i)-max(info.delays)*i,:);
end
%components = double(info.ims(validInds,:))'*info.y;

info.components=components;

if (info.micromodelNumber == 3)
    components = info.microModel.q*(info.macroVars(validInds,:))'*info.y;
    %components = (info.macroVars(firstIm:end,:))'*info.y;
    %components = info.microModel.q(:,info.microModel.microNinds)*info.microModel.b(info.microModel.microNinds,info.microModel.microNinds)*info.macroVars(firstIm:end,:)'*info.y;
    if (info.microModel.invFrame==0)
        info.microModel.invFrame = pinv(info.microModel.frame);
    end
    if (nargin > 1)
        if (showComponents == 1)
            figure(1);
            for i = 1:36
                subplot(6,6,i);
                imagesc(reshape(components(:,i)'*info.microModel.frame,[size(info.ims,2),size(info.ims,3)]));
                axis off;
            end
        end
    end
else
    if (nargin > 1)
        if (showComponents == 1)
            figure(1);
            for i = 1:9
                subplot(3,3,i);
                imagesc(reshape(components(:,i),[size(info.ims,2),size(info.ims,3)]));
                axis off;
            end
        end
    end
end



if (nargin >2)
    if (createVideo == 1)
        if (nargin < 4)
            videoTitle = 'infotest5.mpeg';
        end
        
        MakeQTMovie('start',videoTitle);
        MakeQTMovie('quality',.9);
        %MakeQTMovie('framerate',5);
        H = figure(2);
        set(H,'Position',[200 200 max(600,2*(info.width+100)) max(300,info.height+100)]);
        
        ninds = 1:info.numdims;
        imMin = min(info.ims(:));
        imMax = max(info.ims(:));
        
        for i = 1:1:(size(info.ims,1)-(max(length(info.segments),1))*max(info.delays))
            subplot(1,2,1);
            
            originalIndex = info.originalImageIndices(i)+max(info.delays);
            orig = double(squeeze(info.ims(originalIndex,:,:)));
            imagesc(orig,[imMin,imMax]);
            title('Original','FontSize',20);
            axis off; box on;
            
            subplot(1,2,2);
            reconstruction = reconstructImage(info,i,ninds);
            imagesc(reconstruction,[imMin,imMax]);
            title('DMDC Reconstruction','FontSize',20);
            axis off; box on;
            colormap(gray);
            drawnow;
            
            MakeQTMovie('addfigure');
            
        end
        
        MakeQTMovie('finish');
    end
end

end

