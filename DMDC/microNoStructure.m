function [ macroVars, microModel ] = microNoStructure( macroData, microModel, macrodim )
% macroData is T x N1 x N2 x ... x Nk where T is time and Ni are
% spatial dimensions

fprintf('No Microscopic Structure...           ');
tic;
% % %
%          if (isempty(microModel))
%             spatialdim=size(macroData);
%             spatialdim=prod(spatialdim(2:end));
%             microModel.spatialdim = spatialdim;
%             microModel.macrodim = spatialdim;%macrodim;
%             microModel.projection = randn(spatialdim,macrodim);
%             [microModel.projection,~] = qr(microModel.projection,0);
%             microModel.projection = eye(spatialdim);
%         end
%         macroVars = double(reshape(macroData,[size(macroData,1),microModel.spatialdim]))*microModel.projection;
% %         min(min(macroVars))
% %         max(max(macroVars))
% %         keyboard
% 
%         toc;
% end

% %
% % macrodim = 20;
% % microModel.projection = randn(spatialdim,macrodim);
% % [microModel.projection,~] = qr(microModel.projection,0);
% % if (spatialdim<macrodim)
% %      microModel.projection = randn(macrodim,spatialdim);
% %      [microModel.projection,~] = qr(microModel.projection,0);
% %      microModel.projection = microModel.projection';
% % end
% %
%   macroVars = double(reshape(macroData,[size(macroData,1),microModel.spatialdim]))*microModel.projection;


% % %%%%Ty2
% macrodim = 1000;
 macrodim=100;
 
spatialdim=size(macroData);
spatialdim=prod(spatialdim(2:end));
microModel.spatialdim = spatialdim;
microModel.projection = randn(spatialdim,macrodim);
[microModel.projection,~] = qr(microModel.projection,0);
if (spatialdim<macrodim)
    microModel.projection = randn(macrodim,spatialdim);
    [microModel.projection,~] = qr(microModel.projection,0);
    microModel.projection = microModel.projection';
end

if size(macroData,2)>10
    1
 macroVars = double(reshape(macroData,[size(macroData,1),microModel.spatialdim]))*microModel.projection;
else
    macroVars = double(reshape(macroData,[size(macroData,1),microModel.spatialdim]));%*microModel.projection;
end
 clear macroData
macroVars = 0.5*macroVars./repmat(mean(macroVars.^2,1).^(1/2),size(macroVars,1),1);

toc;
