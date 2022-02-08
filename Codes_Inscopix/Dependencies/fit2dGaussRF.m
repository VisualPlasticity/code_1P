function [thisparams,thisgaussian,bestbeta,coeff,pvalcorr] = fit2dGaussRF(EvokedResponse)

% Create dataset with possible Gaussians
nX = size(EvokedResponse,1);
nY = size(EvokedResponse,2);
nCells = size(EvokedResponse,3);
Centroids = combvec(1:nX,1:nY);
RFSizes = [0.2:0.2:max([nX,nY])/3]; %RF sizes in Full-width half-maximum

ngauss = length(Centroids)*length(RFSizes);
% stimmask = zeros(nX*nY,nX*nY);
% for i = 1:nX*nY
%     stimmask(i,i) = 1;
% end

%Make all Gaussians
disp(['Creating ' num2str(ngauss) ' Gaussians...'])
Gaussians = nan(nX,nY,ngauss);
Gausparams = nan(3,ngauss); %xcenter,Ycenter,FWHM
nrfsizes = length(RFSizes);
countid = 1;
for centerid = 1:size(Centroids,2)
    tmpgauss = nan(nX,nY,nrfsizes);
    tmpparam = nan(3,nrfsizes);
    parfor szid = 1:nrfsizes
        gw = normpdf(1:nX,Centroids(1,centerid),RFSizes(szid));
        gv = normpdf(1:nY,Centroids(2,centerid),RFSizes(szid));
        Z = gv'*gw;
        Z = Z./sum(sum(Z)); %Could normalise by sum
        tmpgauss(:,:,szid) = Z';
        tmpparam(:,szid) = [Centroids(:,centerid); RFSizes(szid)];
    end
    Gaussians(:,:,countid:countid+nrfsizes-1) = tmpgauss;
    Gausparams(:,countid:countid+nrfsizes-1) = tmpparam;
    countid = countid+nrfsizes;
end

%% Find Gaussian with largest weight
disp(['Finding best Gaussian for all ' num2str(nCells) ' units...'])
[b,std,mse] = lscov(reshape(Gaussians,[],ngauss),reshape(EvokedResponse,nX*nY,[]));

[maxval,maxid] = nanmax(b,[],1);

%Save parameters
disp('Saving out parameters...')
thisparams = Gausparams(:,maxid);
thisparams(3,:) = thisparams(3,:); %convert to fwhm
thisgaussian = Gaussians(:,:,maxid);
bestbeta = cell2mat(arrayfun(@(X) b(X),maxid,'UniformOutput',0));
%statistics
ypred = arrayfun(@(X) Gaussians(:,:,X).*b(X),maxid,'UniformOutput',0);
ypred = cat(3,ypred{:});
ypred = reshape(ypred,nX*nY,nCells);
EvokedResponse = reshape(ypred,nX*nY,nCells);
[coeff,pvalcorr] = arrayfun(@(X) corr(ypred(:,X),EvokedResponse(:,X)),1:nCells,'UniformOutput',0);
coeff = cell2mat(coeff);
pvalcorr = cell2mat(pvalcorr);

end