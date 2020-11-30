%5/8/18 Was probably written long before DB
%cycle through files and plot psychophysics collection data
%textFileList = {'EmilyLMTF.txt'};
startpath = 'C:/NO BACKUP/NexFiles/Emily/Emily';
for iterator = 1:length(flist)
    %fileName = textFileList{iterator};
    %flist = flatten(fnamesFromTxt2(fullfile(nexfilepath,'nexfilelists', 'Emily/Emily', fileName)));
    for f = 1:length(flist)
        filestro = nex2stro(findfile(flist{f},startpath));
        infopath = filestro.sum.exptParams;
        pixperdeg(f, :) = infopath.pixperdeg;
        framerate(f, :) = infopath.framerate;
        theta(f, :) = infopath.theta;
        sf(f, :) = infopath.sf;
        sigma(f, :) = infopath.sigma;
        nsigma(f, :) = infopath.nsigmas;
    end
    figure; 
    L = length(flist); x = 1:L;
    semilogy(x, framerate, x, pixperdeg, x, theta, x, sf, x, sigma, x, nsigma); hold on;
    ax = gca; ax.XTickLabel = flist;
    xlabel('nex name'); ylabel('various'); legend('framerate', 'pixperdeg', 'theta', 'sf', 'sigma', 'nsigma', 'Location', 'best');
    title(fileName);
end

%clear all between runs!