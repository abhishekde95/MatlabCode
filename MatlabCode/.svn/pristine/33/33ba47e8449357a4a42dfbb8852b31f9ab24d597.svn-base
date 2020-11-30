function [trialspec]=trialspecs(LMSdir,numReversals,thresh,stepSize,scale,phases)


stimIntensityIndx=[];
queue=[];
for z=1:phases
    if z==1
        for i=1:3
            trialspec(i).predcoord=LMSdir(i,:);
            trialspec(i).round=z;
            trialspec(i).inqueue=1;
        end
    end 
    if z==2
        mirror=sign((fullfact([1 2 2])-1.5)*-1);
        for i=1:4
            mirrorMatrix(:,:,i)=repmat(mirror(i,:)',1,3); %Creates matrixes = in size to coordinates (for multiplying)
            trialspec(length(trialspec)+1).parenttriangle=cat(1,trialspec(1:3).coordinates).*mirrorMatrix(:,:,i);
            trialspec(length(trialspec)).parentoog=cat(2,trialspec(1:3).oog);
            trialspec(length(trialspec)).predcoord=mean(trialspec(length(trialspec)).parenttriangle);
            trialspec(length(trialspec)).round=z;
            trialspec(length(trialspec)).inqueue=1;
        end
    end
    if z>2
        for i=queue
            if ~isempty(trialspec(i).coordinates) && abs(log(trialspec(i).coordinates/trialspec(i).predcoord))<0.3
                trialspec(i).active=0;
            elseif sum([trialspec(i).parentoog])==3
                trialspec(i).active=0;
            else
                trialspec(i).active=1;
                if ~isempty(trialspec(i).coordinates)
                    for n=1:3;a=[1;2;3];a(n,:)=[];
                        trialspec(length(trialspec)+1).parenttriangle=cat(1,trialspec(i).parenttriangle([a'],:),trialspec(i).coordinates);
                        trialspec(length(trialspec)).parentoog=cat(2,trialspec(i).parentoog([a']),trialspec(i).oog);
                        trialspec(length(trialspec)).predcoord=mean(trialspec(length(trialspec)).parenttriangle);
                        trialspec(length(trialspec)).round=z;
                        trialspec(length(trialspec)).inqueue=1;
                    end
                end
            end
        end
    end
    queue=[];
    for i=1:length(trialspec)
        if trialspec(i).inqueue==1
            queue=[queue i];
        end
    end
    if length(queue)>6
        queue=queue(floor([0:5]*(length(queue)/6))+1);
    end
    for n=1:length(queue);i=queue(n);
        [stimIntensityStair, oog] = stairCase2(trialspec(i).predcoord,numReversals,thresh,stepSize,scale);
        trialspec(i).oog=oog;
        if oog==0
            trialspec(i).stimintensity=stimIntensityStair(end);
        elseif isempty(stimIntensityStair) || length(stimIntensityStair)==1
            trialspec(i).stimintensity=1;
        else
            trialspec(i).stimintensity=stimIntensityStair(end-1);
        end
        trialspec(i).coordinates=(trialspec(i).predcoord).*trialspec(i).stimintensity;
        trialspec(i).vectorlength=norm(trialspec(i).coordinates);
        trialspec(i).inqueue=0;
    end
end