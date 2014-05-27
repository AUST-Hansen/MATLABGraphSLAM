function [measurements_out, descriptors_out] = extractFeaturesFromSubmap(cloud,varargin)

    % cloud must be expressed in the frame of reference of the vehicle at
    % the reference time, with the origin of the cloud at the vehicle
    % location.
    
%% Use the inputParser class to validate input arguments.
inp = inputParser;

inp.addRequired('cloud', @(x)isreal(x) && size(x,1) == 3);

inp.addOptional('minRadius', 4.0, @(x)isscalar(x) && x > 0);

inp.addOptional('maxRadius', 4.0, @(x)isscalar(x) && x > 0);

inp.addOptional('numRadii', 4, @(x)isscalar(x) && x > 0);

inp.addOptional('MinNeighbors', 5, @(x)isscalar(x) && x > 0); % multiplier on radius

inp.addOptional('Sparsity', 1, @(x)isscalar(x) && x > 0);

inp.addOptional('MinCurve', .1, @(x)isscalar(x) && x > 0);

inp.addOptional('MaxCurve', .3, @(x)isscalar(x) && x > 0);

inp.addOptional('PFHbins', 5,@(x) x > 0);

inp.addOptional('VertMargin',.2,@(x)isscalar(x) && x > 0 && x<.5);

inp.addOptional('HorizMargin',.2,@(x)isscalar(x) && x > 0 && x<.5);

%validMinimize = {'point','plane','lmapoint'};
%inp.addParamValue('Minimize', 'point', @(x)any(strcmpi(x,validMinimize)));
%inp.addParamValue('Triangulation', [], @(x)isreal(x) && size(x,2) == 3);

inp.addParamValue('Verbose', false, @(x)islogical(x));


inp.parse(cloud,varargin{:});
arg = inp.Results;
clear('inp');

if (arg.minRadius == arg.maxRadius)
    Radius = arg.minRadius;
else
    Radius = linspace(arg.minRadius,arg.maxRadius,arg.numRadii);
end

multi_level = cell(length(Radius),1);
% calculate margins
maxZ = max(cloud(3,:));
minZ = min(cloud(3,:));
maxX = max(cloud(1,:));
minX = min(cloud(1,:));
Zrange = maxZ-minZ;
Xrange = maxX-minX;
topZ = maxZ - arg.VertMargin*Zrange;
botZ = minZ + arg.VertMargin*Zrange;
topX = maxX - arg.HorizMargin*Xrange;
botX = minX + arg.HorizMargin*Xrange;
subSampledCloud = cloud(:,cloud(3,:)<topZ & cloud(3,:)>botZ & cloud(1,:)<topX & cloud(1,:)>botX);
subSampledCloud = subSampledCloud(:,1:arg.Sparsity:end);
normals = zeros(size(cloud));

  %% Precalculate all normals  
    fprintf('calculating normals...\n');
    neighborIndices = rangesearch(cloud',cloud',Radius(end));
    for ii = 1:length(cloud)
        
        % weed out borders and spurious points
        if(length(neighborIndices{ii})<(arg.MinNeighbors))
            continue;
        end
        
        % make neighbor matrix
        P = cloud(:,neighborIndices{ii});
        Po = mean(P,2);
        Preg = P - repmat(Po,1,size(P,2));
        M = Preg*Preg';
        [u s v] = svd(M);
        normal = v(:,3);
        % Make sure normal is pointed back at vehicle. This check will only get weird
        % with really funky geometries, and should never really cause a
        % problem, based on how the cloud is defined (in vehicle ref frame).
        if (normal'*Po > 0)
            normal = -normal;
        end
        normals(:,ii) = normal;
        curvatures(1,ii) = s(3,3)/(s(3,3)+s(2,2));
        curvatures(2,ii) = s(3,3)/(s(3,3)+s(1,1));
        curvatures(3,ii) = [0 0 1]*v(:,3)*sign([1 0 0]*v(:,3));
        
        if(false && arg.Verbose)
            scatter3(Preg(2,:),Preg(1,:),-Preg(3,:),'r');
            axis equal
            hold on;
            quiver3(0,0,0,normal(2),normal(1),-normal(3));
            drawnow();
            pause(.01)
            hold off
        end
        
    end

    %% Cache PFH info
    fprintf('caching pfh info...\n');
sparseIdx = 1;
cacheflag = logical(zeros(length(cloud)));
cache = zeros(length(subSampledCloud),length(subSampledCloud),3);
neighborIndices = rangesearch(cloud',subSampledCloud',Radius(end));    
tic
    for ii = 1:size(subSampledCloud,2)
        
        if (mod(ii,15) == 0)
            fprintf('%d of %d points\n',ii,length(subSampledCloud))
        end
        % weed out borders and spurious points
        try
            if(length(neighborIndices{ii})<arg.MinNeighbors)
                continue;
            end
        catch
            keyboard
        end
        %if and(and(subCurvatures(1,ii) > arg.MinCurve, subCurvatures(1,ii) < arg.MaxCurve), and(subCurvatures(2,ii) > arg.MinCurve, subCurvatures(2,ii) < arg.MaxCurve ))
        % good keypoint. Now describe it
        P = cloud(:,neighborIndices{ii});
        patchNormals = normals(:,neighborIndices{ii});
        descriptor = zeros(arg.PFHbins^3,1);
        for jj = 1:length(P)-1
            for kk = jj+1:length(P)
                
                if(~cacheflag(neighborIndices{ii}(jj),neighborIndices{ii}(kk)) )
                    n1 = patchNormals(:,jj);
                    n2 = patchNormals(:,kk);
                    p1 = P(:,jj);
                    p2 = P(:,kk);
                    % start calculating point feature histogram
                    unitDP = (p2-p1)/norm(p2-p1);
                    u = n1;
                    v = cross(u,unitDP);
                    w = cross(u,v);
                    alpha = acos(v'*n2); % -1:1
                    phi = acos(u'*unitDP); % -1:1
                    theta = atan2(w'*n2,u'*n2); %-pi:pi
                    % store values
                    cache(neighborIndices{ii}(jj),neighborIndices{ii}(kk),1) = alpha;
                    cache(neighborIndices{ii}(jj),neighborIndices{ii}(kk),2) = phi;
                    cache(neighborIndices{ii}(jj),neighborIndices{ii}(kk),3) = theta;
                    % set flag so we don't compute again
                    cacheflag(neighborIndices{ii}(jj),neighborIndices{ii}(kk)) = 1;
                else
                    continue
                end
            end
        end
        %end
    end
    toc

    
    
for i_rad = 1:length(Radius)
    fprintf('r = %d\n',Radius(i_rad))
    % for each point in
    measurements = zeros(size(cloud));
    measindices = zeros(1,size(cloud,2));
    descriptors = zeros(arg.PFHbins^3,size(cloud,2));
    measurementsCounter = 0;
    %subCurvatures = curvatures(:,1:arg.Sparsity:end);
    
    
        fprintf('calculating patches from subsampled cloud...\n');
    for ii = 1:size(subSampledCloud,2)
        
        if (mod(ii,15) == 0)
            fprintf('%d of %d points\n',ii,length(subSampledCloud))
        end
        % weed out borders and spurious points
        if(length(neighborIndices{ii})<arg.MinNeighbors)
            continue;
        end
        %if and(and(subCurvatures(1,ii) > arg.MinCurve, subCurvatures(1,ii) < arg.MaxCurve), and(subCurvatures(2,ii) > arg.MinCurve, subCurvatures(2,ii) < arg.MaxCurve ))
        % good keypoint. Now describe it
        P = cloud(:,neighborIndices{ii});
        patchNormals = normals(:,neighborIndices{ii});
        descriptor = zeros(arg.PFHbins^3,1);
        for jj = 1:length(P)-1
            for kk = jj+1:length(P)
                alpha = cache(neighborIndices{ii}(jj),neighborIndices{ii}(kk),1);
                phi = cache(neighborIndices{ii}(jj),neighborIndices{ii}(kk),2); % -1:1
                theta = cache(neighborIndices{ii}(jj),neighborIndices{ii}(kk),3); %-pi:pi
                % calculate indices into subhistograms
                ialpha = int8(floor((alpha)/(pi)*arg.PFHbins ));
                iphi = int8(floor((phi)/(pi)*arg.PFHbins ));
                itheta = int8(floor((theta+pi)/(2*pi)*arg.PFHbins ));
                descriptorIndex = arg.PFHbins^2*(ialpha) + arg.PFHbins*(iphi) + itheta + 1;
                % add it to the histogram
                descriptor(descriptorIndex) = descriptor(descriptorIndex) + 1;
            end
        end
        if (arg.Verbose)
            plot(descriptor);
            title(num2str(ii));
            drawnow()
            length(P)
            pause(.1)
        end
        % record the descriptor and measurements
        measurementsCounter = measurementsCounter + 1;
        measurements(:,measurementsCounter) = subSampledCloud(:,ii);
        measindices(measurementsCounter) = ii;
        descriptors(:,measurementsCounter) = descriptor/norm(descriptor);
        %end
        
    end
    
    
    % trim
    measurements = measurements(:,1:measurementsCounter);
    measindices = measindices(1:measurementsCounter);
    descriptors = descriptors(:,1:measurementsCounter);
    
    muHist = mean(descriptors,2);
    muHists = repmat(muHist,1,measurementsCounter);
    % Kullback-Leibler distance (divergence)
    fprintf('calculating divergence from mu-histogram...\n');
    KLdivergence = sum( (descriptors -muHists).*log((descriptors+eps)./(muHists +eps)));
    
    % use distance metric to pull out features greater than one std away from
    % the mean Kullback-Leibler distance (divergence) as in "Aligning Point
    % Cloud Views using Persistent Feature Histograms" by Radu Rusu, et.al.
    fprintf('calculating divergence from mu-histogram...\n');
    goodcorners = abs(KLdivergence -mean(KLdivergence)) > 5*std(KLdivergence);
    fprintf('extracting feature candidates...\n');
    measurements = measurements(:,goodcorners);
    measindices = measindices(goodcorners);
    descriptors = descriptors(:,goodcorners);
    
    multi_level{i_rad}.measurements = measurements;
    multi_level{i_rad}.measindices = measindices;
    multi_level{i_rad}.descriptors = descriptors;

end % for i_rad...

% look for persistent features
if length(multi_level)>1
    fprintf('looking for persistent features...\n');
    measurements_out = zeros(3,100);
    descriptors_out = zeros(size(descriptors,1),100);
    persistCount = 1;
    for ii = 1:length(subSampledCloud)
        for jj = 1:length(multi_level)-1
            idx1 = find(multi_level{jj}.measindices==ii);
            idx2 = find(multi_level{jj+1}.measindices==ii);
            if and(~isempty(idx1), ~isempty(idx2))
                % persistence!
                measurements_out(:,persistCount) = multi_level{jj}.measurements(:,idx1);
                descriptors_out(:,persistCount) = multi_level{jj}.descriptors(:,idx1);
                persistCount = persistCount + 1;
                break;
            end
        end
    end
    measurements_out = measurements_out(:,1:persistCount-1);
    descriptors_out = descriptors_out(:,1:persistCount-1);
    persistCount
else
    measurements_out = measurements;
    descriptors_out = descriptors;
end
    if true % DEBUG
        figure(3)
        hold off
        keyboard
        scatter3(cloud(2,:),cloud(1,:),-cloud(3,:),ones(1,size(cloud,2)),curvatures(1,:));hold on;scatter3(subSampledCloud(2,:),subSampledCloud(1,:),-subSampledCloud(3,:),'r+');scatter3(measurements_out(2,:),measurements_out(1,:),-measurements_out(3,:),'g^')
        view(90,0)
        axis equal
        if(~isempty(descriptors_out))
            figure(4)
            plot(descriptors_out - repmat(mean(descriptors_out,2),1,size(descriptors_out,2)))
        end
        keyboard
    end

