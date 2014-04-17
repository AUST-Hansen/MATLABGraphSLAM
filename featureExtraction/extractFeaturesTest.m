clear all; close all; clc;

load wallSegment.mat

% movie?
makeMovie = false;


DEBUG = false;

pointStart = 4000;
pointEnd = 10000; %round(length(pointCloud)/10);
%pointStart = 6000;
%pointEnd = 10000; %round(length(pointCloud)/10);
sparsity = 1;
neighborRadius = 5;
minNeighbors = 15*neighborRadius;

% scatter3(pointCloud(2,pointStart:pointEnd),pointCloud(1,pointStart:pointEnd),-pointCloud(3,pointStart:pointEnd),2*ones(1,size(pointCloud(1,pointStart:pointEnd))));
% hold on;
% scatter3(pointCloud(2,pointStart:sparsity:pointEnd),pointCloud(1,pointStart:sparsity:pointEnd),-pointCloud(3,pointStart:sparsity:pointEnd),2*ones(1,size(pointCloud(1,pointStart:sparsity:pointEnd))),'r');
    
% focus on smaller, subsampled cloud for initial testing
subCloud = pointCloud(:,pointStart:sparsity:pointEnd);

scatter3(subCloud(2,:),subCloud(1,:),-subCloud(3,:),2*ones(1,size(subCloud(1,:))),'r');
axis equal
view(150,10)

normals = zeros(size(subCloud));
curvatures = zeros(2,length(subCloud));

% for each point in 
neighborIndices = rangesearch(subCloud',subCloud',neighborRadius);
for ii = 1:length(subCloud)
    
    % weed out borders and spurious points
    if(length(neighborIndices{ii})<minNeighbors)
        continue;
    end
    
    % make neighbor matrix
    P = subCloud(:,neighborIndices{ii});
    Po = mean(P,2);
    Preg = P - repmat(Po,1,size(P,2));
    M = Preg*Preg';
    [u s v] = svd(M);
    normal = v(:,3);
    normals(:,ii) = normal;
    curvatures(1,ii) = s(3,3)/(s(3,3)+s(2,2));
    curvatures(2,ii) = s(3,3)/(s(3,3)+s(1,1));
    
    if(DEBUG)
        scatter3(Preg(2,:),Preg(1,:),-Preg(3,:),'r');
        hold on;
        axis equal
        quiver3(0,0,0,normal(2),normal(1),-normal(3));
        keyboard
    end
    
end
%%
ff = figure(1);
set(ff,'Position',[0 0 2000 1000])
scatter3(subCloud(2,:),subCloud(1,:),-subCloud(3,:),2*ones(1,size(subCloud(1,:))),'r');
hold on
axis equal
view(0,10)

minCurve = .06;
maxCurve = .2;
subsample = 1;
for ii = 1:subsample:length(subCloud)
    %if and(and(curvatures(1,ii) > minCurve, curvatures(1,ii) < maxCurve), and(curvatures(2,ii) > minCurve, curvatures(2,ii) < maxCurve) ) 
    if and(and(curvatures(1,ii) > minCurve, true), and(curvatures(2,ii) > minCurve, true ))   
        scatter3(subCloud(2,ii),subCloud(1,ii),-subCloud(3,ii),'b');
    end
end


if (makeMovie)
    movie_iter = 0;
    for iMovie = 0:360
        view(iMovie,10)
        drawnow()
        pause(0.01)
        movie_iter = movie_iter+1;
        Feature_Movie(movie_iter) = getframe(ff);
    end
    %Now, make the movie
    
    aviobj = avifile('CloudFeatureBig.avi','compression', 'None', 'fps', 24);
    for k=1:length(Feature_Movie)
        aviobj = addframe(aviobj,Feature_Movie(:,k));
    end
    aviobj = close(aviobj);


    
end
view(160,10)
%%
figure;
plot(curvatures')
legend('curvature1','curvature2')
title('curvature')