

makeMovie = true;


[TR, TT] = icp(subCloud1,subCloud2,'WorstRejection',.5);

plotsparsity = 1;
%subCloud1 = pointCloud1(:,pointStart:plotsparsity:pointEnd);
%subCloud2 = pointCloud2(:,pointStart:plotsparsity:pointEnd);
ff = figure(1);
set(ff,'Position',[0 0 2000 1000])
axis equal;
hold on
scatter3(subCloud1(2,1:sprsty:end),subCloud1(1,1:sprsty:end),-subCloud1(3,1:sprsty:end),3*ones(1,size(subCloud1(1,1:sprsty:end),2)),'c');
scatter3(measurements1(2,:),measurements1(1,:),-measurements1(3,:),'ko');
axis equal;
hold on
dispCloud2 = TR*subCloud2 + repmat(TT,1,size(subCloud2,2));
measCloud2 = TR*measurements2 + repmat(TT,1,size(measurements2,2));
scatter3(dispCloud2(2,1:sprsty:end),dispCloud2(1,1:sprsty:end),-dispCloud2(3,1:sprsty:end),3*ones(1,size(dispCloud2(1,1:sprsty:end),2)),'g^');
scatter3(measCloud2(2,:),measCloud2(1,:),-measCloud2(3,:),'r+');


if (makeMovie)
    movie_iter = 0;
    for iMovie = 0:360
        view(iMovie,20)
        drawnow()
        pause(0.01)
        movie_iter = movie_iter+1;
        Feature_Movie(movie_iter) = getframe(ff);
    end
    %Now, make the movie
    
    aviobj = avifile('NiceRegistration.avi','compression', 'None', 'fps', 24);
    for k=1:length(Feature_Movie)
        aviobj = addframe(aviobj,Feature_Movie(:,k));
    end
    aviobj = close(aviobj);
    
end