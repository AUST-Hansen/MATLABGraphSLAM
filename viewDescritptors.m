
for ii = 1:size(featureDescriptors,2)
    figure(2)
    plot(featureDescriptors(:,ii));
    drawnow()
    pause(.01)
end