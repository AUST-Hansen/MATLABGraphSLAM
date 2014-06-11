function [normals curvatures] = extractCloudNormalStats(cloud)

    subCloud = cloud(:,1:50:end);
    fprintf('calculating normals and curvature...\n');
    neighborIndices = rangesearch(cloud',subCloud',5);
    normals = zeros(size(subCloud));
    curvatures = -ones(1,size(subCloud,2));
    for ii = 1:length(subCloud)
        
        % weed out borders and spurious points
        if(length(neighborIndices{ii}) < 10)
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
        %curvatures(1,ii) = s(3,3)/(s(3,3)+s(2,2));
        %curvatures(2,ii) = s(3,3)/(s(3,3)+s(1,1));
        curvatures(ii) = s(3,3)/(s(3,3)+s(2,2)+s(1,1)); %[0 0 1]*v(:,3)*sign([1 0 0]*v(:,3));

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