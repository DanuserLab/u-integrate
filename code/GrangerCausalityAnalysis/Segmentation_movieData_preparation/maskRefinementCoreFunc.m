function currMask = maskRefinementCoreFunc(currMask, p)
% maskRefinementCoreFunc is a part of refineMovieMasks.m to implement
% imclose, object number and filling holes.
% 2017/05/29

        % ----- Mask Clean Up ------ %
        
seClose = strel('disk',p.ClosureRadius,0);

            %Remove objects that are too small
            if p.MinimumSize > 0
                currMask = bwareaopen(currMask,p.MinimumSize);
            end            
            
            %Perform initial closure operation
            if p.ClosureRadius > 0
                currMask = imclose(currMask,seClose);            
            end            
            
            %%Perform initial opening operation
            %if p.OpeningRadius > 0
            %    currMask = imopen(currMask,seOpen);
            %end
            
        
        
        % ---------- Object Selection -------- %
        
        %Keep only the largest objects
        if ~isinf(p.ObjectNumber)
                
            %Label all objects in the mask
            labelMask = bwlabel(currMask);

            %Get their area
            obAreas = regionprops(labelMask,'Area');       %#ok<MRPBW>

            %First, check that there are objects to remove
            if length(obAreas) > p.ObjectNumber 
                obAreas = [obAreas.Area];
                %Sort by area
                [dummy,iSort] = sort(obAreas,'descend'); %#ok<ASGLU>
                %Keep only the largest requested number
                currMask = false(size(currMask));
                for i = 1:p.ObjectNumber
                    currMask = currMask | labelMask == iSort(i);
                end
            end
        end
        
        % ------ Hole-Filling ----- %
        if p.FillHoles
            
            %If the mask touches the image border, we want to close holes
            %which are on the image border. We do this by adding a border
            %of ones on the sides where the mask touches.
            if any([currMask(1,:) currMask(end,:) currMask(:,1)' currMask(:,end)'])                
                m = size(currMask, 1);
                n = size(currMask, 2);            
                %Add a border of 1s
                tmpIm = vertcat(true(1,n+2),[true(m,1) ...
                                currMask true(m,1)],true(1,n+2));
                
                %Find holes - the largest "hole" is considered to be the
                %background and ignored.
                cc = bwconncomp(~tmpIm,4);                                
                holeAreas = cellfun(@(x)(numel(x)),cc.PixelIdxList);
                [~,iBiggest] = max(holeAreas);                                
                tmpIm = imfill(tmpIm,'holes');
                tmpIm(cc.PixelIdxList{iBiggest}) = false;
                currMask = tmpIm(2:end-1,2:end-1);
             else                        
                 currMask = imfill(currMask,'holes');
            end
        end
        
end

        