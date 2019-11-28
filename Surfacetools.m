% Tools for operating on 3D surfaces
% Works around the surfstruce class
% surfstruct fields:
%         vertices = 0;
%         triangles = 0;
%         dim = 0;
%         vox = 0;
%         orient = 0;
%         color = 0;
%         id = 0;
%         numverts = 0;
%         numtris = 0;

classdef Surfacetools
    
    methods (Static = true)
        
        function vertexNormals  = computeVertexNormals(surfin)
            % numvertsx3 normal vectors of each vertex
            
            if( surfin.isMm == 0)
                surfin = Surfacetools.convertToMm(surfin);
            end
            vertexNormals = zeros(surfin.numverts, 3);
            surfin.vertexNtriangles = Surfacetools.getVertexNtriangles(surfin);
            surfin.triangleNormals  = Surfacetools.computeTriangleNormals(surfin);
            
            %find normals of all points
            for ii = 1 : surfin.numverts
                polys = surfin.vertexNtriangles{ii}; % get neighboring polygons to this point
                % get neighboring polygons to all these polygons
                %for jj = 1 : length(polys)
                %   polys = [polys mesh.triangleNtriangles{polys(jj)}];
                %end
                %polys = unique(polys); % remove redundancy
                normal2 = zeros(1,3);
                for jj = 1 : length(polys)
                    tmp = [  surfin.vertices(surfin.triangles(polys(jj),1),:); ...
                        surfin.vertices(surfin.triangles(polys(jj),2),:); ...
                        surfin.vertices(surfin.triangles(polys(jj),3),:) ...
                        ];
                    v1 = tmp(1,:) - tmp(2,:);
                    v2 = tmp(1,:) - tmp(3,:);
                    
                    normal2 = normal2 + surfin.triangleNormals(polys(jj),:)*0.5*norm(cross(v1,v2));
                    
                end
                if( norm(normal2) > 0.00001)
                    normal2 = normal2 / norm(normal2);
                end
                
                vertexNormals(ii,:) = normal2;
            end
        end
        
        function surfout  = computeVertexNormals2(surfin)
            % numvertsx3 normal vectors of each vertex
            
            flag = 0;
            if( surfin.isMm == 0)
                surfin = Surfacetools.convertToMm(surfin);
                flag = 1;
            end
            vertexNormals = zeros(surfin.numverts, 3);
            surfin.vertexNtriangles = Surfacetools.getVertexNtriangles(surfin);
            surfin.triangleNormals  = Surfacetools.computeTriangleNormals(surfin);
            
            %find normals of all points
            for ii = 1 : surfin.numverts
                polys = surfin.vertexNtriangles{ii}; % get neighboring polygons to this point
                % get neighboring polygons to all these polygons
                %for jj = 1 : length(polys)
                %   polys = [polys mesh.triangleNtriangles{polys(jj)}];
                %end
                %polys = unique(polys); % remove redundancy
                normal2 = zeros(1,3);
                for jj = 1 : length(polys)
                    tmp = [  surfin.vertices(surfin.triangles(polys(jj),1),:); ...
                        surfin.vertices(surfin.triangles(polys(jj),2),:); ...
                        surfin.vertices(surfin.triangles(polys(jj),3),:) ...
                        ];
                    v1 = tmp(1,:) - tmp(2,:);
                    v2 = tmp(1,:) - tmp(3,:);
                    
                    normal2 = normal2 + surfin.triangleNormals(polys(jj),:)*0.5*norm(cross(v1,v2));
                    
                end
                if( norm(normal2) > 0.00001)
                    normal2 = normal2 / norm(normal2);
                end
                
                vertexNormals(ii,:) = normal2;
            end
            surfin.vertexNormals = vertexNormals;
            
            surfout = surfin;
            if(flag)
                surfout = Surfacetools.convertToVox(surfout);
            end
                
        end
        
        function triangleNormals  = computeTriangleNormals(surfin)
            % numtrisx3 normal vectors of each triangle
            triangleNormals = zeros(surfin.numtris, 3);
            for ii = 1:surfin.numtris %find normals of all polygons
                %indices of the points from which the polygon is made
                pointIndex1 = surfin.triangles(ii,1);
                pointIndex2 = surfin.triangles(ii,2);
                pointIndex3 = surfin.triangles(ii,3);
                
                %coordinates of the points
                point1 = surfin.vertices(pointIndex1,:);
                point2 = surfin.vertices(pointIndex2,:);
                point3 = surfin.vertices(pointIndex3,:);
                
                vector1 = point2 - point1;
                vector2 = point3 - point2;
                normal = cross(vector1,vector2);
                if( norm(normal) > 0.00001)
                    normal = normal / norm(normal);
                end
                
                triangleNormals(ii,:)=normal;
            end
        end
        
        function vertexNtriangles  = getVertexNtriangles(surfin)
            % numvertsx1 cell of neighboring triangles of each vertex
            vertexNtriangles = cell(surfin.numverts, 3);
            
            for ii = 1 : surfin.numtris
                %indices of the points from which the polygon is made
                pointIndex1 = surfin.triangles(ii,1);
                pointIndex2 = surfin.triangles(ii,2);
                pointIndex3 = surfin.triangles(ii,3);
                
                %make entry of this polygon as the neighbouring polygon of the three
                %vertex points
                vertexNtriangles(pointIndex1,1)={[vertexNtriangles{pointIndex1,1} ii]};
                vertexNtriangles(pointIndex2,1)={[vertexNtriangles{pointIndex2,1} ii]};
                vertexNtriangles(pointIndex3,1)={[vertexNtriangles{pointIndex3,1} ii]};
            end
            
        end
        
        function triangleNtriangles  = getTriangleNtriangles(surfin)
            % numtrisx1 cell of nieghboring triangles of each triangle
            triangleNtriangles = cell(surfin.numtris, 1);
            
            % find neighbouring polygons of all polygons
            for ii = 1 : surfin.numtris
                polNeighbor = [];
                for jj = 1 : 3
                    polNeighbor = [polNeighbor surfin.vertexNtriangles{surfin.triangles(ii,jj)}];
                end
                polNeighbor = unique(polNeighbor);
                polNeighbor = setdiff(polNeighbor, [ii]);
                triangleNtriangles(ii,1)={[polNeighbor]};
            end
        end
        
        function oneRingNvertices  = getOneRingNvertices(surfin)
            % numvertsx1 cell of nieghboring one ring vertices of each
            % vertex
            oneRingNvertices = cell(surfin.numverts, 1);
            if (isempty(surfin.vertexNtriangles) || ~iscell(surfin.vertexNtriangles))
                surfin.vertexNtriangles = Surfacetools.getVertexNtriangles(surfin);
            end
            
            % find neighbouring vertices 
            for ii = 1 : surfin.numverts
                pol = surfin.vertexNtriangles{ii, 1};
                pts = [];
                for jj = 1:length(pol)
                    pointIndex1 = surfin.triangles(pol(jj),1);
                    pointIndex2 = surfin.triangles(pol(jj),2);
                    pointIndex3 = surfin.triangles(pol(jj),3);
                    pts = [pts pointIndex1];
                    pts = [pts pointIndex2];
                    pts = [pts pointIndex3];
                end
                pts = unique(pts);
                pts = setdiff(pts, [ii]);
                oneRingNvertices{ii} = pts;
            end
                
        end
        
        function twoRingNvertices  = getTwoRingNvertices(surfin)
            % numvertsx1 cell of nieghboring two ring vertices of each
            % vertex
            twoRingNvertices = cell(surfin.numverts, 1);
            if (isempty(surfin.vertexNtriangles))
                surfin.vertexNtriangles = Surfacetools.getVertexNtriangles(surfin);
            end
            triangleNtriangles  = Surfacetools.getTriangleNtriangles(surfin);
            oneRingNvertices = Surfacetools.getOneRingNvertices(surfin);
            % find neighbouring polygons of all polygons
            for ii = 1 : surfin.numverts
                pol = surfin.vertexNtriangles{ii, 1};
                polNeighbor = [];
                for jj = 1:length(pol)
                    polNeighbor = [polNeighbor, triangleNtriangles{pol(jj)}];
                end
                polNeighbor = unique(polNeighbor);
                pts = [];
                for jj = 1:length(polNeighbor)
                    pointIndex1 = surfin.triangles(polNeighbor(jj),1);
                    pointIndex2 = surfin.triangles(polNeighbor(jj),2);
                    pointIndex3 = surfin.triangles(polNeighbor(jj),3);
                    pts = [pts pointIndex1];
                    pts = [pts pointIndex2];
                    pts = [pts pointIndex3];
                end
                pts = unique(pts);
                pts = setdiff(pts, oneRingNvertices{ii});
                pts = setdiff(pts, [ii]);
                twoRingNvertices{ii} = pts;
            end
        end
        
        function nRingNvertices = getNRingNvertices(surfin, n)
            nRingNvertices = cell(surfin.numverts, 1);
%             oneRingNvertices = Surfacetools.getOneRingNvertices(surfin);
%             if n == 1
%                 nRingNvertices = oneRingNvertices;
%                 return;
%             end
%             twoRingNvertices  = Surfacetools.getTwoRingNvertices(surfin);
%             if n == 2
%                 nRingNvertices = twoRingNvertices;
%                 return;
%             end
%             iter = floor(n / 2) - 1;
%             tt = mod(n, 2);
%             for i = 1:iter
%                 oneRingNvertices = 
%             end
%             if tt~=0
%             end
            
        end
        
        function [surf1, surf2] = preventOverlap(surf1, surf2, stepsize)
        end
        
        function surf = readSurfAsMesh(filename)
            % read surface as mesh from a file
            fid = fopen(filename,'rb');
            id = fread(fid,1,'int');
            
            numverts = fread(fid,1,'int');
            numtris = fread(fid,1,'int');
            c1 = fread(fid, 1,'int');
            color = zeros(3,1);
            if (c1 ~= -1)
                color(1,1) = c1;
                color(2:3,1) = fread(fid,2,'int');
                dim = [352 352 350];
                vox = [0.7 0.7 0.7];
                orient = [1 2 3];
            else
                orient = fread(fid, 3, 'int');
                dim = fread(fid,3,'int');
                vox = fread(fid,3,'float');
                color = fread(fid,3, 'int');
            end
            
            vertices = fread(fid,numverts*3,'float');
            triangles = fread(fid,numtris*3,'int');
            fclose(fid);
            
            vertices = reshape(vertices, [3,numverts])';
            triangles = reshape(triangles, [3, numtris])';
            triangles = triangles + ones(size(triangles));
            
            surf = surfstruct;
            
            surf.id = id;
            surf.numverts = length(vertices);
            surf.numtris = length(triangles);
            surf.vertices = vertices;
            surf.triangles = triangles;
            surf.dim = dim;
            surf.vox = vox;
            surf.orient = orient;
            surf.color = color;
        end
        
        function surf = readSurfAsRMesh(filename)
            % read surface as mesh from a file
            fid = fopen(filename,'rb');
            id = fread(fid,1,'int');
            
            numverts = fread(fid,1,'int');
            numtris = fread(fid,1,'int');
            c1 = fread(fid, 1,'int');
            color = zeros(3,1);
            if (c1 ~= -1)
                color(1,1) = c1;
                color(2:3,1) = fread(fid,2,'int');
                dim = [352 352 350];
                vox = [0.7 0.7 0.7];
                orient = [1 2 3];
            else
                orient = fread(fid, 3, 'int');
                dim = fread(fid,3,'int');
                vox = fread(fid,3,'float');
                color = fread(fid,3, 'int');
            end
            
            vertices = fread(fid,numverts*3,'float');
            triangles_tmp = fread(fid, numtris*4, 'uint32');
            triangles_tmp = reshape(triangles_tmp, [4, numtris]);
            triangles = triangles_tmp(2:end, :);
            
%             for i = 1:numtris
%                 % if R-style binary mesh, uncomment the following line:
%                 tmp = fread(fid, 1, 'uint32')
%                 for j = 1:3
%                     if feof(fid)
%                         break;
%                     end
%                     triangles(j,i) = fread(fid, 1, 'uint32');
%                 end
%             end
%             triangles = fread(fid,numtris*3,'int');
            fclose(fid);
            
            vertices = reshape(vertices, [3,numverts])';
            triangles = triangles';
%             triangles = reshape(triangles, [3, numtris])';
            triangles = triangles + ones(size(triangles));
%             vertices(:,2) = 255-vertices(:,2);
            surf = surfstruct;

            surf.id = id;
            surf.numverts = length(vertices);
            surf.numtris = length(triangles);
            surf.vertices = vertices;
            surf.triangles = triangles;
            surf.dim = dim;
            surf.vox = vox;
            surf.orient = orient;
            surf.color = color;
        end
        
        function surf = readSurfAsSTLMesh(filename)
            fid = fopen(filename,'rt');
            s = fgetl(fid);
            mesh = zeros(3,100000);
            count = 0;
            while (length(s) < 8 || (length(s) >= 8  && ~strcmp(s(1:8),'endsolid')))
                if isempty(s)
                    s = fgetl(fid);
                    continue;
                end
                i = 1;
                while (i<length(s))&&(s(i)==' ')
                    i = i+1;
                end
                if (s(i) == 'v')
                    count = count+1;
                    mesh(:,count) = sscanf(s(i+6:end),'%f %f %f');
                end
                s = fgetl(fid);
            end
            fclose(fid);
            mesh = mesh(:,1:count);

            m.vertices = zeros(3,count);
            m.triangles = zeros(3,count/3);

            vertnum=3;
            m.triangles(:,1) = [1;2;3];
            m.vertices(:,1:3) = mesh(:,1:3);
            for i=2:count/3
                for k=1:3
                    j = find( (mesh(1,(i-1)*3+k)==m.vertices(1,1:vertnum))&...
                              (mesh(2,(i-1)*3+k)==m.vertices(2,1:vertnum))&...
                              (mesh(3,(i-1)*3+k)==m.vertices(3,1:vertnum)));
                    if (length(j)==0)
                        vertnum = vertnum+1;
                        m.vertices(:,vertnum) = mesh(:,(i-1)*3+k);
                        m.triangles(k,i) = vertnum;
                    else
                        m.triangles(k,i) = j(1);
                    end
                end
            end
            m.numtris = count/3;
            m.numverts = vertnum;
            m.vertices = m.vertices(:,1:vertnum);
            
            surf = surfstruct;
            surf.numverts = m.numverts;
            surf.numtris = m.numtris;
            surf.vertices = m.vertices';
            surf.triangles = m.triangles';
            surf.isMm = 1;
            
        end
        
        function saveSurfAsMesh(surf, filename, flag)
            surf = Surfacetools.convertToVox(surf);
            % save surface to a file as mesh
            surf.triangles = surf.triangles - ones(size(surf.triangles));
            surf.vertices = surf.vertices';
            surf.triangles = surf.triangles';
            numvrts = length(surf.vertices);
            numtris = length(surf.triangles);
            
            surf.vertices = surf.vertices(:);
            tvertices = surf.vertices;
            surf.triangles = surf.triangles(:);
            
            if nargin < 3
                flag = 0;
            end
            
            v = [1; 2; 3];
%             if flag
%                 for i = 1:length(v)
%                     v(i) = abs(surf.orient(i));
%                 end
%             end
            for j = 1:3
                for i = 1:numvrts-1
                    %            if(orient(j)>0)
                    tvertices(3*i+j) = surf.vertices(3*i+v(j));
                    %            else
                    %                tvertices(3*i+j) = dim(v(j))-1-vertices(3*i+v(j));
                    %            end
                end
            end
                  
            fid = fopen(filename,'wb');
            id=0;
            fwrite(fid, id,'int');
            fwrite(fid, numvrts,'int');
            fwrite(fid, numtris, 'int');
            if flag
                id = -1;
                fwrite(fid, id,'int');
                fwrite(fid, surf.orient,'int');
                fwrite(fid, surf.dim, 'int');
                fwrite(fid, surf.vox,'float');
            end
            fwrite(fid, surf.color,'int');
            fwrite(fid, tvertices,'float');
            fwrite(fid, surf.triangles,'int');
            fclose(fid);
        end
        
        function saveSurfAsRMesh(surf, filename)  
            surf = Surfacetools.convertToVox(surf);
            % save surface to a file as mesh
            surf.triangles = surf.triangles - ones(size(surf.triangles));
            surf.vertices = surf.vertices';
            surf.triangles = surf.triangles';
            numvrts = length(surf.vertices);
            numtris = length(surf.triangles);
            
            triangles_tmp = zeros(4, size(surf.triangles, 2));
            %triangles_tmp(1, :) = 0:size(surf.triangles, 2)-1;
            triangles_tmp(1, :) = 3;
            triangles_tmp(2:end, :) = surf.triangles;
            triangles_tmp = triangles_tmp(:);
            
            surf.vertices = surf.vertices(:);
%             surf.triangles = surf.triangles(:);
    
            fid = fopen(filename,'wb');
            id=0;
            fwrite(fid, id,'int');
            fwrite(fid, numvrts,'int');
            fwrite(fid, numtris, 'int');
            fwrite(fid, surf.color,'int');
            fwrite(fid, surf.vertices,'float');
            fwrite(fid, triangles_tmp,'int');
            fclose(fid);
        end
        
        function saveSurfAsSTLMesh(surf, filename)
            surf =  Surfacetools.convertToMm(surf);
            fid = fopen(filename,'wt');
            fprintf(fid,'solid mesh\n');
            for i=1:surf.numtris
                c = cross(surf.vertices(surf.triangles(i, 2), :)-surf.vertices(surf.triangles(i, 1), :),surf.vertices(surf.triangles(i, 3), :)-surf.vertices(surf.triangles(i, 1), :));
                c = c/norm(c);
                fprintf(fid,'facet normal %f %f %f\n  outer loop\n',c);
                fprintf(fid,'    vertex    %f  %f  %f\n', surf.vertices(surf.triangles(i, 1), :));
                fprintf(fid,'    vertex    %f  %f  %f\n', surf.vertices(surf.triangles(i, 2), :));
                fprintf(fid,'    vertex    %f  %f  %f\n', surf.vertices(surf.triangles(i, 3), :));
                fprintf(fid,'  endloop\nendfacet\n');
            end
            fprintf(fid,'endsolid mesh');
            fclose(fid);
            
        end
        
        function saveSurfAsVTK(surf, filename)
            % save surface to a file as mesh
            surf.triangles = surf.triangles - ones(size(surf.triangles));
            surf.vertices = surf.vertices';
            surf.triangles = surf.triangles';
            numvrts = length(surf.vertices);
            numtris = length(surf.triangles);
           
            triangles_tmp = zeros(4, size(surf.triangles, 2));
            triangles_tmp(2:end, :) = surf.triangles;
            triangles_tmp = triangles_tmp(:);
            
            surf.vertices = surf.vertices';
%             surf.triangles = surf.triangles(:);
    
            fid = fopen(filename,'wb');
            if(fid<0)
                fprintf('could not open file %s for writing\n',filename);
                return
            end

            % write header
            dlmwrite(filename, [ '# vtk DataFile Version 3.0' ],'delimiter','');
            dlmwrite(filename, [ 'DataEntityType: Surface mesh' ],'delimiter','','-append');
            dlmwrite(filename, [ 'ASCII' ],'delimiter','','-append');
            dlmwrite(filename, [ 'DATASET POLYDATA' ],'delimiter','','-append');

            % write points
            dlmwrite(filename, [ 'POINTS ' num2str(numvrts) ' float' ],'delimiter','','-append');
            dlmwrite(filename, surf.vertices,'delimiter',' ','-append');
            %write triangles
            dlmwrite(filename, [ 'POLYGONS ' num2str(numtris) ' ' num2str(4*numtris)  ],'delimiter','','-append');
            dlmwrite(filename, [ ones(numtris,1)*3 surf.triangles'] ,'delimiter',' ','-append');
      
            fclose(fid);
        end
        
        function surfLine = createSurfAsLine(surfin, startPts, endPts, color)
            surfLine = surfin;
            if(nargin > 3)
                surfLine.color = color;
            end
            
            if size(surfin.vertices, 1) == 1
                if size(startPts, 1) == 1
                    startPts = repmat(startPts, [5, 1]);
                    endPts = repmat(endPts, [5, 1]);
                end
                surfLine.numverts = size(startPts, 1) * 2;
                surfLine.numtris = size(startPts, 1) ;
            else 
                surfLine.numverts = 2*surfin.numverts;
                surfLine.numtris = surfin.numverts;
            end
            
            surfLine.vertices = zeros(surfLine.numverts,3);
            surfLine.vertices(1:2:end-1,:) = startPts;
            surfLine.vertices(2:2:end,:) = endPts;
            
            surfLine.triangles = zeros(surfLine.numtris, 2);
            surfLine.triangles(:,1) = [1:2:size(startPts, 1)*2-1]';
            surfLine.triangles(:,2) = [2:2:size(startPts, 1)*2]';
%             surfLine.triangles(:,1) = transpose(1:surfin.numverts);
%             surfLine.triangles(:,2) = surfin.numverts + transpose(1:surfin.numverts);
        end
        
        function sendAsSurfToMeshEditor(surfin, id, caption, colortype)
            
            if(surfin.isMm == 1)
                surfin = Surfacetools.convertToVox(surfin);
            end
            
            if nargin == 1
                ShowInMeshEditor(surfin,'mesh',1000);
            elseif nargin == 2
                ShowInMeshEditor(surfin,'mesh',id);
            elseif nargin == 3
                ShowInMeshEditor(surfin,'mesh',id,caption);
            elseif nargin == 4
                surfin.color = colortype;
                ShowInMeshEditor(surfin,'mesh',id,caption);
            else
                ShowInMeshEditor(surfin,'mesh',0,'surfin',[255 0 0]);
            end
            
        end
        
        function sendAsNormalsToMeshEditor(surfin, normLength, id, caption,colortype)
            
            if(surfin.isMm == 0)
                surfin = Surfacetools.convertToMm(surfin);
            end
            
            startPts = surfin.vertices;
            endPts = surfin.vertices + normLength*surfin.vertexNormals;
            
            surfin = Surfacetools.createSurfAsLine(surfin,startPts, endPts);
            
            if(surfin.isMm == 1)
                surfin = Surfacetools.convertToVox(surfin);
            end
            
            if nargin == 2
                Surfacetools.sendAsLinesToMeshEditor(surfin, 1000 );
            elseif nargin == 3
                Surfacetools.sendAsLinesToMeshEditor(surfin, id );
            elseif nargin == 4
                Surfacetools.sendAsLinesToMeshEditor(surfin, id, caption);
            elseif nargin == 5
                surfin.color = colortype;
                Surfacetools.sendAsLinesToMeshEditor(surfin, id, caption);
            else
                Surfacetools.sendAsLinesToMeshEditor(surfin, 0, 'surfin', [255 0 0]);
            end
            
        end
        
        function sendAsLinesToMeshEditor(surfin, id, caption,colortype)
            if nargin == 1
                ShowInMeshEditor(surfin,'lines',1000);
            elseif nargin == 2
                ShowInMeshEditor(surfin,'lines',id);
            elseif nargin == 3
                ShowInMeshEditor(surfin,'lines',id,caption);
            elseif nargin == 4
                surfin.color = colortype;
                ShowInMeshEditor(surfin,'lines',id,caption);
            else
                ShowInMeshEditor(surfin,'lines',0,'surfin',[255 0 0]);
            end
            
        end
        
        function sendAsPointsToMeshEditor(surfin, id, caption,colortype)
            
            if(surfin.isMm == 1)
                surfin = Surfacetools.convertToVox(surfin);
            end
            
            surfin.vertices = surfin.vertices';
            %Showing lines in mesh editor
            if nargin == 1
                ShowInMeshEditor(surfin,'points',1000);
            elseif nargin == 2
                ShowInMeshEditor(surfin,'points',id);
            elseif nargin == 3
                ShowInMeshEditor(surfin,'points',id,caption);
            elseif nargin == 4
                surfin.color = colortype;
                ShowInMeshEditor(surfin,'points',id,caption);
            else
                ShowInMeshEditor(surfin,'points',0,'surfin',[255 0 0]);
            end
            
        end
        
        function labelim = createLabelMap(surflist, surfthal, order)
            if nargin < 3
                order = [4, 2, 7, 9, 10, 11, 22, 3, 23, 1, 6, 5, 18, 17, 16, 15, 12, 13, 14, 21, 20, 19, 8];
%                 {'AVD'; 'AM'; 'LD'; 'MTT'; 'MD'; 'CM'; 'Pf'; 'CL'; 'CeM'; 'Pv'; 'Hb'; 'VA'; 'VLa'; 'VLp'; 'VM'; 'VPL'; 'VPM'; 'VPI'; 'PuMI'; 'PuA'; 'PuL'; 'Li'; 'LP'}
%                     1      2      3      4      5      6      7      8      9      10    11    12    13      14    15     16      17     18     19       20      21   22    23
            end
            dim = surfthal.dim;
            labelim = zeros(dim(1), dim(2), dim(3));
            labelim = labelim(:);
            for i = 1:length(surflist)
                j = order(i);
                mask = Surfacetools.convertMeshToMask(surflist{j});
                labelim(intersect(find(mask.data > 127), find(labelim == 0))) = j * 10;
            end
            for i = 1:length(surflist)
                j = order(i);
                mask = Surfacetools.convertMeshToMask(surflist{j});
                labelim(intersect(find(mask.data > 50), find(labelim == 0))) = j * 10;
            end
            mask = Surfacetools.convertMeshToMask(surfthal);
            labelim(find(mask.data < 127)) = 0;
        end
        
        function labelim = createLabelMapDistMap(surflist, surfthal)
            bbox = zeros(2, 3);
            bbox(1, :) = floor(min(surfthal.vertices)) - 5;
            bbox(2, :) = ceil(max(surfthal.vertices)) + 5;
            bdim = bbox(2, :) - bbox(1, :) + 1;
            fpoint = bbox(1, :) - 1;
            dim = surfthal.dim;
            threshold = 7;
            
            distmaps = zeros(length(surflist), bdim(1)*bdim(2)*bdim(3));
            for i = 1:length(surflist)
                surftmp = Surfacetools.smooth(surflist{i});
                surftmp = Surfacetools.subtractFirstPoint(surftmp, fpoint);
                surftmp.dim = bdim;
                distmap = Surfacetools.computeDistanceMap(surftmp, threshold);
                distmaps(i, :) = distmap.data;
            end
            
            [minV, minInd] = min(distmaps);
            blabel = minInd * 10;
            clear minV minInd distmaps;
            
            surfthal2 = Surfacetools.subtractFirstPoint(surfthal, fpoint);
            surfthal2.dim = bdim;
            mask = Surfacetools.convertMeshToMask(surfthal2); % create bdim image if converting manually
            blabel(find(mask.data < 127)) = 0; % both mask and blabel should have bdim
            blabel = reshape(blabel, bdim);
            
            labelim = zeros(dim(1), dim(2), dim(3));
            labelim(bbox(1,1):bbox(2, 1), bbox(1, 2):bbox(2, 2), bbox(1,3):bbox(2, 3)) = blabel;
            labelim = labelim(:);
        end
        
        function surfout = addFirstPoint(surfin, fpoint)
            surfout = surfin;
            surfout.vertices = surfout.vertices + repmat(transpose(fpoint(:)), length(surfout.vertices), 1);
        end
        
        function surfout = subtractFirstPoint(surfin, fpoint)
            surfout = surfin;
            surfout.vertices = surfout.vertices - repmat(transpose(fpoint(:)), length(surfout.vertices), 1);
        end
        
        function surfout = updateSurface(surfin, dim, vox, orient, color)
            %surfout = updateSurface(surfin, newDim, newVox, newOrient, newColor, fPoint)
            if nargin == 2
                surfout = Surfacetools.setSurfProporties(surfin, dim, surfin.vox, surfin.orient, surfin.color);
            elseif nargin == 3
                surfout = Surfacetools.setSurfProporties(surfin, dim, vox, surfin.orient,surfin.color);
            elseif nargin == 4
                surfout = Surfacetools.setSurfProporties(surfin, dim, vox, orient,surfin.color);
            else
                surfout = Surfacetools.setSurfProporties(surfin, dim, vox, orient,color);
            end
            
            surfout.vertices = surfout.vertices.*repmat(transpose(surfin.vox(:)), length(surfout.vertices),1);
            surfout.vertices = surfout.vertices./repmat(transpose(surfout.vox(:)), length(surfout.vertices),1);
            
        end
        
        function surflist = separateSurfs(surfin, numvrtslist, numtrislist,colorlist)
            surflist = cell(length(numvrtslist),1);
            clrlist = cell(length(numvrtslist),1);
            if nargin < 4
                for ii = 1:length(numvrtslist)
                    clrlist{ii} =  surfin.color;
                end
            else
                for ii = 1:length(numvrtslist)
                    clrlist{ii} =  colorlist{ii};
                end
            end
            
            vrtStep = 1;
            trisStep = 1;
            for ii = 1:length(numvrtslist)
                surflist{ii} = surfin;
                vertices = surfin.vertices(vrtStep:numvrtslist{ii} + (vrtStep-1),:);
                triangles = surfin.triangles(trisStep:numtrislist{ii} + (trisStep-1),:);
                triangles = triangles - (vrtStep - 1);
                surflist{ii}.vertices = vertices;
                surflist{ii}.numverts = length(vertices);
                surflist{ii}.triangles = triangles ;
                surflist{ii}.numtris = length(triangles);
                surflist{ii} = Surfacetools.setSurfProporties(surflist{ii}, surfin.dim, surfin.vox, surfin.orient, clrlist{ii});
                surflist{ii}.isMm = surfin.isMm;
                vrtStep = vrtStep + numvrtslist{ii};
                trisStep = trisStep +  numtrislist{ii};
            end
            
        end
        
        % take surflist and return combsurf
        function combsurf = combineSurfs(surflist)
            % initialise surface to combine
            combsurf = surfstruct;
            
            numvertslist = zeros(length(surflist), 1);
            numtrislist = zeros(length(surflist), 1);
            for ii = 1:length(surflist)
                numvertslist(ii,1) = length(surflist{ii}.vertices);
                numtrislist(ii,1) = length(surflist{ii}.triangles);
            end
            combsurf.vertices = zeros(sum(numvertslist), 3);
            combsurf.triangles = zeros(sum(numtrislist), 3);
            
            cumnumverts = cumsum(numvertslist);
            cumnumtris = cumsum(numtrislist);
            
            % concatenate surface vertices
            for ii = 1:length(surflist)
                surflist{ii} = Surfacetools.convertToVox(surflist{ii});
                if(ii == 1)
                    combsurf.vertices(1:numvertslist(ii),:) = surflist{ii}.vertices;
                    combsurf.triangles(1:numtrislist(ii),:) = surflist{ii}.triangles;
                else
                    combsurf.vertices(cumnumverts(ii-1) + 1:cumnumverts(ii),:) = surflist{ii}.vertices;
                    combsurf.triangles(cumnumtris(ii-1) + 1:cumnumtris(ii),:) = surflist{ii}.triangles + cumnumverts(ii-1)*ones(size(surflist{ii}.triangles));
                end
            end
            combsurf.numverts = sum(numvertslist);
            combsurf.numtris = sum(numtrislist);
            combsurf.dim = surflist{1}.dim;
            combsurf.vox = surflist{1}.vox;
            combsurf.orient = surflist{1}.orient;
            combsurf.color = surflist{1}.color;
        end
        
        function [Dice, SurfErr1, SurfErr2] = compareTwoSurfs(surf1, surf2, options)
            if nargin < 3 && nargout == 3
                options.Dice = 'TRUE';
                options.SurfaceError = 'TRUE';
            else if nargin < 3 && nargout == 1
                    options.Dice = 'TRUE';
                    options.SurfaceError = 'FALSE';
                else if nargin < 3 && nargout == 2
                        options.Dice = 'TRUE';
                        options.SurfaceError = 'TRUE';
                    end
                end
            end
                
            if strcmpi(options.Dice, 'TRUE')
                Dice = Surfacetools.calDiceCoeff(surf1, surf2);
            else
                Dice = -1;
            end
            
            if strcmpi(options.SurfaceError, 'TRUE')
                SurfErr1 = Surfacetools.calSurfaceError(surf1, surf2);
                if nargout > 2
                    SurfErr2 = Surfacetools.calSurfaceError(surf2, surf1);
                else
                    SurfErr2 = 0;
                end
            else
                SurfErr1 = 0;
                SurfErr2 = 0;
            end
            
        end
        
        function Dice = calDiceCoeff(surf1, surf2)
            Surfacetools.saveSurfAsRMesh(surf1, 'TMP1.mesh');
            system(['RMesh2Mask.exe TMP1.mesh ' int2str(surf1.dim(1)) ' ' int2str(surf1.dim(2)) ' ' int2str(surf1.dim(3)) ' '  ...
                num2str(surf1.vox(1)) ' ' num2str(surf1.vox(2)) ' ' num2str(surf1.vox(3)) ' ' ' TMP1.mask']);
            Surfacetools.saveSurfAsRMesh(surf2, 'TMP2.mesh');
            system(['RMesh2Mask.exe TMP2.mesh ' int2str(surf2.dim(1)) ' ' int2str(surf2.dim(2)) ' ' int2str(surf2.dim(3)) ' '  ...
                num2str(surf2.vox(1)) ' ' num2str(surf2.vox(2)) ' ' num2str(surf2.vox(3)) ' ' ' TMP2.mask']);         
            mask1 = imgstruct;
            mask1.type = 1;
            mask1.dim = surf1.dim;
            mask1 = Imagetools.readImage(mask1, 'TMP1.mask');
            mask2 = imgstruct;
            mask2.type = 1;
            mask2.dim = surf2.dim;
            mask2 = Imagetools.readImage(mask2, 'TMP2.mask');
            Dice = Imagetools.calDiceCoeff(mask1, mask2);
        end
        
        function SurfErr = calSurfaceError(surfin, refsurf)
            SurfErr = zeros(refsurf.numverts, 1);
            if (isempty(surfin.vertexNtriangles) || ~iscell(surfin.vertexNtriangles))
                surfin.vertexNtriangles = Surfacetools.getVertexNtriangles(surfin);
            end
            for i = 1:refsurf.numverts
                SurfErr(i) = Surfacetools.pointToSurfaceDistance(surfin, refsurf.vertices(i, :));
            end
        end     
        
        function distmap = computeDistanceMap(surfin, threshold)
            if nargin < 2
                threshold = 10;
            end
            Surfacetools.saveSurfAsMesh(surfin, 'TMP.mesh', 0);
            system(['Mesh2DistMap.exe TMP.mesh ' int2str(surfin.dim(1)) ' ' int2str(surfin.dim(2)) ' ' int2str(surfin.dim(3)) ' '  ...
                num2str(surfin.vox(1)) ' ' num2str(surfin.vox(2)) ' ' num2str(surfin.vox(3)) ' ' int2str(threshold) ' TMP.im']);
            distmap = imgstruct;
            distmap.type = 4;
            distmap.dim = surfin.dim;
            distmap.dimx = distmap.dim(1);
            distmap.dimxy =  distmap.dim(1)* distmap.dim(2);
            distmap.dimxyz =  distmap.dim(1)* distmap.dim(2)* distmap.dim(3);
            distmap.vox = surfin.vox;
            distmap.orient = surfin.orient;
            distmap = Imagetools.readImage(distmap, 'TMP.im');
        end
        
        function com = computeCenterOfMass(surfin)
            surfout = Surfacetools.convertToMm(surfin);
            com =mean(surfout.vertices);
        end
        
        function [surfoutlist, surfout] = groupSurfsNew(surflist, threshold)
             if nargin < 2
                threshold = -0.25;
             end
             if length(surflist) == 1
                 surfoutlist = surflist;
                 surfoutlist{1}.deactivateFlag = zeros(surfoutlist{1}.numverts, 1);
                 surfout = surflist{1};
                 return;
             end
            stepsize = 0.7;
            dist_threshold = 5;
            updatedistmap = Surfacetools.computeDistanceMap(surflist{1}, dist_threshold);
            bd = [min(surflist{1}.vertices(:,1)), min(surflist{1}.vertices(:,2)), min(surflist{1}.vertices(:,3)); ...
                max(surflist{1}.vertices(:,1)), max(surflist{1}.vertices(:,2)), max(surflist{1}.vertices(:,3))];
            for i = 2:length(surflist)
                distmap = Surfacetools.computeDistanceMap(surflist{i}, dist_threshold);
                updatedistmap.data(distmap.data<updatedistmap.data) = distmap.data(distmap.data<updatedistmap.data);
                bd = [min(bd(1, 1), min(surflist{i}.vertices(:,1))), min(bd(1, 2), min(surflist{i}.vertices(:,2))), min(bd(1, 3), min(surflist{i}.vertices(:,3))); ...
                max(bd(2, 1), max(surflist{i}.vertices(:,1))), max(bd(2, 2), max(surflist{i}.vertices(:,2))), max(bd(2, 3), max(surflist{i}.vertices(:,3)))];
            end
            bd(1, 1) = max(round(bd(1, 1)) - 4, 1);
            bd(1, 2) = max(round(bd(1, 2)) - 4, 1);
            bd(1, 3) = max(round(bd(1, 3)) - 4, 1);
            bd(2, 1) = min(round(bd(2, 1)) + 6, updatedistmap.dim(1));
            bd(2, 2) = min(round(bd(2, 2)) + 6, updatedistmap.dim(2));
            bd(2, 3) = min(round(bd(2, 3)) + 6, updatedistmap.dim(3));
            openim = Imagetools.open(updatedistmap, bd);
            openim = Imagetools.flipim(openim);
            updatedistmap = Imagetools.calDistMap(openim, 2);
            surfoutlist = cell(length(surflist), 1);
            numverts = 0;
            alldeactivateFlag = [];
            for i = 1:length(surflist)
                tmpsurf = surflist{i};
                vertexNormals  = Surfacetools.computeVertexNormals(tmpsurf);
                dists = zeros(tmpsurf.numverts, 1);
                tmpsurf.deactivateFlag = zeros(tmpsurf.numverts, 1);
                for j = 1:tmpsurf.numverts
                    dists(j) = Surfacetools.pointToSurfaceDistance(tmpsurf, tmpsurf.vertices(j, :), updatedistmap);
                    if dists(j) < threshold 
                        tmpsurf.deactivateFlag(j) = 1;
                        p1 = tmpsurf.vertices(j, :) +  stepsize * vertexNormals(j, :);
                        p2 = tmpsurf.vertices(j, :) + 8 * stepsize * vertexNormals(j, :);
                        dist1 = Surfacetools.pointToSurfaceDistance(tmpsurf, p1, updatedistmap);
                        dist2 = Surfacetools.pointToSurfaceDistance(tmpsurf, p2, updatedistmap);
                        if dist1 > 0 && dist2 > 0
                            tmpsurf.deactivateFlag(j) = 0;
                        end
                    end
                end
                surfoutlist{i} = tmpsurf;
                numverts = numverts + surfoutlist{i}.numverts;
                alldeactivateFlag = [alldeactivateFlag; tmpsurf.deactivateFlag];
            end
            surfout = Surfacetools.combineSurfs(surfoutlist);
            mapping = zeros(numverts, 1);
            curmap = 1;
            for i = 1:length(mapping)
                if alldeactivateFlag(i)
                    continue;
                end
                mapping(i) = curmap;
                curmap = curmap + 1;
            end
            surfout.vertices = surfout.vertices(~alldeactivateFlag, :);
            surfout.numverts = size(surfout.vertices, 1);
            newtriangles = [];
            for i = 1:size(surfout.triangles, 1)
                if alldeactivateFlag(surfout.triangles(i, 1)) || alldeactivateFlag(surfout.triangles(i, 2)) || alldeactivateFlag(surfout.triangles(i, 3)) 
                    continue;
                end
                newtriangles = [newtriangles; mapping(surfout.triangles(i, 1)), mapping(surfout.triangles(i, 2)), mapping(surfout.triangles(i, 3))];
            end
            surfout.triangles = newtriangles;
            surfout.numtris = size(surfout.triangles, 1);  
            % Surfacetools.sendAsSurfToMeshEditor(surfout, 100);
        end
        
        function [surfoutlist, surfout] = groupSurfsExt(surfref, surflist, threshold)
             if nargin < 3
                threshold = 0.25;
             end
            stepsize = 0.7;
            dist_threshold = 10;
            updatedistmap = Surfacetools.computeDistanceMap(surfref, dist_threshold);
            surfoutlist = cell(length(surflist), 1);
            numverts = 0;
            alldeactivateFlag = [];
            for i = 1:length(surflist)
                tmpsurf = surflist{i};
                dists = zeros(tmpsurf.numverts, 1);
                for j = 1:tmpsurf.numverts
                    if tmpsurf.deactivateFlag(j) == 1
                        continue;
                    end
                    dists(j) = Surfacetools.pointToSurfaceDistance(tmpsurf, tmpsurf.vertices(j, :), updatedistmap);
                    if dists(j) > threshold 
                        tmpsurf.deactivateFlag(j) = 1;
                    end
                end
                surfoutlist{i} = tmpsurf;
                numverts = numverts + surfoutlist{i}.numverts;
                alldeactivateFlag = [alldeactivateFlag; tmpsurf.deactivateFlag];
            end
            surfout = Surfacetools.combineSurfs(surfoutlist);
            mapping = zeros(numverts, 1);
            curmap = 1;
            for i = 1:length(mapping)
                if alldeactivateFlag(i)
                    continue;
                end
                mapping(i) = curmap;
                curmap = curmap + 1;
            end
            surfout.vertices = surfout.vertices(~alldeactivateFlag, :);
            surfout.numverts = size(surfout.vertices, 1);
            newtriangles = [];
            for i = 1:size(surfout.triangles, 1)
                if alldeactivateFlag(surfout.triangles(i, 1)) || alldeactivateFlag(surfout.triangles(i, 2)) || alldeactivateFlag(surfout.triangles(i, 3)) 
                    continue;
                end
                newtriangles = [newtriangles; mapping(surfout.triangles(i, 1)), mapping(surfout.triangles(i, 2)), mapping(surfout.triangles(i, 3))];
            end
            surfout.triangles = newtriangles;
            surfout.numtris = size(surfout.triangles, 1);  
            % Surfacetools.sendAsSurfToMeshEditor(surfout, 100);
        end
        
        function [surfoutlist, surfout] = groupSurfs(surflist, threshold)
             if nargin < 2
                threshold = 1;
             end
             if length(surflist) == 1
                 surfoutlist = surflist;
                 surfoutlist{1}.deactivateFlag = zeros(surfoutlist{1}.numverts, 1);
                 surfout = surflist{1};
                 return;
             end
             surfoutlist = cell(length(surflist), 1);
            dist_threshold = 5;
            distmap1 = Surfacetools.computeDistanceMap(surflist{1}, dist_threshold);
            updatedistmap = distmap1;
            allnumverts = zeros(length(surflist), 1);
            for i = 1:length(surflist)-1
                allnumverts(i) = surflist{i}.numverts;
                if i == 1
                    [surfoutlist{i}, surfoutlist{i+1}, updatedistmap, surfout] = Surfacetools.groupTwoSurfs(surflist{i}, surflist{i+1}, updatedistmap, threshold);
                else
                    [tmp, surfoutlist{i+1}, updatedistmap, surfout] = Surfacetools.groupTwoSurfs(surfout, surflist{i+1}, updatedistmap, threshold);
                    % update deactivateFlag for previously grouped structures surfout
                    alldeactivateFlag = zeros(sum(allnumverts(1:i)), 1);
                    startIndex = 1;
                    for j = 1:i
                        alldeactivateFlag(startIndex:startIndex+allnumverts(j)-1) =surfoutlist{j}.deactivateFlag; 
                        startIndex = startIndex+allnumverts(j);
                    end
                    count = 1;
                    for j = 1:length(alldeactivateFlag)
                        if alldeactivateFlag(j) == 0
                            alldeactivateFlag(j) = tmp.deactivateFlag(count);
                            count = count + 1;
                        end
                    end
                    count = 1;
                    for j = 1:i
                        surfoutlist{j}.deactivateFlag = alldeactivateFlag(count:count+surfoutlist{j}.numverts-1);
                        count = count+surfoutlist{j}.numverts;
                    end
                end
            end            
        end
        
        % group two surfaces with adjacent boundaries
        function [surfout1, surfout2, newdistmap, surfout] = groupTwoSurfs(surf1, surf2, distmap1, threshold)        
            if nargin < 3
                threshold = 1;
            end
            
            dist_threshold = 5;
            distmap2 = Surfacetools.computeDistanceMap(surf2, dist_threshold);
            
            newdistmap = distmap1;
            newdistmap.data(distmap2.data<distmap1.data) = distmap2.data(distmap2.data<distmap1.data);
            
            dists1 = zeros(surf1.numverts, 1);
            dists2 = zeros(surf2.numverts, 1);
            
            surf1 = Surfacetools.convertToVox(surf1);
            surf2 = Surfacetools.convertToVox(surf2);

            for i = 1:surf1.numverts
                dists1(i) = Surfacetools.pointToSurfaceDistance(surf2, surf1.vertices(i, :), distmap2);
            end

            for i = 1:surf2.numverts
                dists2(i) = Surfacetools.pointToSurfaceDistance(surf1, surf2.vertices(i, :), distmap1);
            end
            
            surfout1 = surf1;
            surfout2 = surf2;
            surfout1.deactivateFlag = dists1<threshold;
            surfout2.deactivateFlag = dists2<threshold;
            
            if nargout <=3
                return;
            end
            
            % compute surfout
            % boundary1 = Surfacetools.findBoundaryPoints(surfout1);
            % boundary2 = Surfacetools.findBoundaryPoints(surfout2);
            surfout = Surfacetools.combineSurfs({surf1, surf2});
            % delete deactivated vertices and associated triangles
            mapping = zeros(surf1.numverts+ surf2.numverts, 1);
            alldeactivateFlag = [surfout1.deactivateFlag; surfout2.deactivateFlag];
            curmap = 1;
            for i = 1:length(mapping)
                if alldeactivateFlag(i)
                    continue;
                end
                mapping(i) = curmap;
                curmap = curmap + 1;
            end
            surfout.vertices = surfout.vertices(~alldeactivateFlag, :);
            surfout.numverts = size(surfout.vertices, 1);
            newtriangles = [];
            for i = 1:size(surfout.triangles, 1)
                if alldeactivateFlag(surfout.triangles(i, 1)) || alldeactivateFlag(surfout.triangles(i, 2)) || alldeactivateFlag(surfout.triangles(i, 3)) 
                    continue;
                end
                newtriangles = [newtriangles; mapping(surfout.triangles(i, 1)), mapping(surfout.triangles(i, 2)), mapping(surfout.triangles(i, 3))];
            end
            surfout.triangles = newtriangles;
            surfout.numtris = size(surfout.triangles, 1);  
        end
        
        % apply grouping
        function [surfoutlist, surfout] = applyGrouping(surflist, refsurflist)
            surfoutlist = surflist;
            for i = 1:length(surflist)
                surfoutlist{i}.deactivateFlag = refsurflist{i}.deactivateFlag;
                surfoutlist{i}.isMm = surflist{i}.isMm;
            end
            if nargout < 2
                return;
            end
            surfout = Surfacetools.combineSurfs(surflist);
            % delete deactivated vertices and associated triangles
            mapping = zeros(surfout.numverts, 1);
            alldeactivateFlag = zeros(surfout.numverts, 1);
            count = 1;
            for i = 1:length(surflist)
                alldeactivateFlag(count:count+surflist{i}.numverts-1) = refsurflist{i}.deactivateFlag;
                count = count+surflist{i}.numverts;
            end
            curmap = 1;
            for i = 1:length(mapping)
                if alldeactivateFlag(i)
                    continue;
                end
                mapping(i) = curmap;
                curmap = curmap + 1;
            end
            surfout.vertices = surfout.vertices(~alldeactivateFlag, :);
            surfout.numverts = size(surfout.vertices, 1);
            newtriangles = [];
            for i = 1:size(surfout.triangles, 1)
                if alldeactivateFlag(surfout.triangles(i, 1)) || alldeactivateFlag(surfout.triangles(i, 2)) || alldeactivateFlag(surfout.triangles(i, 3)) 
                    continue;
                end
                newtriangles = [newtriangles; mapping(surfout.triangles(i, 1)), mapping(surfout.triangles(i, 2)), mapping(surfout.triangles(i, 3))];
            end
            surfout.triangles = newtriangles;
            surfout.numtris = size(surfout.triangles, 1);  
            surfout.isMm = surflist{1}.isMm;
            
        end
        
        % update grouped structures vertices for each structure in the
        % group, with deactivatedFlag = 0
        function surfoutlist = updateGroupSurfsVertices(surflist, surfin)
            surfoutlist = surflist;
            allnumverts = zeros(length(surflist), 1);
            for i = 1:length(surflist)
                allnumverts(i) = surflist{i}.numverts;
            end
            count = 1; i = 1;
            for j = 1:sum(allnumverts)
                if j > sum(allnumverts(1:i))
                    i = i+1;
                end
                mapindex = j - sum(allnumverts(1:i-1));
                if surfoutlist{i}.deactivateFlag(mapindex) == 0
                    surfoutlist{i}.vertices(mapindex, :) = surfin.vertices(count, :);
                    count = count + 1;
                end
            end      
        end
 
        % recent update
        function surfout = refineSurfVertices(surfin, refsurflist, flip, distval)
            if nargin < 3
                flip = 1;
            end
            if nargin < 4
                distval = 5;
            end
            surfout = surfin;
            
            distmap = Surfacetools.computeDistanceMap(refsurflist{1}, distval);
            for j = 2:length(refsurflist)
                distmap1 = Surfacetools.computeDistanceMap(refsurflist{j}, distval);
                distmap(distmap1.data < distmap.data) = distmap1.data(distmap1.data < distmap.data);
            end
            distmap.data(distmap.data > 10) = 10;
            vertexNormals  = Surfacetools.computeVertexNormals(surfin);
            
            for i = 1:surfin.numverts
                dist = Surfacetools.pointToSurfaceDistance(refsurf, surfin.vertices(i, :), distmap);
                if dist * flip > 0
                    continue;
                end
                for j = 0.05:0.05:2
                    surfout.vertices(i, :) = surfin.vertices(i, :) - vertexNormals(i, :) * j;
                    dist = Surfacetools.pointToSurfaceDistance(refsurf, surfout.vertices(i, :), distmap);
                    if dist * flip > 0
                        break;
                    end
                end
            end   
        end
        
        function surfout = refineSurfVerticesDistMap(surfin, refsurf, distmap, flipDist, flipNormal)
            if nargin < 4
                flipDist = 1;
            end
            if nargin < 5
                flipNormal = 1;
            end
            surfout = surfin;
            if min(distmap.data) > 0
                return;
            end
            
            distmap.data(distmap.data > 10) = 10;
            vertexNormals  = Surfacetools.computeVertexNormals(surfin);
            
            for i = 1:surfin.numverts
                dist = Surfacetools.pointToSurfaceDistance(refsurf, surfin.vertices(i, :), distmap);
                if dist * flipDist > 0
                    continue;
                end
                for j = 0.05:0.05:2
                    surfout.vertices(i, :) = surfin.vertices(i, :) - flipNormal * vertexNormals(i, :) * j;
                    dist = Surfacetools.pointToSurfaceDistance(refsurf, surfout.vertices(i, :), distmap);
                    if dist * flipDist > 0
                        break;
                    end
                end
            end   
        end
        
        function mask = convertMeshToMask(surf)
            Surfacetools.saveSurfAsRMesh(surf, 'TMP.mesh'); % mesh created has wrong dims and vox
            system(['RMesh2Mask.exe TMP.mesh ' int2str(surf.dim(1)) ' ' int2str(surf.dim(2)) ' ' int2str(surf.dim(3)) ' '  ...
                num2str(surf.vox(1)) ' ' num2str(surf.vox(2)) ' ' num2str(surf.vox(3)) ' ' ' TMP.mask']);
            mask = imgstruct;
            mask.type = 1;
            mask.dim = surf.dim;
            mask.vox = surf.vox;
            mask.orient = surf.orient;
            mask = Imagetools.readImage(mask, 'TMP.mask');
        end
        
        function newsurf = convertMaskToMesh(surf, mask)
            fid = fopen('TMP.mask', 'wb');
            fwrite(fid, mask, 'uint8');
            fclose(fid);
             system(['Mask2Mesh.exe TMP.mask '  int2str(surf.dim(1)) ' ' int2str(surf.dim(2)) ' ' int2str(surf.dim(3)) ' '  ...
                num2str(surf.vox(1)) ' ' num2str(surf.vox(2)) ' ' num2str(surf.vox(3)) ' ' ' TMP.mesh']);
            newsurf = Surfacetools.readSurfAsMesh('TMP.mesh');
            newsurf = Surfacetools.setSurfProporties(newsurf,surf.dim, surf.vox, surf.orient);
            %system('del TMP.mask');
            %system('del TMP.mesh');            
        end  

        function surfout = setSurfProporties(surfin, dim, vox,orient, color)
            
            surfout = surfin;
            if nargin == 1
                surfout.dim = surfin.dim;
                surfout.vox = surfin.vox;
                surfout.orient = surfin.orient;
                surfout.color = surfin.color;
            elseif nargin == 2
                surfout.dim = dim;
            elseif nargin == 3
                surfout.dim = dim;
                surfout.vox = vox;
            elseif nargin == 4
                surfout.dim = dim;
                surfout.vox = vox;
                surfout.orient = orient;
            else
                surfout.dim = dim;
                surfout.vox = vox;
                surfout.orient = orient;
                surfout.color = color;
            end
            
        end
        
        function surfout = transformSurface(surfin, transformation)
            % input surface
            surfout = surfin;
            
            % transform surface
            surfout.vertices = transfpts(surfin.vertices, transformation);
            
        end
        
        function surfout = transformSurfaceMM(surfin, transformation)
            % input surface
            surfout = surfin;
            surfin = Surfacetools.convertToMm(surfin);
            
            % transform surface
            surfout.vertices = transfpts(surfin.vertices, transformation);
            surfout = Surfacetools.convertToVox(surfout);
            
        end
        
        function surfout = affDefSurf1toSurf2(surf1, surf2, transformation)
            % input surface
            surfout = surf1;
            
            % convert both to mm
            surfin = Surfacetools.convertToMm(surf1);
            
            % transform surface
            surfout.vertices = transfpts(surfin.vertices, transformation);
            
            % copy surf2 propoerties to surfout
            surfout  = Surfacetools.setSurfProporties(surfout, surf2.dim, surf2.vox, surf2.orient);
            
            % convert both to mm
            surfout = Surfacetools.convertToVox(surfout);
            
        end
        
        function surfout = flipAxis(surfin, dir)
            if nargin == 1
                dir = 0;
            end
            surfout = Surfacetools.convertToVox(surfin);
            dir = dir + 1;
            surfout.vertices(:, dir) = surfout.dim(dir) - surfout.vertices(:, dir);
        end
        
        function [dist, PP0] = pointToSurfaceDistance(surfin, pt, distmap)
            
            if nargin == 3
                % calculate closest distance according to the given
                % distance map
                dist1 = Imagetools.getIntensityAt(distmap, pt);
            end
            
            if (isempty(surfin.vertexNtriangles) || ~iscell(surfin.vertexNtriangles))
                surfin.vertexNtriangles = Surfacetools.getVertexNtriangles(surfin);
            end
            pt = pt(:); pt = pt';
            dists = surfin.vertices - repmat(pt, [surfin.numverts, 1]);
            dists = sqrt(sum(dists.^2, 2));
            [mindist, ptid] = min(dists);
            ptid = ptid(1);
            
            polys = surfin.vertexNtriangles{ptid}; 
            triangles_to_check = zeros(3*length(polys), 3);
            for i = 1:length(polys)
                tmp = [  surfin.vertices(surfin.triangles(polys(i),1),:); ...
                    surfin.vertices(surfin.triangles(polys(i),2),:); ...
                    surfin.vertices(surfin.triangles(polys(i),3),:) ...
                    ];
                triangles_to_check(3*i-2:3*i, :) = tmp;
            end
            distances_to_check = zeros(length(polys), 1);
            closestpt_to_check = zeros(length(polys), 3);
            for i = 1:length(polys)
                [distances_to_check(i), closestpt_to_check(i, :)] = Surfacetools.pointToTriangleDistance(pt, triangles_to_check(3*i-2:3*i, :));
            end
            [dist, ptid] = min(distances_to_check);         
            if nargout>1
              PP0 = closestpt_to_check(ptid(1), :);
            end 
            if nargin == 3
                dist = dist1;
            end
            
        end
        
        % Calculate the distance of a given point pt from a triangle triangle.
        % pt [] 1*3. triangle: [pt1; pt2; pt3] 3*3
        % additionally could return the closest point PP0 to pt on the
        % triangle.
        function [dist, PP0] = pointToTriangleDistance(pt, triangle)
            
            % rewrite triangle in normal form
            B = triangle(1,:);
            E0 = triangle(2,:)-B;
            %E0 = E0/sqrt(sum(E0.^2)); %normalize vector
            E1 = triangle(3,:)-B;
            %E1 = E1/sqrt(sum(E1.^2)); %normalize vector


            D = B - pt;
            a = dot(E0,E0);
            b = dot(E0,E1);
            c = dot(E1,E1);
            d = dot(E0,D);
            e = dot(E1,D);
            f = dot(D,D);

            det = a*c - b*b; % do we have to use abs here?
            s   = b*e - c*d;
            t   = b*d - a*e;

            % Terible tree of conditionals to determine in which region of the diagram
            % shown above the projection of the point into the triangle-plane lies.
            if (s+t) <= det
              if s < 0
                if t < 0
                  %region4
                  if (d < 0)
                    t = 0;
                    if (-d >= a)
                      s = 1;
                      sqrDistance = a + 2*d + f;
                    else
                      s = -d/a;
                      sqrDistance = d*s + f;
                    end
                  else
                    s = 0;
                    if (e >= 0)
                      t = 0;
                      sqrDistance = f;
                    else
                      if (-e >= c)
                        t = 1;
                        sqrDistance = c + 2*e + f;
                      else
                        t = -e/c;
                        sqrDistance = e*t + f;
                      end
                    end
                  end %of region 4
                else
                  % region 3
                  s = 0;
                  if e >= 0
                    t = 0;
                    sqrDistance = f;
                  else
                    if -e >= c
                      t = 1;
                      sqrDistance = c + 2*e +f;
                    else
                      t = -e/c;
                      sqrDistance = e*t + f;
                    end
                  end
                end %of region 3 
              else
                if t < 0
                  % region 5
                  t = 0;
                  if d >= 0
                    s = 0;
                    sqrDistance = f;
                  else
                    if -d >= a
                      s = 1;
                      sqrDistance = a + 2*d + f;% GF 20101013 fixed typo d*s ->2*d
                    else
                      s = -d/a;
                      sqrDistance = d*s + f;
                    end
                  end
                else
                  % region 0
                  invDet = 1/det;
                  s = s*invDet;
                  t = t*invDet;
                  sqrDistance = s*(a*s + b*t + 2*d) ...
                              + t*(b*s + c*t + 2*e) + f;
                end
              end
            else
              if s < 0
                % region 2
                tmp0 = b + d;
                tmp1 = c + e;
                if tmp1 > tmp0 % minimum on edge s+t=1
                  numer = tmp1 - tmp0;
                  denom = a - 2*b + c;
                  if numer >= denom
                    s = 1;
                    t = 0;
                    sqrDistance = a + 2*d + f; % GF 20101014 fixed typo 2*b -> 2*d
                  else
                    s = numer/denom;
                    t = 1-s;
                    sqrDistance = s*(a*s + b*t + 2*d) ...
                                + t*(b*s + c*t + 2*e) + f;
                  end
                else          % minimum on edge s=0
                  s = 0;
                  if tmp1 <= 0
                    t = 1;
                    sqrDistance = c + 2*e + f;
                  else
                    if e >= 0
                      t = 0;
                      sqrDistance = f;
                    else
                      t = -e/c;
                      sqrDistance = e*t + f;
                    end
                  end
                end %of region 2
              else
                if t < 0
                  %region6 
                  tmp0 = b + e;
                  tmp1 = a + d;
                  if (tmp1 > tmp0)
                    numer = tmp1 - tmp0;
                    denom = a-2*b+c;
                    if (numer >= denom)
                      t = 1;
                      s = 0;
                      sqrDistance = c + 2*e + f;
                    else
                      t = numer/denom;
                      s = 1 - t;
                      sqrDistance = s*(a*s + b*t + 2*d) ...
                                  + t*(b*s + c*t + 2*e) + f;
                    end
                  else  
                    t = 0;
                    if (tmp1 <= 0)
                        s = 1;
                        sqrDistance = a + 2*d + f;
                    else
                      if (d >= 0)
                          s = 0;
                          sqrDistance = f;
                      else
                          s = -d/a;
                          sqrDistance = d*s + f;
                      end
                    end
                  end
                  %end region 6
                else
                  % region 1
                  numer = c + e - b - d;
                  if numer <= 0
                    s = 0;
                    t = 1;
                    sqrDistance = c + 2*e + f;
                  else
                    denom = a - 2*b + c;
                    if numer >= denom
                      s = 1;
                      t = 0;
                      sqrDistance = a + 2*d + f;
                    else
                      s = numer/denom;
                      t = 1-s;
                      sqrDistance = s*(a*s + b*t + 2*d) ...
                                  + t*(b*s + c*t + 2*e) + f;
                    end
                  end %of region 1
                end
              end
            end

            % account for numerical round-off error
            if (sqrDistance < 0)
              sqrDistance = 0;
            end

            dist = sqrt(sqrDistance);

            if nargout>1
              PP0 = B + s*E0 + t*E1;
            end            
            
        end
        
        function area = SurfaceArea(surfin)
            area = 0;
            for polNo = 1 : length(surfin.triangles)
                vertexNo = surfin.triangles(polNo,:);
                vertices = surfin.vertices(vertexNo,:);
                area =area + 0.5 * norm(cross((vertices(2,:)-vertices(1,:)), (vertices(3,:)-vertices(1,:))));
            end
        end
        
        function vol = SurfaceVolume(surfin)
            mask = Surfacetools.convertMeshToMask(surfin);
            vol = sum(mask.data>127) *surfin.vox(1)*surfin.vox(2)*surfin.vox(3);
        end
        
        function boundary = findBoundaryPoints(surfin)
            
            if (~ismatrix(surfin.deactivateFlag))
                flag = 0;
                return;
            end
            
            nvert = surfin.numverts;
            nface=surfin.numtris;

            A=sparse(nvert,nvert);
            for i=1:nface
%                 if verb
%                     progressbar(i,nface);
%                 end
                f=surfin.triangles(i,:);
                if (surfin.deactivateFlag(f(1)) || surfin.deactivateFlag(f(2))  || surfin.deactivateFlag(f(3))  )
                    continue;
                end
                A(f(1),f(2))=A(f(1),f(2))+1;
                A(f(1),f(3))=A(f(1),f(3))+1;
                A(f(3),f(2))=A(f(3),f(2))+1;
            end
            A=A+A';

            for i=1:nvert
                u=find(A(i,:)==1);
                if ~isempty(u)
                    boundary=[i u(1)];
                    break;
                end
            end

            s=boundary(2);
            i=2;
            while(i<=nvert)
                u=find(A(s,:)==1);
                if length(u)~=2
                    warning('problem in boundary');
                end
                if u(1)==boundary(i-1)
                    s=u(2);
                else
                    s=u(1);
                end
                if s~=boundary(1)
                    boundary=[boundary s];
                else
                    break;
                end
                i=i+1;
            end

            if i>nvert
                warning('problem in boundary');
            end         
        end

        function surfout = convertToMm(surfin)
            surfout = surfin;
            if(surfin.isMm == 0)
                surfout.isMm = 1;
                surfout.vertices = surfout.vertices .* repmat(transpose(surfin.vox(:)), surfin.numverts, 1);
            end
        end
        
        function surfout = convertToVox(surfin)
            surfout = surfin;
            if(surfin.isMm == 1)
                surfout.isMm = 0;
                surfout.vertices = surfout.vertices ./ repmat(transpose(surfin.vox(:)), surfin.numverts, 1);
            end
        end
        
        function grads = computeGradients(surfin, img, options)
            if nargin < 3
                options.numSteps = 4;
                options. stepSz = 0.5;
            end
            surfin = Surfacetools.computeVertexNormals2(surfin);       
            surfin.vertexNormals = Imagetools.normalDir(img)*surfin.vertexNormals;
            grads = zeros(surfin.numverts, 1);
            for i = 1:surfin.numverts
                vertexPos = surfin.vertices(i,:);
                vertexNormal = surfin.vertexNormals(i,:);
                currProfile = Imagetools.getProfileAt(img, vertexPos, vertexNormal, options.numSteps, options.stepSz);
                currProfile2 = currProfile(end:-1:1);
                grads(i) = max(abs(currProfile(1:options.numSteps) - currProfile2(1:options.numSteps)));
            end
        end
        
        function newcosts = averageCosts(surfin, costs)
            newcosts = zeros(length(costs), 1);
            if length(costs) ~= surfin.numverts
                disp('Costs length incorrect');
                return;
            end
            oneRingNvertices  =  Surfacetools.getOneRingNvertices(surfin);
            for ii = 1:length(costs)
                newcosts(ii) = mean(costs([ii, oneRingNvertices{ii,1}]));
            end            
        end
        
        function writemeshcolor(costs, pathname)
            
            min_cost = min(costs);
            max_cost = max(costs);
            
%             load objcolormap;
%             objcolormap(1,:) = objcolormap(2,:);
%             mycolormap = objcolormap;
            mycolormap = colormap;
            mycostmap = round((costs-min_cost)/(max_cost - min_cost)*length(mycolormap))-1;
            mycostmap(mycostmap<0)=0;
            
            fp = fopen(pathname, 'wb');
            fwrite(fp, length(mycolormap), 'uint32');
            fwrite(fp, length(costs), 'uint32');
            fwrite(fp, mycolormap', 'double');
            fwrite(fp, mycostmap, 'uint32');
            fclose(fp);
        end
        
        function writeNormalisedColorMap(costs, pathname, min_cost, max_cost)
            
            costs = (costs - min_cost) / (max_cost - min_cost);
            
%             load objcolormap;
%             objcolormap(1,:) = objcolormap(2,:);
%             mycolormap = objcolormap;
            mycolormap = colormap;
            mycostmap = round((costs-min_cost)/(max_cost - min_cost)*length(mycolormap))-1;
            mycostmap(mycostmap<0)=0;
            
            fp = fopen(pathname, 'wb');
            fwrite(fp, length(mycolormap), 'uint32');
            fwrite(fp, length(costs), 'uint32');
            fwrite(fp, mycolormap', 'double');
            fwrite(fp, mycostmap, 'uint32');
            fclose(fp);
        end
        
        function [transf, fre, freComps] = registerSurfceToSurface(Surf1, Surf2, mode)
            if nargin < 3
                mode = 'RIGID';
            end
            
            % convert Surf1 to vox;
            Surf1 = Surfacetools.convertToMm(Surf1);
            Surf2 = Surfacetools.convertToMm(Surf2);
            
            if strcmp(mode, 'SIMILARITY')
                [transf, fre, freComps] = EstimateSimilarityTransformation(Surf1.vertices, Surf2.vertices);
            elseif strcmp(mode, 'P_RIGID')
                [transf, fre, freComps] = EstimateProperRigidTransformation(Surf1.vertices, Surf2.vertices);
            elseif strcmp(mode, 'P_SIMILARITY')
                [transf, fre, freComps] = EstimateProperSimilarityTransformation(Surf1.vertices, Surf2.vertices);
            else
                [transf, fre, freComps] = EstimateRigidTransformation(Surf1.vertices, Surf2.vertices);
            end
        end
        
        function surfout = MeshToMeshRegistration(surfs, surft, options)
            % register surfs to surft
            default_options = struct;
            default_options.num_level = 2;
            default_options.num_iter = [15, 5];
            default_options.smooth = [0, 0.1];
            default_options.curvature = [1, 0.1];
            default_options.finalstep = 0; % ICP
            if nargin < 3
                options = default_options;
            end
            field_names = fieldnames(default_options);
            for i = 1:length(field_names)
                if ~isfield(options, field_names{i})
                    options.(field_names{i}) = default_options.(field_names{i});
                end
            end
            
            Surfacetools.saveSurfAsMesh(surfs, 'source.mesh', 1);
            Surfacetools.saveSurfAsMesh(surft, 'target.mesh', 1);
            command = ['MeshReg.exe source.mesh target.mesh ' int2str(options.num_level) ' [' ];
            for i = 1:options.num_level
                command = [command int2str(options.num_iter(i))];
                if i ~= options.num_level
                    command = [command ','];
                end
            end
            command = [command '] ['];
            for i = 1:options.num_level
                command = [command num2str(options.smooth(i))];
                if i ~= options.num_level
                    command = [command ','];
                end
            end
            command = [command '] ['];
            for i = 1:options.num_level
                command = [command num2str(options.curvature(i))];
                if i ~= options.num_level
                    command = [command ','];
                end
            end
            command = [command '] ' int2str(options.finalstep) ' register.mesh'];
            system(command);
            surfout = Surfacetools.readSurfAsMesh('register.mesh');            
            
        end
        
        function surfout = smooth(surfin, numIter)
            if(nargin < 2)
                numIter = 2;
            end
            Surfacetools.saveSurfAsMesh(surfin, 'tmp.mesh');
            system(['Smooth.exe TMP.mesh '  int2str(numIter) ' ' ' TMP2.mesh']);
            surfout = Surfacetools.readSurfAsMesh('TMP2.mesh');
            surfout = Surfacetools.setSurfProporties(surfout,surfin.dim, surfin.vox, surfin.orient);
        end
        
        function surfout = affinedeformsurf(surfin, AFFINE_TRANSF_FILE, TRANSF_TYPE)
            surfout = surfin;
            transf = readtransf(AFFINE_TRANSF_FILE, TRANSF_TYPE);            
            surfout.vertices = ApplyGeneralTransformation(surfin.vertices, transf);
        end
        
        function surfout = nrdeformsurf(surfin, srcimg, NR_DEFORMATION_FILE, fquite)
            if(nargin < 4) 
                fquite = 1;
            end
            surfout = surfin;
            surfout.vertices = nrdeformpointsstot(surfin.vertices, srcimg, NR_DEFORMATION_FILE, fquite);
        end
        
        function transform = estimateTPSTrans(surf1, surf2, type, lamda, support)
            if nargin < 3
                type = 'Regular';
            end
            if nargin < 4
                lamda = 0.5;  % 0.1
            end
            if nargin < 5
                support = 3; % 1
            end
            
            pt1 = surf1.vertices;
            pt2 = surf2.vertices;

            x1 = pt1 (:,1); y1 = pt1 (:,2); z1 = pt1 (:,3);
            x2 = pt2 (:,1); y2 = pt2 (:,2); z2 = pt2 (:,3);

            N = length(x1);
            A = zeros(N+4, N+4);
            I_f = eye(N);
            err1 = zeros(length(x1),1);
            err2 = zeros(length(x1),1);
            err = err2;
            err(find(err<err1)) = err1(find(err<err1));
            W_inv = diag(err);
            target = zeros(N+4, 3);
            target(1:N, :) = [x2 y2 z2];

            A(1:N, 1) = ones(N,1);
            A(1:N, 2) = x1;
            A(1:N, 3) = y1;
            A(1:N, 4) = z1;

            ref_points = [x1 y1 z1];
            dist = zeros(N, N);
            for i = 1:N
                dist(:,i) = (ref_points(:,1)-ref_points(i,1)).^2 + (ref_points(:,2)- ref_points(i,2)).^2 + (ref_points(:,3)- ref_points(i,3)).^2;
            end
            % different Radial Basis Functions:
            % normal one:
            % R = r^2log(r^2) (r^2: dist)
            if strcmpi(type, 'Regular')
                dist(find(dist==0)) = 1;
                A(1:N,5:N+4) = dist.*log(dist);
            % from Rohr's paper:
            % R = -1/(8*pi)*r
            % A(1:N,4:N+3) = -1/(8*pi)*sqrt(dist);
            % from Rohr's paper:
            % R = -1/(8*pi)*(4*pi*sigma*G(r/sigma) + (r^2+sigma^2/r)*erf(r/sqrt(2)/sigma))
            % error function erf(r) = (2/sqrt(pi)) integral (exp(-epsilon^2) from 0 to r)
            else if strcmpi(type, 'Gaussian')
                dist(find(dist==0)) = 0.001;
                sigma = support;
                A(1:N,5:N+4) = -1/(8*pi)*(  4*pi*sigma/((sqrt(2*pi))^3)*exp(-dist/sigma^2/2)  + 2/sqrt(pi)* (pi^(1/2)*(sqrt(dist) + (sigma^2)./sqrt(dist)).*erf(sqrt(dist)/(sqrt(2)*sigma)))/2 );
                else if strcmpi(type, 'Wendland')
                        r = sqrt(dist) / support;
                        phi= ((1-r).^4).*(4*r+1);
                        phi(r>=1) = 0;
                        A(1:N,5:N+4) = phi;
                    end
                end
            end
            A(1:N,5:N+4)  = A(1:N,5:N+4)  + lamda*W_inv; %I_f;
            A(N+1,5:N+4) = ones(1, N);
            A(N+2,5:N+4) = x1';
            A(N+3,5:N+4) = y1';
            A(N+4,5:N+4) = z1';
            A(isnan(A)) = 0;

%             transform = inv(A)*target;
            transform = A\target;
        end
        
        function surfout = transformSurfaceTPS(surfin, refsurf, transform, type, support)
            if nargin < 4
                type = 'Regular';
            end
            if nargin < 5
                support = 1;
            end
            
            len = surfin.numverts;
            N = refsurf.numverts;
            A = zeros(len+4, N+4);
            A(1:len, 1) = ones(len,1);
            A(1:len, 2) = surfin.vertices(:,1);
            A(1:len, 3) = surfin.vertices(:,2);
            A(1:len, 4) = surfin.vertices(:,3);
            points = surfin.vertices;
            ref_points = refsurf.vertices;
            dist = zeros(1, N, 'single');
            % for i = 1:N
            %     dist(:,i) = (points(:,1)-ref_points(i,1)).^2 + (points(:,2)- ref_points(i,2)).^2 +(points(:,3)- ref_points(i,3)).^2;
            % end
            if strcmpi(type, 'Regular')
                for i = 1:len
                    curpt = points(i,:);
                    dist = sum((repmat(curpt, [length(ref_points), 1]) - ref_points).^2, 2);
            %         for j = 1:N
            %             dist(j) = (points(i,1)-ref_points(j,1)).^2 + (points(j,2)- ref_points(j,2)).^2 +(points(i,3)- ref_points(j,3)).^2;
            %         end
                    dist(find(dist==0)) = 1;
                    A(i,5:N+4) = dist.*log(dist);
                % A(1:len,4:N+3) = -1/(8*pi)*sqrt(dist);
                end
            else if strcmpi(type, 'Gaussian')
                    for i = 1:len
                        curpt = points(i,:);
                        dist = sum((repmat(curpt, [length(ref_points), 1]) - ref_points).^2, 2);
                        dist(find(dist==0)) = 0.001;
                        sigma = support;
                        A(i,5:N+4) = -1/(8*pi)*(  4*pi*sigma/((sqrt(2*pi))^3)*exp(-dist/sigma^2/2)  + 2/sqrt(pi)* (pi^(1/2)*(sqrt(dist) + (sigma^2)./sqrt(dist)).*erf(sqrt(dist)/(sqrt(2)*sigma)))/2 );
                    end
                    else if strcmpi(type, 'Wendland')
                            for i = 1:len
                                curpt = points(i,:);
                                dist = sum((repmat(curpt, [length(ref_points), 1]) - ref_points).^2, 2);
                                r = sqrt(dist) / support;
                                phi = ((1-r).^4).*(4*r+1);
                                phi(r>=1) = 0;
                                A(i,5:N+4) = phi;
                            end
                        end
                end
            end
            A(len+1,5:N+4) = ones(1, N);
            A(len+2,5:N+4) = ref_points(:,1)';  % control points
            A(len+3,5:N+4) = ref_points(:,2)';
            A(len+4,5:N+4) = ref_points(:,3)';
            A(isnan(A)) = 0;

            new_points = A*transform;
            surfout = surfin;
            surfout.vertices = new_points(1:len, :);            
        end
        
    end % methods
end % surfacetools
