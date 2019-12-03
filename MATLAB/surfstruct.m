classdef surfstruct
    properties
        vertices = 0;
        triangles = 0;
        dim = 0;
        vox = 0;
        orient = 0;
        color = 0;
        id = 0;
        numverts = 0;
        numtris = 0;
        isMm = 0;
        vertexNormals = 0;
        vertexNtriangles = 0;
        triangleNormals = 0;
        vertexProfiles = 0;
        deactivateFlag = 0;
        
    end
    methods
        function obj = surfstruct(vertices, triangles, dim, vox, orient, color)
            obj.id = 0;
            obj.isMm = 0;
            obj.vertexNormals = 0;
            obj.vertexNtriangles = 0;
            obj.triangleNormals = 0;
            obj.vertexProfiles = 0;
            obj.deactivateFlag = 0;
            if nargin < 6
                obj.color = [255,0,0];
            end
            if nargin < 5
                orient = [1,2,3];
                color = [255,0,0];
            end
            if nargin < 3
                dim = [100, 100, 100];
                vox = [1,1,1];
            end
            if nargin == 0
                vertices = zeros(1,3);
                triangles = zeros(1,3);
                dim = [0, 0, 0];
                vox = [0, 0, 0];
                orient = [0, 0, 0];
                color = [0, 0, 0];
            end
            obj.vertices = vertices;
            obj.triangles = triangles;
            obj.dim = dim;
            obj.vox = vox;
            obj.orient = orient;
            obj.color = color;
            obj.numverts = length(obj.vertices);
            obj.numtris = length(obj.triangles);
        end
    end
end

