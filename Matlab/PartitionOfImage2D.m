classdef PartitionOfImage2D < handle
    properties
        Image
        PaddedImage
        PatchSize
        Margin
        SubPatchSize
        PC
        NumPatches
        Output
        Count
        NR
        NC
        NRPI
        NCPI
        Mode
        W
    end
    methods
        function PI2D = PartitionOfImage2D(image,patchSize,margin,mode)
% Example:
%
% I = im2double(imread('rice.png'));
% PI2D = PartitionOfImage2D(I,32,8,'accumulate');
% PI2D.createOutput(1);
% for i = 1:PI2D.NumPatches
%     disp([i PI2D.NumPatches])
%     P = PI2D.getPatch(i);
%     Q = imgaussfilt(P,2);
%     PI2D.patchOutput(i,Q);
% end
% J = normalize(PI2D.getValidOutput());
% K = normalize(imgaussfilt(I,2));
% imshow([I J K])
%
% Marcelo Cicconet, Jan 23 2019
            
            PI2D.Image = image;
            PI2D.PatchSize = patchSize;
            PI2D.Margin = margin;
            subPatchSize = patchSize-2*margin;
            PI2D.SubPatchSize = subPatchSize;

            W = ones(patchSize,patchSize);
            W([1,end],:) = 0;
            W(:,[1,end]) = 0;
            for i = 1:2*margin-1
                v = i/(2*margin);
                W([i+1,patchSize-i],i+1:patchSize+1-i-1) = v;
                W(i+1:patchSize+1-i-1,[i+1,patchSize-i]) = v;
            end
            PI2D.W = W;

            if length(size(image)) == 2
                [nr,nc] = size(image);
            elseif length(size(image)) == 3 % multi-channel image
                [nr,nc,nz] = size(image); % rows, cols, channels
            end
            PI2D.NR = nr;
            PI2D.NC = nc;

            npr = ceil(nr/subPatchSize); % number of patch rows
            npc = ceil(nc/subPatchSize); % number of patch cols

            nrpi = npr*subPatchSize+2*margin; % number of rows in padded image 
            ncpi = npc*subPatchSize+2*margin; % number of cols in padded image 

            PI2D.NRPI = nrpi;
            PI2D.NCPI = ncpi;

            if length(size(image)) == 2
                PI2D.PaddedImage = zeros(nrpi,ncpi);
                PI2D.PaddedImage(margin+1:margin+nr,margin+1:margin+nc) = image;
            elseif length(size(image)) == 3
                PI2D.PaddedImage = zeros(nrpi,ncpi,nz);
                PI2D.PaddedImage(margin+1:margin+nr,margin+1:margin+nc,:) = image;
            end
 
            PI2D.PC = []; % patch coordinates [r0,r1,c0,c1,z0,z1]
            for i = 0:npr-1
                r0 = i*subPatchSize+1;
                r1 = r0+patchSize-1;
                for j = 0:npc-1
                    c0 = j*subPatchSize+1;
                    c1 = c0+patchSize-1;
                    PI2D.PC = [PI2D.PC; [r0,r1,c0,c1]];
                end
            end

            PI2D.NumPatches = size(PI2D.PC,1);
            PI2D.Mode = mode; % 'replace' or 'accumulate'
        end
        
        function P = getPatch(PI2D,i)
            r0 = PI2D.PC(i,1); r1 = PI2D.PC(i,2);
            c0 = PI2D.PC(i,3); c1 = PI2D.PC(i,4);
            
            if length(size(PI2D.PaddedImage)) == 2
                P = PI2D.PaddedImage(r0:r1,c0:c1);
            elseif length(size(PI2D.PaddedImage)) == 3
                P = PI2D.PaddedImage(r0:r1,c0:c1,:);
            end
        end
        
        function createOutput(PI2D,nChannels)
            if nChannels == 1
                PI2D.Output = zeros(PI2D.NRPI,PI2D.NCPI);
            else
                PI2D.Output = zeros(PI2D.NRPI,PI2D.NCPI,nChannels);
            end
            if strcmp(PI2D.Mode,'accumulate')
                PI2D.Count = zeros(PI2D.NRPI,PI2D.NCPI);
            end
        end
        
        function patchOutput(PI2D,i,P)
            r0 = PI2D.PC(i,1); r1 = PI2D.PC(i,2);
            c0 = PI2D.PC(i,3); c1 = PI2D.PC(i,4);
            
            if strcmp(PI2D.Mode,'accumulate')
                PI2D.Count(r0:r1,c0:c1) = PI2D.Count(r0:r1,c0:c1)+PI2D.W;
            end
            if length(size(P)) == 2
                if strcmp(PI2D.Mode,'accumulate')
                    PI2D.Output(r0:r1,c0:c1) = PI2D.Output(r0:r1,c0:c1)+P.*PI2D.W;
                elseif strcmp(PI2D.Mode,'replace')
                    PI2D.Output(r0:r1,c0:c1) = P;
                end
            elseif length(size(P)) == 3
                if strcmp(PI2D.Mode,'accumulate')
                    for i = 1:size(P,3)
                        PI2D.Output(r0:r1,c0:c1,i) = PI2D.Output(r0:r1,c0:c1,i)+P(:,:,i).*PI2D.W;
                    end
                elseif strcmp(PI2D.Mode,'replace')
                     PI2D.Output(r0:r1,c0:c1,:) = P;
                end
            end
        end
        
        function R = getValidOutput(PI2D)
            margin = PI2D.Margin;
            nr = PI2D.NR; nc = PI2D.NC;
            if strcmp(PI2D.Mode,'accumulate')
                C = PI2D.Count(margin+1:margin+nr,margin+1:margin+nc);
            end
            if length(size(PI2D.Output)) == 2
                if strcmp(PI2D.Mode,'accumulate')
                    R = PI2D.Output(margin+1:margin+nr,margin+1:margin+nc)./C;
                elseif strcmp(PI2D.Mode,'replace')
                    R = PI2D.Output(margin+1:margin+nr,margin+1:margin+nc);
                end
            elseif length(size(PI2D.Output)) == 3
                if strcmp(PI2D.Mode,'accumulate')
                    for i = 1:size(PI2D.Output,3)
                        PI2D.Output(margin+1:margin+nr,margin+1:margin+nc,i) = PI2D.Output(margin+1:margin+nr,margin+1:margin+nc,i)./C;
                    end
                end
                R = PI2D.Output(margin+1:margin+nr,margin+1:margin+nc,:);
            end
        end
    end
end