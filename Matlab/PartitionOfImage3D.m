classdef PartitionOfImage3D < handle
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
        NZ
        NRPI
        NCPI
        NZPI
        Mode
        W
    end
    methods
        function PI3D = PartitionOfImage3D(image,patchSize,margin,mode)
% Example:
%
% load mri
% V = double(squeeze(D))/255;
% PI3D = PartitionOfImage3D(V,32,8,'accumulate');
% PI3D.createOutput(1);
% for i = 1:PI3D.NumPatches
%     disp([i PI3D.NumPatches])
%     P = PI3D.getPatch(i);
%     Q = imgaussfilt3(P,2);
%     PI3D.patchOutput(i,Q);
% end
% W = normalize(PI3D.getValidOutput());
% X = normalize(imgaussfilt3(V,2));
% z = round(size(V,3)/2);
% imshow([V(:,:,z) W(:,:,z) X(:,:,z)])
%
% Marcelo Cicconet, Jan 23 2019
            
            PI3D.Image = image;
            PI3D.PatchSize = patchSize;
            PI3D.Margin = margin;
            subPatchSize = patchSize-2*margin;
            PI3D.SubPatchSize = subPatchSize;

            W = ones(patchSize,patchSize,patchSize);
            W([1,end],:,:) = 0;
            W(:,[1,end],:) = 0;
            W(:,:,[1,end]) = 0;
            for i = 1:2*margin-1
                v = i/(2*margin);
                W([i+1,patchSize-i],i+1:patchSize+1-i-1,i+1:patchSize+1-i-1) = v;
                W(i+1:patchSize+1-i-1,[i+1,patchSize-i],i+1:patchSize+1-i-1) = v;
                W(i+1:patchSize+1-i-1,i+1:patchSize+1-i-1,[i+1,patchSize-i]) = v;
            end
            PI3D.W = W;

            if length(size(image)) == 3
                [nr,nc,nz] = size(image);
            elseif length(size(image)) == 4 % multi-channel image
                [nr,nc,nz,nw] = size(image); % rows, cols, planes, channels
            end
            PI3D.NR = nr;
            PI3D.NC = nc;
            PI3D.NZ = nz;

            npr = ceil(nr/subPatchSize); % number of patch rows
            npc = ceil(nc/subPatchSize); % number of patch cols
            npz = ceil(nz/subPatchSize); % number of patch planes

            nrpi = npr*subPatchSize+2*margin; % number of rows in padded image 
            ncpi = npc*subPatchSize+2*margin; % number of cols in padded image 
            nzpi = npz*subPatchSize+2*margin; % number of plns in padded image 

            PI3D.NRPI = nrpi;
            PI3D.NCPI = ncpi;
            PI3D.NZPI = nzpi;

            if length(size(image)) == 3
                PI3D.PaddedImage = zeros(nrpi,ncpi,nzpi);
                PI3D.PaddedImage(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz) = image;
            elseif length(size(image)) == 4
                PI3D.PaddedImage = zeros(nrpi,ncpi,nzpi,nw);
                PI3D.PaddedImage(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz,:) = image;
            end
 
            PI3D.PC = []; % patch coordinates [r0,r1,c0,c1,z0,z1]
            for i = 0:npr-1
                r0 = i*subPatchSize+1;
                r1 = r0+patchSize-1;
                for j = 0:npc-1
                    c0 = j*subPatchSize+1;
                    c1 = c0+patchSize-1;
                    for iZ = 0:npz-1
                        z0 = iZ*subPatchSize+1;
                        z1 = z0+patchSize-1;
                        PI3D.PC = [PI3D.PC; [r0,r1,c0,c1,z0,z1]];
                    end
                end
            end

            PI3D.NumPatches = size(PI3D.PC,1);
            PI3D.Mode = mode; % 'replace' or 'accumulate'
        end
        
        function P = getPatch(PI3D,i)
            r0 = PI3D.PC(i,1); r1 = PI3D.PC(i,2);
            c0 = PI3D.PC(i,3); c1 = PI3D.PC(i,4);
            z0 = PI3D.PC(i,5); z1 = PI3D.PC(i,6);
            
            if length(size(PI3D.PaddedImage)) == 3
                P = PI3D.PaddedImage(r0:r1,c0:c1,z0:z1);
            elseif length(size(PI3D.PaddedImage)) == 4
                P = PI3D.PaddedImage(r0:r1,c0:c1,z0:z1,:);
            end
        end
        
        function createOutput(PI3D,nChannels)
            if nChannels == 1
                PI3D.Output = zeros(PI3D.NRPI,PI3D.NCPI,PI3D.NZPI);
            else
                PI3D.Output = zeros(PI3D.NRPI,PI3D.NCPI,PI3D.NZPI,nChannels);
            end
            if strcmp(PI3D.Mode,'accumulate')
                PI3D.Count = zeros(PI3D.NRPI,PI3D.NCPI,PI3D.NZPI);
            end
        end
        
        function patchOutput(PI3D,i,P)
            r0 = PI3D.PC(i,1); r1 = PI3D.PC(i,2);
            c0 = PI3D.PC(i,3); c1 = PI3D.PC(i,4);
            z0 = PI3D.PC(i,5); z1 = PI3D.PC(i,6);
            
            if strcmp(PI3D.Mode,'accumulate')
                PI3D.Count(r0:r1,c0:c1,z0:z1) = PI3D.Count(r0:r1,c0:c1,z0:z1)+PI3D.W;
            end
            if length(size(P)) == 3
                if strcmp(PI3D.Mode,'accumulate')
                    PI3D.Output(r0:r1,c0:c1,z0:z1) = PI3D.Output(r0:r1,c0:c1,z0:z1)+P.*PI3D.W;
                elseif strcmp(PI3D.Mode,'replace')
                    PI3D.Output(r0:r1,c0:c1,z0:z1) = P;
                end
            elseif length(size(P)) == 4
                if strcmp(PI3D.Mode,'accumulate')
                    for i = 1:size(P,4)
                        PI3D.Output(r0:r1,c0:c1,z0:z1,i) = PI3D.Output(r0:r1,c0:c1,z0:z1,i)+P(:,:,:,i).*PI3D.W;
                    end
                elseif strcmp(PI3D.Mode,'replace')
                     PI3D.Output(r0:r1,c0:c1,z0:z1,:) = P;
                end
            end
        end
        
        function R = getValidOutput(PI3D)
            margin = PI3D.Margin;
            nr = PI3D.NR; nc = PI3D.NC; nz = PI3D.NZ;
            if strcmp(PI3D.Mode,'accumulate')
                C = PI3D.Count(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz);
            end
            if length(size(PI3D.Output)) == 3
                if strcmp(PI3D.Mode,'accumulate')
                    R = PI3D.Output(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz)./C;
                elseif strcmp(PI3D.Mode,'replace')
                    R = PI3D.Output(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz);
                end
            elseif length(size(PI3D.Output)) == 4
                if strcmp(PI3D.Mode,'accumulate')
                    for i = 1:size(PI3D.Output,4)
                        PI3D.Output(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz,i) = PI3D.Output(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz,i)./C;
                    end
                end
                R = PI3D.Output(margin+1:margin+nr,margin+1:margin+nc,margin+1:margin+nz,:);
            end
        end
    end
end