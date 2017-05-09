function writeROI(inCoords,inData,varargin)
% This is a general function that will write nifti files from given data
% and its associated voxel coordinates. Its utility is in several different
% thresholding/manipulation features. The necessary arguments are a vector
% of coordinates (inCoords) and a vector of data (inData). The vector of
% data may be left empty if a vector of p-values is supplied to write out
% instead be used instead (see below). Additional arguments can be added in
% the following order: 
% 
%% Required Inputs:
% inCoords : this is your vector of coordinates to write out in a nifti
% file. Each row corresponds to a unique voxel. Each voxel should
% correspond to its data point in inData (i.e. rows must match). If the
% vector has a single column, the script will assume that the values of
% each row correspond to indices of a 3D matrix, specifically of some
% template that you provide (or you can default to 2mm/1mm space; see
% below). Alternatively, you can feed in a 3 column vector where the first
% column will be considered the x-coordinate (in voxel space NOT mm), the
% second column will be the y-coordinate, and the third column will be the
% z-coordinate.
%
% inData : this is your vector of data. Each row is a voxel, and should
% match the voxel identity (i.e. row) in inCoords. Columns are treated as
% seperate maps that you want to write. 
%
%% Optional Inputs:
% outPaths : full directory to write data to, plus whatever base name you
% want to give your file. If you have multiple columns in inData, you can
% provide multiple filepaths to write out. If you only have a single
% column, this variable can be a string, otherwise each path should be a
% row within a cell. If you do not provide as many file names as you have
% columns in inData, then the last file name you provided will be used for
% all of the missing names. Files are appended based on the type of
% thresholding that you chose. If this variable is left empty, your files
% will be placed in the working directory, but if you don't distinguish
% multiple input files (i.e. inData cols), your data will be overwritten. 
%
% templatePath : this either a path to your template file, or the already
% imported nifti structure (using load_untoch_nii). The size of the image,
% and the nifti header are used to write out your new files (i.e. your
% voxels MUST correspond to this image dimensions/resolution). If this is
% left empty, a 2mm MNI template will be used. If this is set to 1, then a
% 1mm MNI template will be used. 
%
% adjustedSwitch : if this is set to 'true', then coordinates will be
% adjusted for matlab indexing. In nifti, voxel space indexing starts at
% zero, but this corresponds to row 1 in matlab. This will essentially just
% add 1 to your coordinates. Otherwise, set to 'false'.
%
% threshVal : if this is not empty, your inData will be thresholded by this
% given value. Otherwise, the script will simply write out all of inData.
% You can specify whether you'd like to threshold only things above this
% value, or everything between this value and its inverse (see threshTail). 
%
% pVals : this is a vector identical to inData in size that provides p -
% values for each observed row/col combination. If this is provided along
% with inData, the threshold will be applied to these p-values, but the the
% inData values of the voxels surviving the threshold will be written into
% the output nifti file. Alternatively, inData can be left empty, in which
% case the p-values themselves will be written out (threshold will still be
% applied if specified).
%
% threshTail : this determines whether only values above the threshold will
% be considered, or whether you don't want to consider anything lying between
% the threshold and its inverse. If this is set to 'two', then any value
% above the supplied threshold, or below the reciprocal of that threshold
% will be removed. If this is set to 'one' only values above the threshold
% will be considered.
%
% top : This will remove all voxels not in the top X% of inData. If
% top.Switch is set to 'true', then top.Percent is a whole number that
% corresponds to the top X% of inData that you want to keep. Otherwise, set
% to 'false'.
%
% Examples:
% writeROI(inCoords,inData,outPaths,templatePath,adjustedSwitch,threshVal,pVals,threshTail,topSwitch)
% writeROI([50, 80, 25; 25 25 25; 50 50 50],[3,2,12,55;3,1,4,8;3,10,435,4],{[pwd '/Test1']; [pwd '/Test2']; [pwd '/Test3']},[],'true',1,[],[],'false')
%
%% Alex Notes: 
% - If using default pwd for outFiles, make sure that files aren't being
% overwritten by appending them based on their loop. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Alex Teghipco (alex.teghipco@uci.edu)
% Last update 5/9/17


%% Defaults
outPaths = {pwd};
scriptPath = which('writeROI');
scriptPath = scriptPath(1:end-11);
templatePath = load([scriptPath '/2mmTemplate.mat'],'template');
adjustedSwitch = 'false';
threshVal = [];
pVals = [];
threshTail = 'false';
top.Switch = 'false';

%% Import user-specified settings
fixedargn = 2;
if nargin > (fixedargn + 0)
    if ~isempty(varargin{1})
        outPaths = varargin{1};
    end
end
if nargin > (fixedargn + 1)
    if ~isempty(varargin{2})
        templatePath = varargin{2};
    end
end
if nargin > (fixedargn + 2)
    if ~isempty(varargin{3})
        adjustedSwitch = varargin{3};
    end
end
if nargin > (fixedargn + 3)
    if ~isempty(varargin{4})
        threshVal = varargin{4};
    end
end
if nargin > (fixedargn + 4)
    if ~isempty(varargin{5})
        pVals = varargin{5};
    end
end
if nargin > (fixedargn + 5)
    if ~isempty(varargin{6})
        threshTail = varargin{6};
    end
end
if nargin > (fixedargn + 6)
    if ~isempty(varargin{7})
        top = varargin{7};
    end
end

%% Setup
% If templatePath is a structure we will assume it is in standard nifti
% structure format
if isstruct(templatePath) == 1
    warning('Assuming your template path is actually a preloaded nifti structure')
    template = templatePath;
else
    %if templatePath is 1 then load up a pre-assembled 1mm template
    if templatePath == 1
        warning('Using the 1mm MNI template')
        templatePath = load([scriptPath '/1mmTemplate.mat','template']);
    end
    %if templatePath is a file, then check if it is unzipped.
    if exist(templatePath,'file') == 2
        disp('Your nifti input was unzipped for import ... ')
        templatePath = checkGZ(templatePath,'yes');
        template = load_untouch_nii(templatePath);
    end
end
clear templatePath

%If you give a 2D matrix, assume that it is a reshaped subset of the
%3d image matrix. Get all of the associated coordinate identities. 
if(size(inCoords,2)) == 1
    warning('Assuming you have given me a reshaped 3d  matrix whose linear indices are a subset of the template ...')
    if exist('templatePath','var') == 0
        %niftiMatTemp = zeros([size(template.img,1)*size(template.img,2)*size(template.img,3),1]);
        [inCoordTemp(:,1),inCoordTemp(:,2),inCoordTemp(:,3)] = ind2sub(size(template.img),inCoords);
        inCoords = inCoordTemp;
    end
    clear inDataTemp
end

%If adjustedSwitch is turned on, shift voxels up by 1. 
switch adjustedSwitch
    case 'true'
        warning('Adjusting voxel coordinate position by 1')
        inCoords = inCoords+1;
    case 'fase'
        warning('Make sure that you have already adjusted for zero starting position of voxel space.')
end

%If only one output file name is given and as a string, convert it to a
%cell for easier output file name processing
if ischar(outPaths) == 1 
    outPaths = cellstr(outPaths);
end

%If you intended to write out more files than you gave names, then make all
%missing file names the last known file name. It assumes missing files
%always come after known files.
if isempty(inData) == 0
    files = cell(size(inData,2),1);
    if size(outPaths,1) < size(inData,2)
        lastName = size(outPaths,1);
        missingEnd = size(inData,2);
        for i = lastName:missingEnd
            if i >= lastName
                files{i} = outPaths{lastName};
            end
        end
    else
        files = outPaths;
    end
    clear outPaths
end
if isempty(inData) == 1
    files = cell(size(pVals,2),1);
    if size(outPaths,1) < size(pVals,2)
        lastName = size(outPaths,1);
        missingEnd = size(pVals,2);
        for i = lastName:missingEnd
            if i >= lastName
                files{i} = outPaths{lastName};
            end
        end
    else
        files = outPaths;
    end
    clear outPaths
end

%If top.Switch is turned on, keep only the top X% of the data. 
switch top.Switch
    case 'true'
        disp(['Keeping only the top ' num2str(top.Percent) ' % of data ... ']);
        if isempty(inData) == 0
            for i = 1:size(inData,2)
                sortedReal = sort(inData(:,i));
                topPercentage = round(size(sortedReal,1)*((top.Percent)/100));
                topVal = sortedReal(end-topPercentage);
                idx = find(inData(:,i) < topVal);
                inData(idx,i) = 0;
            end
            for i = 1:size(files,1)
                files{i,1} = [files{i,1} '_Top_' num2str(top.Percent) '_Percent'];
            end
        else
           warning('Input data is empty...if you are trying to keep the top 5% of p-values, this functionality is not supported') %
        end
end

%% Main loop: if inData is not empty, then threshold as specified (seperate loop for p-values)
if isempty(inData) == 0
    for i = 1:size(inData,2)
        niftiMat = zeros(size(template.img,1),size(template.img,2),size(template.img,3)); %generate an empty vector of template size.
        for j = 1:size(inCoords,1)
            if isempty(threshVal) == 0 && isempty(pVals) == 1
                switch threshTail
                    case 'one'
                        if inData(j,i) >= threshVal
                            niftiMat(inCoords(j,1),inCoords(j,2),inCoords(j,3)) = inData(j,i);
                        end
                    case 'two'
                        if abs(inData(j,i)) >= threshVal
                            niftiMat(inCoords(j,1),inCoords(j,2),inCoords(j,3)) = inData(j,i);
                        end
                    case 'false'
                end
            end
            if isempty(threshVal) == 0 && isempty(pVals) == 0
                if pVals(j,i) <= threshVal
                    niftiMat(inCoords(j,1),inCoords(j,2),inCoords(j,3)) = inData(j,i);
                end
            end
            if isempty(threshVal) == 1 && isempty(pVals) == 1
                niftiMat(inCoords(j,1),inCoords(j,2),inCoords(j,3)) = inData(j,i);
            end
        end
        %% append filename based on steps
        if isempty(threshVal) == 0 && isempty(pVals) == 1
            switch threshTail
                case 'one'
                    files{i,1} = [files{i,1} '_OneTailThreshold_' num2str(threshVal)];
                case 'two'
                    files{i,1} = [files{i,1} '_TwoTailThreshold_' num2str(threshVal)];
                case 'false'
            end
        end
        if isempty(threshVal) == 0 && isempty(pVals) == 0
            files{i,1} = [files{i,1} '_PCorrected' num2str(threshVal)];
        end
        %% now write out the file
        template.img = double(niftiMat);
        save_untouch_nii(template, [files{i,1} '.nii']);
    end
end

%% Main loop: if inData is empty, assuming you want to threshold p-values.
if isempty(inData) == 1
    inData = pVals;
    for i = 1:size(inData,2)
        niftiMat = zeros(size(template.img,1),size(template.img,2),size(template.img,3)); %generate an empty vector of template size.
        for j = 1:size(inCoords,1)
            if isempty(threshVal) == 0 
                if inData(j,i) <= threshVal
                    niftiMat(inCoords(j,1),inCoords(j,2),inCoords(j,3)) = inData(j,i);
                end
            end
            if isempty(threshVal) == 1
                niftiMat(inCoords(j,1),inCoords(j,2),inCoords(j,3)) = inData(j,i);
            end
        end
        %% append filename based on steps
        if isempty(threshVal) == 1
            files{i,1} = [files{i,1} '_PValues' num2str(threshVal)];
        end
        if isempty(threshVal) == 0
            files{i,1} = [files{i,1} '_PValues_Thresh_' num2str(threshVal)];
        end
        %% now write out the file
        template.img = double(niftiMat);
        save_untouch_nii(template, [files{i,1} '.nii']);
    end
end

