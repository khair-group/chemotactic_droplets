classdef multidroplets
    %object for the analysis of multicomponent particles samples where the
    %particles can be distinguished by grayscale intensity
    %
    %Requires helper codes:
    % - get_numerical_input.m
    % - get_textual_input.m
    % - track.m
    % - bpass.m
    % - pkfind.m
    %
    %Credits:
    % - track written by J.C.Crocker in 1993 and translated to matlab by D.Blair in 2005
    % - pkfind written by Eric Dufresne in 2005
    % - bpass written by David Grier and John Crocker in 1993. Translated to Matlab by D.Blair.
    %       Copyright 1997 D.Grier & John Crocker. Contributions by E.R.D & J.W.M
    %       
    %
    %Shorthand user manual
    %1) Create object of type multidroplets: objname=multidroplets
    %2) Find image processing parameters: objname=optimizeimprocparams(objname)
    %       -Repeat to save multiple sets of parameters
    %3) Find particle detection parameters: objname=optimizefindparams(objname)
    %       -Repeat for each set of parameters from (2)
    %4) Process image & detect particles: objname=find_pos(objname)
    %       -Repeat for each set of parameters from (2)
    %5) Remove duplicate particle detections: objname=remove_overlap_pos(objname)
    %6) Split  positions by particle type based on a threshold grayscale
    %value: objname=split_types(objname,threshold,R);
    %7)  Track particles: objname=track_obj_split(objname)
    %       -Repeat to track both sets of the split positions.
    %
    %
    %Modification history
    %June2023: created by Pepijn Moerman
    properties
        filename    %Name of the file that this object analyzes
        foldername  %Name of the folder in which that file is located
        NOF         %Number of frames
        improcparams%List of image processing parameter sets {set 1, set 2..}
        findparams  %List of particle finding parameter sets {set 1, set 2..}
        trparams    %Tracking parameter sets for pos_split1 and post_split2
        splitparams %Threshold and R used for distinguishing bewteen the two types of particles: {T R}
        pos_all     %List of positions for each parameter set: {set1frame1[x y] set2frame1[x y]; set1frame2[x y]..}} .
        pos_clean   %List of all positions that do not overlap [x y frame]
        pos_split   %List of positions split by type: {type 1 [x y frame], type 2 [x y frame]}.
        tr_split    %Tracks of each particle, split by type: {type 1 [x y frame ID], type 2 [x y frame ID]}
    end
    methods
        function obj=multidroplets(fname)
            %connect the object to a tiffstack. Save file and foldername
            if nargin==0
                disp('Select folder that contains the tiffstack(s)')
                [filen, foldername_temp] = uigetfile('*.tif','Select tif stack');
                fname=[foldername_temp filen];
            else
                startIndex=max(regexp(fname,'\'));
                foldername_temp=fname(1:startIndex);
                filen=fname(startIndex+1:end);
            end
            
            %Open the first image; save NaN if the file cannot be opened
            try A = imread(fname,1);
            catch
                disp(['No tifstack exists with name: ' foldername_temp filen '.'])
                disp('The imagestack is now "NaN". Be aware that some functions might not work')
                foldername_temp=NaN;
                filen=NaN;
            end
            
            %Save the number of frames in the stack
            info = imfinfo(fname);
            obj.NOF = numel(info);
            
            %Save the filename, and foldername
            obj.foldername=foldername_temp;
            if isnan(filen)
                obj.filename=NaN;
            else
                obj.filename=regexprep(filen,'.tif','');
            end
            
            %Prep elements of this object that have predetermined stucture
            obj.pos_all={};
        end
        function obj=optimize_improcparams(obj,i)
            %This function uses the ith image in the stack to interactively
            %select image preprocessing parameters. The parameters are
            %saved in a cell of type
            %{1. "Invert Contrast"  'y/n'                 }
            %{2. "Binarize"         'y/n'   treshold      }
            %{3. "Bandpass"         'y/n'   high       low}
            %{4. "Gaussblur"        'y/n'   length        }
            %The processing parameters can be saved to set 1 or set 2;
            %Set_n determines which set to save to. Can be either 1 or 2.
            %i determines the frame on which the image analysis is used
            
            if nargin<2
                i=1; %If not specified, the 1st image is used to find image preprocessing parameters
            else
                if i>obj.NOF
                    disp(['The tiffstack does not have a frame ' num2str(i) '. Selecting frame 1 instead.'])
                    i=1;
                end
            end
            
            %read the selected image in the tiffstack
            fname=[obj.foldername obj.filename '.tif'];
            A = imread(fname, 1);
            
            %Assert image is 8-bit
            if length(size(A))>2
                disp('Image is not of type 8-bit. Converting to 8bit first.')
                %Convert to grayscale
                A = rgb2gray(A);
            end
            
            %Stretch contrast to go from 1-255;
            A = uint8((double(A)-double(min(A(:))))/(double(max(A(:)))-double(min(A(:))))*255);
            size_A=size(A);
            
            % Check: Invert image contrast
            message= 'The tracking code finds light objects in a dark background. Do you want to invert the contrast (y/n): ' ;
            acceptables=['y' 'n'];
            figure
            hold on
            colormap('gray'),imagesc(A);
            axis([0 size_A(2) 0 size_A(1)])
            hold off
            invert_check = get_textual_input(message,acceptables);
            if invert_check=='y'
                A=255-A;
            end
            
            % Check: apply bandpass filter;
            %Select a high and low bandpass filter value
            message='Do you want to apply a bandpass filter (y/n): ';
            acceptables=['y' 'n'];
            bpass_check=get_textual_input(message,acceptables);
            lowrange=1;
            highrange=2;
            
            ready='n';
            while bpass_check=='y' && ready~='y'
                Ai = bpass(A,lowrange,highrange);
                figure
                hold on
                colormap('gray'),imagesc(Ai);
                axis([0 size_A(2) 0 size_A(1)])
                hold off
                message='Do you want to continue with these values? If so: press "y". If not: press "l" to change the low value filter or "h" to change the high value filter: ';
                acceptables=['y' 'l' 'h'];
                ready=get_textual_input(message,acceptables);
                if ready=='l'
                    message= 'What new value for the low range filter you want to use? ' ;
                    lowrange= get_numerical_input(message);
                elseif ready=='h'
                    message= 'What new value for the high range filter you want to use? ' ;
                    highrange= get_numerical_input(message);
                end
            end
            if bpass_check=='y'
                A=Ai;
            end
            
            % Check: binarize image;
            %Select a threshold value for binarization
            message='Do you want to binarize the image (y/n): ';
            acceptables=['y' 'n'];
            binarize_check=get_textual_input(message,acceptables);
            if binarize_check=='y'
                A = imbinarize(A);
                A=double(A);
                figure
                hold on
                colormap('gray'),imagesc(A);
                axis([0 size_A(2) 0 size_A(1)])
                hold off
            end
            
            % Check: gaussian blur;
            %Select a high and low bandpass filter value
            message=('Do you want apply gaussian blur (y/n): ');
            acceptables=['y' 'n'];
            blur_check=get_textual_input(message,acceptables);
            accept='n';
            blur_length=1;
            while blur_check=='y' && accept~='y'
                Ai=imgaussfilt(A,blur_length);
                figure
                hold on
                colormap('gray'),imagesc(Ai);
                axis([0 size_A(2) 0 size_A(1)])
                hold off
                message='To proceed press "y", to change the blur lenght, press "n": ';
                acceptables=['y' 'n'];
                accept=get_textual_input(message,acceptables);
                if accept=='n'
                    disp(['The current blur length is ' num2str(blur_length) '.'])
                    message= 'What new value for the guassian filter do you want to use? ' ;
                    blur_length= get_numerical_input(message);
                end
            end
            if blur_check=='y'
                A=Ai;
            end
            
            %Save the processing parameters (by overriding or adding to set)
            n_sets=length(obj.improcparams);
            if n_sets>0
                message='Do you want to override a set (y/n): ';
                acceptables=['y' 'n'];
                overridecheck=get_textual_input(message,acceptables);
                if overridecheck=='n'
                    obj.improcparams{n_sets+1}={'1.' 'Invert Contrast' invert_check [] [];
                        '2.' 'Bandpass' bpass_check lowrange highrange;
                        '3.' 'Binarize' binarize_check [] [];
                        '4.' 'Blur'     blur_check blur_length []};
                elseif overridecheck=='y'
                    if n_sets==1
                        obj.improcparams{1}={'1.' 'Invert Contrast' invert_check [] [];
                            '2.' 'Bandpass' bpass_check lowrange highrange;
                            '3.' 'Binarize' binarize_check [] [];
                            '4.' 'Blur'     blur_check blur_length []};
                    else
                        message=['Which set do you want to override (max ' num2str(n_sets) '): '];
                        overriden=get_numerical_input(message);
                        while overriden>n_sets
                            message=['That set does not exist and cannot be overriden. Select a value equal to or less than ' num2str(n_sets) '): '];
                            overriden=get_numerical_input(message);
                        end
                        obj.improcparams{overriden}={'1.' 'Invert Contrast' invert_check [] [];
                            '2.' 'Bandpass' bpass_check lowrange highrange;
                            '3.' 'Binarize' binarize_check [] [];
                            '4.' 'Blur'     blur_check blur_length []};
                    end
                end
            else
                obj.improcparams{1}={'1.' 'Invert Contrast' invert_check [] [];
                    '2.' 'Bandpass' bpass_check lowrange highrange;
                    '3.' 'Binarize' binarize_check [] [];
                    '4.' 'Blur'     blur_check blur_length []};
            end
            close all
        end
        function im=preprocess_image(obj,set_n,i)
            %Takes image i from the stack and preprocesses it according to
            %the parameters stored in set_n
            
            if nargin<3
                error('not enough input parameters')
            end
            if set_n>length(obj.improcparams)
                error('there are no image processing parameters for this set number. Select a smaller value');
            end
            fname=[obj.foldername obj.filename '.tif'];
            A=imread(fname,i);
            processparams=obj.improcparams{set_n};
            
            %Assert image is 8-bit
            if length(size(A))>2
                %Convert to grayscale
                A = rgb2gray(A);
            end
            %Stretch contrast to go from 1-255;
            A = uint8((double(A)-double(min(A(:))))/(double(max(A(:)))-double(min(A(:))))*255);
            %Invert contrast
            if processparams{1,3}=='y'
                A=255-A;
            end
            %Bandpass
            if processparams{2,3}=='y'
                lowrange=processparams{2,4};
                highrange=processparams{2,5};
                A = bpass(A,lowrange,highrange);
            end
            %Binarize
            if processparams{3,3}=='y'
                A = imbinarize(A);
            end
            %Blur
            if processparams{4,3}=='y'
                blur_length=processparams{4,4};
                A=imgaussfilt(A,blur_length);
            end
            im=A;
        end
        function obj=optimize_findparams(obj,set_n,i)
            %Code used to select a set of parameters to find particles in your tiffstack
            %saved in a cell of type
            %{1. "Hough transform"  'y/n'   cutoff_low   cutoff_high   senstivity}
            %{2. "Gausspeakfind"    'y/n'   featuresize  sensitivity             }
            %The particle finding parameters can be saved to set 1 or set 2;
            %Set_n determines which set to save to. Can be either 1 or 2.
            %i determines the frame on which the image analysis is used
            %IMPORTANT: particle find params set 1 go with
            %imagepreprocessing parameters set 1 and idem dito for set 2
            if nargin<3
                i=1; %If not specified, the 1st image is used to find image preprocessing parameters
            else
                if i>obj.NOF
                    disp(['The tiffstack does not have a frame ' num2str(i) '. Selecting frame 1 instead.'])
                    i=1;
                end
            end
            if nargin<2
                message='Which preprocessing set do you want to use (also the set number you save to): ';
                set_n=get_numerical_input(message);
            end
            if length(obj.improcparams)<1
                error('Optimize image processing parameters first by running optimize_improcparams');
            end
            while set_n>length(obj.improcparams)
                message=['You only have ' num2str(length(obj.improcparams)) 'prepocessing parameter sets. Select a number equal or smaller: '];
                set_n=get_numerical_input(message);
            end
            
            %Choose between hough transform and gaussian peakfind
            message= 'Choose between hough transform (h) and gaussian peakfind (g): ' ;
            acceptables=['h' 'g'];
            find_type = get_textual_input(message,acceptables);
            
            %Preprocess image i according to the chose set of image
            %processing parameters
            A=preprocess_image(obj,set_n,i);
            size_A=size(A);
            
            %optimize particle finding parameters based on selection
            if find_type=='h'
                hough_check='y';
                gauss_check='n';
                % Set default parameters
                sensitivity=0.9;
                r_lower=1;
                r_upper=10;
                threshold=1;
                feature_size=10;
                
                %identify particles in the image
                accept='n';
                centers=[];
                radii=[];
                while accept ~= 'y'
                    [centers,radii] = imfindcircles(A,[round(r_lower),round(r_upper)],'Sensitivity',sensitivity);
                    if isempty(centers)
                        pkl=0;
                    else
                        pkl=length(centers(:,1));
                    end
                    %show the identified particles
                    figure
                    hold on
                    colormap('gray'),imagesc(A);
                    axis([0 size_A(2) 0 size_A(1)]);
                    if ~isempty(centers)
                        viscircles(centers(:,1:2),radii(:,1));
                    end
                    hold off
                    
                    message='To continue with these parameters press "y". Else, press "s" to change sensitivity, press "l" to change lower bound or press "u" to change upper bound: ';
                    acceptables=['y' 's' 'l' 'u'];
                    accept=get_textual_input(message,acceptables);
                    if accept=='s'
                        message= 'What new value for sensitivity do you want to use? High is sensitive. Value between 0 and 1: ' ;
                        sensitivity= get_numerical_input(message);
                        if sensitivity>1
                            sensitivity=1;
                        elseif sensitivity<0
                            sensitivity=0;
                        end
                    elseif accept=='l'
                        message= 'What new value for the lower bound of radii you want to use? ' ;
                        r_lower= get_numerical_input(message);
                        if round(r_upper)<round(r_lower)
                            r_lower=r_upper-1;
                            if r_lower<1
                                r_lower=1;
                                r_upper=r_lower+1;
                            end
                            disp(['Your value for r_lower was lower than r_upper. We used ' num2str(round(r_lower)) ' instead.'])
                        end
                    elseif accept=='u'
                        message= 'What new value for the upper bound of radii you want use? ' ;
                        r_upper= get_numerical_input(message);
                        if round(r_upper)<round(r_lower)
                            r_upper=r_lower+1;
                            disp(['Your value for r_upper was lower than r_lower. We used ' num2str(round(r_upper)) ' instead.'])
                        end
                    end
                end
                
            elseif find_type=='g'
                gauss_check='y';
                hough_check='n';
                % Take processed image and apply peakfind
                sensitivity=0.9;
                r_lower=1;
                r_upper=10;
                threshold=1;
                feature_size=10;
                accept='n';
                while accept ~= 'y'
                    pk = pkfnd(A,threshold,feature_size);
                    if ~isempty(pk)
                        pk = cntrd(A,pk,2*feature_size+1);
                    end
                    
                    if isempty(pk)
                        pkl=0;
                    else
                        pkl=length(pk(:,1));
                    end
                    radiplot=zeros(pkl,1)+feature_size/2;
                    figure
                    hold on
                    colormap('gray'),imagesc(A);
                    if ~isempty(pk)
                        viscircles(pk(:,1:2),radiplot);
                    end
                    axis([0 size_A(2) 0 size_A(1)]);
                    hold off
                    
                    message='To continue press "y". Else, press "f" to change the feature size or "t" to change threshold: ';
                    acceptables=['y' 'f' 't'];
                    accept=get_textual_input(message,acceptables);
                    if accept=='f'
                        message= 'What feature size value do you want: ' ;
                        feature_size= get_numerical_input(message);
                    elseif accept=='t'
                        message= 'What threshold do you want: ' ;
                        threshold= get_numerical_input(message);
                    end
                end
            end
            %save the parameters
            obj.findparams{set_n}={'1.' "Hough transform"  hough_check   r_lower   r_upper   sensitivity;
                '2.' "Gausspeakfind"    gauss_check   feature_size  threshold []};
            close all
        end
        function obj=find_pos(obj,set_n)
            %finds positions according to the image pre processing parameters and find
            %parameters that have been optimzed before and stored in set_n.
            %saves the positions in element n of the pos cell.
            if nargin<2
                sizen=size(obj.pos_all);
                set_n=sizen(2)+1; %Default is the set that hasn't been chosen yet
                if set_n>length(obj.improcparams)
                    set_n=1;%If everything has already been analyzed, just analyze set 1
                end
                disp(['Analyzing with parameter set ' num2str(set_n) '.']);
            end
            if length(obj.findparams)<1
                error('Optimize parameters first by running optimize_find_params');
            end
            
            fname=[obj.foldername obj.filename '.tif'];
            findparameters=obj.findparams{set_n};
            if findparameters{1,3}=='y' %then use hough transform
                r_lower=findparameters{1,4};
                r_upper=findparameters{1,5};
                sensitivity=findparameters{1,6};
                for i=1:obj.NOF
                    A=imread(fname,i);
                    A=preprocess_image(obj,set_n,i);
                    pos = imfindcircles(A,[round(r_lower),round(r_upper)],'Sensitivity',sensitivity);
                    if ~isempty(pos)
                        pos(find(pos(:,1)==0),:)=[]; %removes particles at origin that are tracking artifact
                    else
                       disp(['No particles were found in frame ' num2str(i) '.'])
                    end
                    obj.pos_all{i,set_n}=pos;
                end
            elseif findparameters{2,3}=='y' %then use gaussianpeakfind
                feature_size=findparameters{2,4};
                threshold=findparameters{2,5};
                for i=1:obj.NOF
                    A=imread(fname,i);
                    A=preprocess_image(obj,set_n,i);
                    pos=pkfnd(A,threshold,feature_size);
                    if ~isempty(pos)
                        pos(find(pos(:,1)==0),:)=[]; %removes particles at origin that are tracking artifact
                    else
                       disp(['No particles were found in frame ' num2str(i) '.'])
                    end
                    obj.pos_all{i,set_n}=pos;
                end
            end
        end
        function obj=remove_overlap_pos(obj,threshold)
            %compares the positions between different sets in pos_all and
            %removes any particle positions that correspond to the same
            %particle, i.e. particles that are closer than within a
            %distance "threshold" (in pixels);
            if nargin<2
                threshold=1;%default is 1;
            end
            sizepos=size(obj.pos_all);
            nsets=sizepos(2);
            if nsets<2
                error('This function only works if there is more than one set of particle positions');
            end
            
            %pos_clean_temp=[];
            for i=1:obj.NOF %for each frame
                pos_i=[]; %make an empty matrix pos_i
                for n=1:nsets %then for each set of particle positions
                    pos_i=[pos_i; obj.pos_all{i,n}]; %fill pos_i with the positions of that set
                end%loop over sets
                
                %find the distance matrix for frame i
                nparts=length(pos_i(:,1));
                distance=NaN(nparts,nparts);
                for j=1:nparts%length(centers(:,1))
                    for k=1:nparts%length(centers(:,1))
                        %Get distance matrix
                        distance(j,k)=sqrt((pos_i(j,1)-pos_i(k,1))^2+(pos_i(j,2)-pos_i(k,2))^2);
                    end
                end
                
                %find duplicates
                numneighs=NaN(nparts,1);
                neighbours={};
                for p=1:nparts
                    neighbourstest=find(distance(p,:)<threshold & distance(p,:)>0);
                    neighbours_temp=neighbourstest(find(neighbourstest~=p));
                    neighbours{p}=neighbours_temp;
                    numneighs(p,1)=length(neighbourstest);
                end
                
                %remove duplicates
                pos_new=[];
                for p=1:nparts
                    if numneighs(p,1)==0 %if a particle has no duplicates
                        pos_new=[pos_new; pos_i(p,:)]; %add it to the new pos
                    elseif numneighs(p,1)>0 % otherwise, if a particle has at least one duplicate
                        %Do not add it to the new pos list, but instead
                        for q=1:length(neighbours{p}) %for each duplicate
                            neighbourID=neighbours{p}(q); %Find the ID of the duplicate
                            neighbours{neighbourID}=neighbours{neighbourID}(2:end); %remove the ID of particle p from the list of neighbours of that duplicate
                            numneighs(neighbourID)=numneighs(neighbourID)-1; %and reduce the number of neighbours of the duplicate by 1
                        end
                    end
                end
                pos_new(:,3)=i;
                obj.pos_clean{i,1}=pos_new; %pos_clean_temp=[pos_clean_temp ; pos_new];
            end%loop over frames
            %obj.pos_clean=pos_clean_temp;
        end
        function plot_positions_all(obj,i)
            if nargin<2
                i=1; %default is image 1
            end
            sizeposall=size(obj.pos_all);
            nsets=sizeposall(2);
            if nsets==0
                error('You have no positions in the table obj.pos_all to plot yet');
            end
            fname=[obj.foldername obj.filename '.tif'];
            A=imread(fname,i);
            size_A=size(A);
            figure
            hold on
            colormap('gray'),imagesc(A);
            axis([0 size_A(2) 0 size_A(1)]);
            cmap=jet(nsets);
            for n=1:nsets
                centers=obj.pos_all{i,n};
                radiplot=zeros(length(centers(:,1)),1)+5;
                viscircles(centers(:,1:2),radiplot(:,1),'Color',cmap(n,:));
            end
            hold off
        end
        function plot_positions_clean(obj,i)
            if nargin<2
                i=1; %default is image 1
            end
            if length(obj.pos_clean)==0
                error('You have no positions in the table obj.pos_clean to plot yet');
            end
            fname=[obj.foldername obj.filename '.tif'];
            A=imread(fname,i);
            size_A=size(A);
            centers=obj.pos_clean{i};
            radiplot=zeros(length(centers(:,1)),1)+5;
            figure
            hold on
            colormap('gray'),imagesc(A);
            axis([0 size_A(2) 0 size_A(1)]);
            if ~isempty(centers)
                viscircles(centers(:,1:2),radiplot(:,1));
            end
            hold off
        end
        function obj=split_types(obj,threshold,R)
            %splits all particles in each frame into two sets based on
            %their average gray value. Above a treshold gray value is set 1
            %and below is set 2
            %The sets are saved as matrices in the elements of the cell
            %output
            %The gray value of a particle is the average within a square of
            %edge R.
            %For the image: blue is above threshold, red is below threshold
            if isempty(obj.pos_clean)
                disp('You need to identify your particles before you can distinguish types')
                error('Find your positions and run the code remove_overlap_pos first')
            end
            if nargin<2
                threshold=50;
            end
            if nargin<3
                R=2;
            end
            fname=[obj.foldername obj.filename '.tif'];
            
            for t=1:obj.NOF
                pos1=[];
                pos2=[];
                %get the right image frame
                imtemp=imread(fname,t);
                meanim=mean(mean(imtemp));
                sizeim=size(imtemp);
                %add a border of size R around the image
                im=ones(sizeim(1)+2*R,sizeim(2)+2*R);
                im=im.*meanim;
                im(1+R:sizeim(1)+R,1+R:sizeim(2)+R)=imtemp;
                
                %select pos file of right frame here
                post=obj.pos_clean{t}; %pos(find(pos(:,3)==t),:);
                NOP=length(post(:,1));
                
                for p=1:NOP
                    xp=post(p,1);
                    yp=post(p,2);
                    area=im(round(yp):round(yp+2*R),round(xp):round(xp+2*R));
                    grayval=mean(mean(area));
                    if grayval>threshold
                        pos1=[pos1 ; post(p,:)+R];
                    else
                        pos2=[pos2 ; post(p,:)+R];
                    end
                end%end loop particles
                obj.pos_split{t,1}=pos1;
                obj.pos_split{t,2}=pos2;
                obj.splitparams={threshold R};
                %plot the positions of the two types of particles in
                %different colors. Only for first frame, to check
                if t==1
                    if length(pos1)>1
                        radii1=ones(length(pos1(:,1)),1)*10;
                    end
                    if length(pos2)>1
                        radii2=ones(length(pos2(:,1)),1)*10;
                    end
                    figure
                    hold on
                    colormap('gray'),imagesc(im);
                    if length(pos1)>1
                        viscircles(pos1(:,1:2), radii1, 'Color','b' )
                    end
                    if length(pos2)>1
                        viscircles(pos2(:,1:2), radii2, 'Color','r' )
                    end
                    hold off
                end
            end%end loop frames
        end
        function obj=track_obj_split(obj,set_n)
            %Takes pos_split and tracks the particle positions for each set
            %of positions separately
            %outputs a n x 4 matrix organized as [x y f ID]
            if nargin<2
                sizen=size(obj.tr_split);
                set_n=sizen(2)+1; %Default is the set that hasn't been chosen yet
                if set_n>length(obj.pos_split)
                    set_n=1;%If everything has already been analyzed, just analyze set 1
                end
                disp(['Analyzing particle position set ' num2str(set_n) '.']);
            end
            
            if isempty(obj.pos_split)
                error('obj.pos_split is still empty. Run split_types_pos first');
            end
            
            % Set parameters for trackingsoftware
            maxdisp=get_numerical_input('What is the max displacement? ');
            param.mem=get_numerical_input('How many frames can a particle go missing? ');
            param.dim=2; % dimensionality of data
            param.good=get_numerical_input('What is the minimal number of frames for a track to be accepted? ');
            param.quiet=1; % 0 = text, 1 = no text
            
            %prep the positions
            pos=[];
            for i=1:obj.NOF
                pos=[pos; obj.pos_split{i,set_n}];
            end
            
            moveon=1;
            try trace=track(pos,maxdisp,param);
            catch
                disp('Tracking failed. Try using different parameters');
                moveon=0;
            end
            if moveon==1
                obj.tr_split{1,set_n}=trace;
                obj.trparams{1,set_n}=param;
                obj.trparams{2,set_n}=maxdisp;
                %obj=complete_trace(obj);
            end
        end
        function plot_tracks(obj)
            if isempty(obj.tr_split)
                error('First track the particles using track_obj_split')
            end
            
            fname=[obj.foldername obj.filename '.tif'];
            A=imread(fname,1);
            size_A=size(A);
            cmap=jet(length(obj.tr_split));
            cmap=[0 0 1; 1 0 0; cmap];
            
            figure
            hold on
            box on
            colormap('gray'),imagesc(A);
            axis([0 size_A(2) 0 size_A(1)]);
            for n=1:length(obj.tr_split)
                traj=obj.tr_split{n};
                for j=1:max(traj(:,4))
                    trj=traj(find(traj(:,4)==j),1:3);
                    plot(trj(:,1),trj(:,2),'Color',cmap(n,:),'LineWidth',1);
                end
            end
        end
    end%end methods
end%end class