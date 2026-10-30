function automateTrack_particle_Fluo_DIC_iSCATandDF_clean(file, outputFileName, outDir,WindowsOrMac)
%%%%%
% inputs are the path to the .nd2 file (containing the data) + name and path
% to where the output data will be saved in a folder
% 
% The output is a .mat file with the iSCAT, C-SCAT, DIC, fluoresnece, and
% position information for each particle above the threshold
%%%%%

% Check that the file exists
% if ~isfile(file)
%     error('File does not exist: %s', file);
% end
if ~isfolder(file)
    error('Folder does not exist: %s', file);
end

% Attempt to open the file
try
    %     result = bfopen(file);
    fileList = dir(fullfile(file, '*.mat'));
    for k =  1:length(fileList) % [16,14,12,10,9,8,7,6,5,4,3]%
        %     load(fileList(k).name)
        if WindowsOrMac==1 %    windows
            load([fileList(k).folder '\' fileList(k).name])
        else
            load([fileList(k).folder '/' fileList(k).name]) % on MacOS, "\" needs to be changed to "/"
        end
    end
catch ME
    rethrow(ME);
end

if k==9
    fluo_video_array=zeros(size(fluo_video_array_part1,1),size(fluo_video_array_part1,2),size(fluo_video_array_part1,3)+size(fluo_video_array_part2,3));
    Cscat_red_video_array=zeros(size(fluo_video_array_part1,1),size(fluo_video_array_part1,2),size(fluo_video_array_part1,3)+size(fluo_video_array_part2,3));
    Cscat_blue_video_array=zeros(size(fluo_video_array_part1,1),size(fluo_video_array_part1,2),size(fluo_video_array_part1,3)+size(fluo_video_array_part2,3));
    DIC_video_array=zeros(size(fluo_video_array_part1,1),size(fluo_video_array_part1,2),size(fluo_video_array_part1,3)+size(fluo_video_array_part2,3));

    fluo_video_array(:,:,1:size(fluo_video_array_part1,3))=fluo_video_array_part1;
    fluo_video_array(:,:,(size(fluo_video_array_part1,3)+1):end)=fluo_video_array_part2;

    Cscat_red_video_array(:,:,1:size(fluo_video_array_part1,3))=Cscat_red_array_part1;
    Cscat_red_video_array(:,:,(size(fluo_video_array_part1,3)+1):end)=Cscat_red_array_part2;

    Cscat_blue_video_array(:,:,1:size(fluo_video_array_part1,3))=Cscat_blue_array_part1;
    Cscat_blue_video_array(:,:,(size(fluo_video_array_part1,3)+1):end)=Cscat_blue_array_part2;

    DIC_video_array(:,:,1:size(fluo_video_array_part1,3))=DIC_video_array_part1;
    DIC_video_array(:,:,(size(fluo_video_array_part1,3)+1):end)=DIC_video_array_part2;
end

% side length for a ROI
roiS = 32; 

% number of simultanous measured channels
im_s = 4; 

% % order of the data in the .nd2 files 
% iSCAT1=0;
% fluo_s=1;
% iSCAT2=2;
% DIC_s=3;

%%%%% initializing arrays used to save the data later
fluo_stack=zeros(size(Cscat_red_video_array,1),size(Cscat_red_video_array,2),size(Cscat_red_video_array,3));
iSCAT1_stack=fluo_stack;
iSCAT2_stack=fluo_stack;
iSCAT2_bg_stack=fluo_stack;
DIC_stack=fluo_stack;
DIC_stack_bg=fluo_stack;
Field1_sqrtCorr_stack=fluo_stack;
Field2_sqrtCorr_stack=fluo_stack;
dark_current_iSCAT1=400000;
dark_current_iSCAT2=400000; % from data min mean value around 80-100. 
bg_vec = zeros(1,size(Cscat_red_video_array,3));
bg_vec2 = zeros(1,size(Cscat_red_video_array,3));




if 1==1% tf==1 % for consistency, only galvano data is included as scanning mode affects properties such as signal-to-noise
    firstloop=1; % related to the how the data is saved when the .nd2 file contains more than recording

    for mm=1 % loop over all videos in the file
        it = 0;

        for j = 1 : size(fluo_stack,3) % loop over all frames each recording
            % load the data for each frame
            im_fluo = double(fluo_video_array(:,:,j));
            im_iSCAT2 = double(Cscat_red_video_array(:,:,j));
            iSCAT_1 = double(Cscat_blue_video_array(:,:,j));
            im1_DIC = double(DIC_video_array(:,:,j));

            % dark current estimation for C-SCAT, and background estimation for DIC
            if it == 0 %
                for jj = 1  : size(fluo_stack,3)
                    DIC_stack_bg(:,:,jj) =double(DIC_video_array(:,:,jj));

                    im1_tmp = double(Cscat_blue_video_array(:,:,jj));
                    im2_tmp = double(Cscat_red_video_array(:,:,jj));

                    dark_current_iSCAT1=min([dark_current_iSCAT1,mean(im1_tmp(:))]);
                    dark_current_iSCAT2=min([dark_current_iSCAT2,mean(im2_tmp(:))]);
                end
                DIC_bg=median(DIC_stack_bg,3); % since DIC is a transmission methods, the average background signal can be estiamted using the full stack
            end

            % subtract dark current and obtain C-SCAT
            tmp_field1=iSCAT_1-dark_current_iSCAT1;
            tmp_field1(tmp_field1<0)=0;
            im_field_1=sqrt(tmp_field1);
            tmp_field2=im_iSCAT2-dark_current_iSCAT2;
            tmp_field2(tmp_field2<0)=0;
            im_field_2=sqrt(tmp_field2);
            iSCAT_1=iSCAT_1-dark_current_iSCAT1;
            % have iSCAT processing as reference 
            im_iSCAT2=im_iSCAT2-dark_current_iSCAT2;

            it = it + 1;

            % in some videos it might be needed to remove scanning line artifacts
            % aa=mean(im1,2);  aa3=aa./mean(im1(:)); im1=im1./aa3;
            % aa=mean(im_iSCAT2,2);  aa3=aa./mean(im_iSCAT2(:)); im_iSCAT2=im_iSCAT2./aa3;
            % aa=mean(im_field_1,2);  aa3=aa./mean(im_field_1(:)); aa3(aa3<1)=1; im_field_1=im_field_1./aa3;
            % aa=mean(im_field_2,2);  aa3=aa./mean(im_field_2(:)); aa3(aa3<1)=1; im_field_2=im_field_2./aa3;


            % the background signal for all image modalities except for DIC are estimated by convolving the image with a wide Gaussian
            mean_bg = imgaussfilt(iSCAT_1, 20);% mean(Bg_frames_iSCAT, 3);
            bg_iSCAT2=imgaussfilt(im_iSCAT2, 20);
            bg_iSCAT2_toField=im_iSCAT2;
            bg_fluo= imgaussfilt(im_fluo, 25);
            bg_field_1=imgaussfilt(im_field_1, 20);
            bg_field_2=imgaussfilt(im_field_2, 20); 

            % save the average background per frame
            bg_vec(mm,j) = mean(iSCAT_1(:));
            bg_vec2(mm,j) = mean(im_iSCAT2(:));

            % for each image modality, subtract background and normalize
            % (in the case of iSCAT and DIC processing)
            iSCAT_1 = real(iSCAT_1)- mean_bg;
            iSCAT_1 = (iSCAT_1 ./ (mean_bg)); 
            iSCAT_1=real(iSCAT_1);
            iSCAT_1(isnan(iSCAT_1)) = 0;
            iSCAT_1(isinf(iSCAT_1)) = 0;


            im_DIC = real(im1_DIC) - DIC_bg*mean(im1_DIC(:))/mean(DIC_bg(:));
            im_DIC = (im_DIC ./ DIC_bg);
            im_DIC(isnan(im_DIC)) = 0;
            im_DIC(isinf(im_DIC)) = 0;
            aa2_DIC = real(im_DIC);
            aa3_DIC = imgaussfilt(aa2_DIC, 1);


            im_fluo = im_fluo - bg_fluo;
            im_fluo = imgaussfilt(im_fluo, 1); % convolve to improve SNR

            im_iSCAT2  = (im_iSCAT2 - bg_iSCAT2)./(bg_iSCAT2);
            im_iSCAT2=real(im_iSCAT2);
            im_iSCAT2(isnan(im_iSCAT2)) = 0;
            im_iSCAT2(isinf(im_iSCAT2)) = 0;

            im_field_1=im_field_1-bg_field_1;
            im_field_2=im_field_2-bg_field_2; 

  
            % save each processed image/modality in a stack
            iSCAT1_stack(:,:,it)=iSCAT_1;
            iSCAT2_stack(:,:,it)=im_iSCAT2;
            DIC_stack(:,:,it)=aa3_DIC;
            fluo_stack(:,:,it)=im_fluo;
            iSCAT2_bg_stack(:,:,it)=bg_iSCAT2_toField*norm_ratio_field^2;
            Field1_sqrtCorr_stack(:,:,it)=im_field_1;
            Field2_sqrtCorr_stack(:,:,it)=im_field_2*norm_ratio_field;
        end
    
        % the fluoresence signal is used to create a mask indicating where
        % particles are located in the subseqent processing
        % median of a maxk([],3) projection is used to minimize sensetivity
        % to particles that are suspended in the solution 
        aaa_fluo_proj_tmp=maxk(fluo_stack,5,3);
        fluo_proj=median(aaa_fluo_proj_tmp,3);

        th=200; % fluoresence signal threshold
        % the fluoresnece projection and all image stacks are then sent to
        % a function which estiamtes the signal on the single-particle level
        [Positions_table_tmp, ~] = QuantSingleImageProj_clean(fluo_proj,roiS,th,iSCAT1_stack,iSCAT2_stack,fluo_stack,DIC_stack,Field1_sqrtCorr_stack,Field2_sqrtCorr_stack,px_size,mm,iSCAT2_bg_stack);

        % save all the single-particle data in a large array. The if-case
        % is used to merge data from multiple videos in case that the .nd2
        % file contains multiple recordings in the same file
        Positions_tmp = table2array(Positions_table_tmp);
        if firstloop==1
            Positions=Positions_tmp;
        else
            Positions=[Positions_tmp;Positions];
        end
        firstloop=0;

    end
    % save all single-particle information in a .mat file
    save(fullfile(outDir, [outputFileName '.mat']), 'Positions','Positions_table_tmp','bg_vec','bg_vec2','px_size','px_size_Z')
end

end