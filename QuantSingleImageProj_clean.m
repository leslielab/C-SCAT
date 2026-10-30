function [Positions, ROImat] = QuantSingleImageProj_clean(FluoImage,roiS,th,iSCAT1_stack,iSCAT2_stack,fluo_stack,DIC_stack,Field1_stack,Field2_stack,pixelsize,videoNr,bg_iSCAT2)
%%%%%
% FluoImage is used to identify potential particle canidates (fluorescence projection image), where the local maxima needs to larger than th
% The data processing is based on selecting a region of interest (ROI) for which all the processing and signal estimation is done in all modalites
% 
% For each particle, the postition (x,y,z), fluoresence signal, DIC signal,
% iSCAT estiamte, C-SCAT estimate, signal-to-background, background, and
% video number are saved
% Some of the used ROIs are also saved in case of post inspection is wanted
%%%%%

%%%%% potential particle identification based on local maxima approach
%%%%% (image dilation of a non-binary image)

SE_tmp=strel('disk',9); % local maxima withing a disc of radius 9
ROImat=cell(1,1); % to save the ROIs

% the local maxima over threshold code
dS=FluoImage;
dS3=imdilate(dS,SE_tmp);
Indices=(find(and((dS3==dS),dS>th)));
[~, Inds]=sort(dS(Indices),'descend');
Inds = (Inds);

% creating a table for which all the data will be saved
Positions = table();
x = zeros(length(Inds),1);
y = x;
iSCAT1_save = x;
iSCAT2_save = x;
Field1_save = x;
Field2_save = x;
DIC_save = x;
Fluo_save = x;
Rsq_DIC_save = x;
Z_slice_save = x;
local_iSCAT_bg_save = x;
sig_bg_ratio=x;
Rsq_FIELD_save=x;

px0 = 0.0493;

% x and y grids that are used during the fits
x_fit=1:roiS;
x_fit=x_fit-mean(x_fit);
y_fit=x_fit;
[X_f,Y_f]=meshgrid(x_fit,y_fit);
inc_index=0;

gaussEqn = 'a*exp(- (sqrt((x-b)^2+(y-c)^2 ))^2/d^2  )';% function used to fit the data 

for k=1:length(Inds)% loop through all particles 

    %%%%% for each particle candidate, extract a ROI window for each modality of interest 
    I=Indices(Inds(k));
    [rq, cq] = ind2sub(size(dS), I);

    if or(cq+roiS/2>size(FluoImage,2)-(roiS/2+1), rq+roiS/2>size(FluoImage,1)-(roiS/2+1)) % in case the particle is to close to the border of the image, exclude that particle
        continue; % In MATLAB, continue skips any remaining statements in the body of the loop for the current iteration. The code continues execution from the next iteration (here particle candidate)
    elseif or(cq-roiS/2<(roiS/2),rq-roiS/2<(roiS/2)) 
        continue;
    end
    idx = (cq-roiS/2+1):(cq+roiS/2);
    idy = (rq-roiS/2+1):(rq+roiS/2);

    % for each particle, the frame for which the central pixels for the
    % fluoresence signal is maximum is used as the z-position
    FluoSmall_stack=fluo_stack(idy, idx,:);
    [~,ind1] = max(FluoSmall_stack,[],3);
    test_ind1=round(mean(ind1(roiS/2-2:roiS/2+3,roiS/2-2:roiS/2+3),"all"));

    % if the z-position occurs to close to the end interval of the z-stack, exclude that particle
    if or(test_ind1<3,size(FluoSmall_stack,3)-test_ind1<3)
        continue;
    end
    % save the fluoresence ROI with the highest signal for later processing
    FluoSmall=fluo_stack(idy, idx,test_ind1);

    % for the iSCAT and DIC processing, the maximum projection of an
    % absolute image within +/- 2 frames from the z-position is used to
    % estimate the highest particle signal 
    iSCAT1_frame=test_ind1;
    iSCAT1_small=iSCAT1_stack(idy, idx,iSCAT1_frame-2:iSCAT1_frame+2);
    iSCAT1_small=max(abs(iSCAT1_small),[],3);

    DIC_frame=test_ind1;
    DIC_small=DIC_stack(idy, idx,DIC_frame-2:DIC_frame+2);
    DIC_small=max(abs(DIC_small),[],3);

    iSCAT2_frame=test_ind1;
    iSCAT2_small=iSCAT2_stack(idy, idx,iSCAT2_frame-2:iSCAT2_frame+2);
    bg_iSCAT2_small=bg_iSCAT2(idy, idx,iSCAT2_frame-2:iSCAT2_frame+2);
    iSCAT2_small=max(abs(iSCAT2_small),[],3);


    % a reduced ROI is here used to estimate the local background around
    % the particle candidate 
    red_fit_s=9;
    sh_x=0;
    sh_y=0;
    X_f_red=X_f(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);
    X_f_red=X_f_red-mean(X_f_red(:));
    Y_f_red=Y_f(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);
    Y_f_red=Y_f_red-mean(Y_f_red(:));

    border_sel=zeros(size(X_f_red));
    border_sel(1:2,:)=1;
    border_sel(end-1:end,:)=1;
    border_sel(:,1:2)=1;
    border_sel(:,end-1:end)=1;
    border_sel=border_sel>0.5;

    % for the C-SCAT processing, the maximum projection of an
    % absolute image (after local background subraction) within +/- 2 frames from the z-position is used to
    % estimate the highest particle signal 

    field1_frame=test_ind1;
    Field1_small=Field1_stack(idy, idx,field1_frame-2:field1_frame+2);
    for jj=1:size(Field1_small,3)
        field2_im_tmp=squeeze(Field1_small(:,:,jj));
        median_bg=median(field2_im_tmp(border_sel));
        Field1_small(:,:,jj)=Field1_small(:,:,jj)-median_bg;
    end
    Field1_small=max(abs(Field1_small),[],3);

    field2_frame=test_ind1;
    Field2_small=Field2_stack(idy, idx,field2_frame-2:field2_frame+2);
    for jj=1:size(Field2_small,3) 
        field2_im_tmp=squeeze(Field2_small(:,:,jj));
        median_bg=median(field2_im_tmp(border_sel));
        Field2_small(:,:,jj)=Field2_small(:,:,jj)-median_bg;
    end
    Field2_small=max(abs(Field2_small),[],3);





    
    % reduce the ROI for the Gaussian fits so that goodness of fit fior the Gaussian fit is more
    % sensetive to the particle and less so of the background 
    iSCAT1_red=iSCAT1_small(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);
    iSCAT2_red=iSCAT2_small(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);
    DIC_red=DIC_small(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);
    Fluo_red=FluoSmall(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);
    Field1_red=Field1_small(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);
    Field2_red=Field2_small(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y);

    % convolve image to supress noise
    pix_cor=0.058/pixelsize; % let the bounds scale with the pixel size
    sigma=2*pix_cor;
    iSCAT1_red = imgaussfilt(iSCAT1_red,sigma);
    iSCAT2_red = imgaussfilt(iSCAT2_red,sigma);
    DIC_red = imgaussfilt(DIC_red,sigma);
    Fluo_red = imgaussfilt(Fluo_red,sigma);
    Field1_red=imgaussfilt(Field1_red,sigma);
    Field2_red_no_conv=Field2_red;
    Field2_red=imgaussfilt(Field2_red,sigma);

    % subtract background signal using the border of the ROI
    iSCAT1_red=iSCAT1_red-median(iSCAT1_red(border_sel));
    iSCAT2_red=iSCAT2_red-median(iSCAT2_red(border_sel));
    DIC_red=DIC_red-median(DIC_red(border_sel));
    Fluo_red=Fluo_red-median(Fluo_red(border_sel));
    Field1_red=Field1_red-median(Field1_red(border_sel));
    Field2_red_no_conv=Field2_red_no_conv-median(Field2_red_no_conv(border_sel));
    Field2_red=Field2_red-median(Field2_red(border_sel));


    % fit the data using a Guassian
    [f1,a1] = fit([X_f_red(:),Y_f_red(:)],iSCAT1_red(:),gaussEqn,'Start', [max(iSCAT1_red(:))*0.8,0,0,2*pix_cor], 'Lower', [max(iSCAT1_red(:))*0.1,-4,-4,1.5]*pix_cor, 'Upper', [max(iSCAT1_red(:))*1.0,4*pix_cor,4*pix_cor,5*pix_cor]);
    [f2,a2] = fit([X_f_red(:),Y_f_red(:)],iSCAT2_red(:),gaussEqn,'Start', [max(iSCAT2_red(:))*0.8,0,0,2*pix_cor], 'Lower', [max(iSCAT2_red(:))*0.1,-2.5,-2.5,1.5]*pix_cor, 'Upper', [max(iSCAT2_red(:))*1.0,2.5*pix_cor,2.5*pix_cor,5*pix_cor]);
    [f3,a3] = fit([X_f_red(:),Y_f_red(:)],DIC_red(:),gaussEqn,'Start', [max(DIC_red(:))*0.8,0,0,2.1*pix_cor], 'Lower', [max(DIC_red(:))*0.1,-4,-4,1.8]*pix_cor, 'Upper', [max(DIC_red(:))*1.0,4*pix_cor,4*pix_cor,4*pix_cor]);
    [f4,a4] = fit([X_f_red(:),Y_f_red(:)],Fluo_red(:),gaussEqn,'Start', [max(Fluo_red(:))*0.8,0,0,3*pix_cor], 'Lower', [0,-3,-3,2.5]*pix_cor, 'Upper', [max(Fluo_red(:))*1,3*pix_cor,3*pix_cor,4.75*pix_cor]);
    [f6,a6] = fit([X_f_red(:),Y_f_red(:)],Field1_red(:),gaussEqn,'Start', [max(Field1_red(:))*0.8,0,0,3*pix_cor], 'Lower', [max(Field1_red(:))*0.1,-3,-3,1.5]*pix_cor, 'Upper', [max(Field1_red(:))*1.3,3*pix_cor,3*pix_cor,5*pix_cor]);
    [f7_No_conv,~] = fit([X_f_red(:),Y_f_red(:)],Field2_red_no_conv(:),gaussEqn,'Start', [max(Field2_red(:))*0.8,0,0,2.5*pix_cor], 'Lower', [max(Field2_red(:))*0.0,-3,-3,2.0]*pix_cor, 'Upper', [max(Field2_red(:))*2,3*pix_cor,1.5*pix_cor,3.3*pix_cor]);
    [f7,a7] = fit([X_f_red(:),Y_f_red(:)],Field2_red(:),gaussEqn,'Start', [max(Field2_red(:))*0.8,0,0,4.7*pix_cor], 'Lower', [max(Field2_red(:))*0.0,-3,-3,4.2]*pix_cor, 'Upper', [max(Field2_red(:))*1.3,3*pix_cor,3*pix_cor,6*pix_cor]);

    % subpixel position for the particles
    x(k) = cq+f4.b;
    y(k) = rq+f4.c;

    % save the Guassian fit integrals as well as the R-square values
    % if the -square valueis to low for the scattering, then it is
    % noise/uncertain fit and will not be included
    Field2_save(k) = f7.a*(a7.rsquare>0.3)*(f7.d)^2*pi*(pixelsize/px0)^2;
    iSCAT1_save(k)=abs(f1.a)*(a1.rsquare>0.3)*(f1.d)^2*pi*(pixelsize/px0)^2;
    iSCAT2_save(k)=abs(f2.a)*(a2.rsquare>0.3)*(f2.d)^2*pi*(pixelsize/px0)^2;
    Field1_save(k) = f6.a*(a6.rsquare>0.3)*(f6.d)^2*pi*(pixelsize/px0)^2;
    Rsq_DIC_save(k)=a3.rsquare;
    Rsq_FIELD_save(k)=a7.rsquare;
    DIC_save(k)=abs(f3.a)*(a3.rsquare>-0.1)*(f3.d)^2*pi*(pixelsize/px0)^2;
    Fluo_save(k)=abs(f4.a)*(a4.rsquare>0.5)*(f4.d)^2*pi*(pixelsize/px0)^2;
    Z_slice_save(k)=test_ind1;


    % estimate the local background to estimate the signal to background ratio
    bg_iSCAT2_small_tmp=bg_iSCAT2_small(roiS/2-red_fit_s+sh_x:roiS/2+red_fit_s+sh_x,roiS/2-red_fit_s+sh_y:roiS/2+red_fit_s+sh_y,:);
    sel_tmp=zeros(size(bg_iSCAT2_small_tmp,1),size(bg_iSCAT2_small_tmp,2),size(bg_iSCAT2_small_tmp,3));
    sel_tmp(1:2,:,:)=1;
    sel_tmp(end-1:end,:,:)=1;
    sel_tmp(:,1:2,:)=1;
    sel_tmp(:,end-1:end,:)=1;
    sel_tmp=sel_tmp>0.5;

    bg_tmp=mean(bg_iSCAT2_small_tmp(sel_tmp));
    if bg_tmp<0.1
        bg_tmp=0.1;
    end
    local_iSCAT_bg_save(k)=bg_tmp;%
    % the signal in the signal to background ratio is the C-SCAT signal with any prior convolution as a convolution affects the overall signal magnitude 
    sig_bg_ratio(k)=f7_No_conv.a/sqrt(local_iSCAT_bg_save(k));

    % save the ROIs
    inc_index=inc_index+1;
    ROImat{inc_index}=[Field2_small,DIC_small,Field2_small];

end

% save the data, but ignore the particles that are skipped (whenever the
% code executes 'continue')
y = y(x ~= 0);
iSCAT1_save = iSCAT1_save(x ~= 0);
iSCAT2_save = iSCAT2_save(x ~= 0);
Field1_save = Field1_save(x ~= 0);
Field2_save = Field2_save(x ~= 0);

DIC_save = DIC_save(x ~= 0);
Fluo_save = Fluo_save(x ~= 0);
Rsq_DIC_save=Rsq_DIC_save(x ~= 0);
Rsq_FIELD_save=Rsq_FIELD_save(x ~= 0);
Z_slice_save=Z_slice_save(x ~= 0);
local_iSCAT_bg_save=local_iSCAT_bg_save(x ~= 0);
sig_bg_ratio=sig_bg_ratio(x ~= 0);
x = x(x ~= 0);

Positions.X = x(:);
Positions.Y = y(:);
Positions.iSCAT1_save = iSCAT1_save(:);
Positions.iSCAT2_save = iSCAT2_save(:);
Positions.DIC_save = DIC_save(:);
Positions.Fluo_save = Fluo_save(:);
Positions.RsqDIC = Rsq_DIC_save(:);
Positions.Z_slice = Z_slice_save(:);
Positions.local_iSCAT2_bg=local_iSCAT_bg_save(:);
Positions.videoNr = Z_slice_save(:)*0+videoNr;
Positions.s_r_ratio=sig_bg_ratio(:);
Positions.Field1 = Field1_save(:);
Positions.Field2 = Field2_save(:);
Positions.Field2_gof = Rsq_FIELD_save(:);



end