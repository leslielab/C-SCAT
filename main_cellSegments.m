function main_cellSegments(WindowsOrMac)

%% step one is to load the images as a 3D array 
% WindowsOrMac=0; % windows=1, mac=0, this to correct for '/'
currentFolder = pwd;

if WindowsOrMac==1 %
    load([currentFolder '\data\scattering_actin_part1.mat'])
    load([currentFolder '\data\scattering_actin_part2.mat'])
    addpath(genpath([currentFolder '\data']))
else
    load([currentFolder '/data/scattering_actin_part1.mat'])
    load([currentFolder '/data/scattering_actin_part2.mat'])
    addpath(genpath([currentFolder '/data']))
end

scattering_mat=zeros(size(scattering_video_array_firsthalf,1),size(scattering_video_array_firsthalf,2),size(scattering_video_array_firsthalf,3)*2);
scattering_mat(:,:,1:size(scattering_video_array_firsthalf,3))=scattering_video_array_firsthalf;
scattering_mat(:,:,(size(scattering_video_array_firsthalf,3)+1):end)=scattering_video_array_secondhalf;

if WindowsOrMac==1 %
    load([currentFolder '\data\fluo_actin_part1.mat'])
    load([currentFolder '\data\fluo_actin_part2.mat'])
else
    load([currentFolder '/data/fluo_actin_part1.mat'])
    load([currentFolder '/data/fluo_actin_part2.mat'])
end

fluoresence_mat=zeros(size(scattering_video_array_firsthalf,1),size(scattering_video_array_firsthalf,2),size(scattering_video_array_firsthalf,3)*2);
fluoresence_mat(:,:,1:size(scattering_video_array_firsthalf,3))=fluo_video_array_firsthalf;
fluoresence_mat(:,:,(size(scattering_video_array_firsthalf,3)+1):end)=fluo_video_array_secondhalf;


%% the cell segment analysis consists of manual selection, where that inforamtion is saved so that all regons can be processed afterwards

% first we need to quantify the dark current signal in the image by finding
% the lowest average value in the z-stack
bg_dark_use=100; % just a temporary value
for kk=1:size(scattering_mat,3)
    frame_scat=scattering_mat(:,:,kk);
    avg_bg=mean(frame_scat(:));
    if avg_bg<bg_dark_use
        bg_dark_use=avg_bg;
    end
end
%

frame_start=111; % a manually found starting point for which the the strucutre is approximately in focus
frames_search=4; % the focus range for which the signal will be estimated


px_size=49.3; % nn, pixel size in the image
One_micron_in_pixels=1000/px_size; % conversion parameter to obtain the signal corresponding to one 1 micrometer
Finding_segment_firstStep=0; % parameter to regulate the step through process


power_measurement=31.6; % illumination power for the specific measurement
power_beads=75; % illumination power for the reference beads measurement
signal_beads=19; % (integrated field, length unit corresponding to pixels squared) signal for the reference bead, used to relate signal to a polarizability estimate
pol_bead=3*(4*pi/3*55^3)*(1.58^2-1.33^2)/(1.58^2+2*1.33^2); % theoretical polarizability for the reference beads
form_factor_segment=0.8; % theoretically estimated correction, come from the width of the illumication beam

% parameters for the selected example segment
if Finding_segment_firstStep==1
    frame_vec=frame_start;
    selected_region_done=0;
else
    frame_vec=(-frames_search:frames_search)+frame_start;
    selected_region_done=1;
    x_sel=340:510;
    y_sel=470:640;

    rotation_done=1;
    final_selection_done=1;
    rotation=30;

    x_sel_small=50:200;
    y_sel_small=85:125;

    scat_signal=zeros(1,1);
    fluo_signal=zeros(1,1);
end

for jj=frame_vec % loop frames around focus

    % for the fluoresence we use a maxium projection (less background
    % signal than for scattering)
    frame_avg_test_fluo(:,:,1)=double(fluoresence_mat(:,:,jj));
    frame_avg_test_fluo(:,:,2)=double(fluoresence_mat(:,:,jj+1));
    frame_avg_test_fluo(:,:,3)=double(fluoresence_mat(:,:,jj+2));
    frame_avg_test_fluo(:,:,4)=double(fluoresence_mat(:,:,jj+3));
    frame_avg_test_fluo(:,:,5)=double(fluoresence_mat(:,:,jj-1));
    frame_avg_test_fluo(:,:,6)=double(fluoresence_mat(:,:,jj-2));
    frame_avg_test_fluo(:,:,7)=double(fluoresence_mat(:,:,jj-3));
    fluo_image=max(frame_avg_test_fluo,[],3);

    % start with subtracting dark current
    frame_scat=double(scattering_mat(:,:,jj))-bg_dark_use;
    frame_scat(frame_scat<0)=0;
    frame_scat_sqrt=sqrt(frame_scat); % then take the square root to perform the C-SCAT processing


    frame_bg=imgaussfilt(frame_scat_sqrt,25); % wide Gaussian to estimate background signal
    Cscat_im_sub=frame_scat_sqrt-frame_bg; % subtract background
    Cscat_im_sub=imgaussfilt(Cscat_im_sub,0.3); % convolve with a Gaussian to supress noise
    Cscat_im_sub=Cscat_im_sub *2/sqrt(1000)*sqrt(power_beads/power_measurement); % normalize the data and take into consideration differences between cell and reference bead data

    % plot the previously included segments (to avoid that the same segment
    % is included twice)
    if jj==frame_start
        figure(1001); clf
        imagesc(Cscat_im_sub)
        colormap('gray')
        hold on
        plot([504,600,618,535,504],[840,670,716,898,840],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([518,581,605,550,518],[598,479,493,616,598],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([477,578,616,525,477],[473,321,344,504,473],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([646,699,724,671,646],[451,354,373,476,451],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([777,838,857,803,777],[562,429,448,583,562],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([717,753,783,752,717],[940,835,843,967,940],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([565,600,635,607,565],[875,772,788,887,875],'color',[0.93,0.69,0.13],LineWidth=2)
        hold off

        figure(1002); clf
        imagesc(fluo_image)
        hold on
        plot([504,600,618,535,504],[840,670,716,898,840],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([518,581,605,550,518],[598,479,493,616,598],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([477,578,616,525,477],[473,321,344,504,473],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([646,699,724,671,646],[451,354,373,476,451],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([777,838,857,803,777],[562,429,448,583,562],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([717,753,783,752,717],[940,835,843,967,940],'color',[0.93,0.69,0.13],LineWidth=2)
        plot([565,600,635,607,565],[875,772,788,887,875],'color',[0.93,0.69,0.13],LineWidth=2)
        hold off
    end
    % step one is to create a rectangular subimage with the segment of
    % interest
    if selected_region_done==0
        disp('Select a new rectangular box around a straight segement part')
    end

    if selected_region_done==1
        im_small_field=Cscat_im_sub(x_sel,y_sel); % select the smaller region
        im_small_fluo=fluo_image(x_sel,y_sel);

        im_afterRot=imrotate(im_small_field,rotation,'bilinear'); % rotate the segment such that it is parallell to one axis (used during averaging=
        im_afterRot_fluo=imrotate(im_small_fluo,rotation,'bilinear');


        figure(1231); colormap('gray');
        imagesc(im_afterRot_fluo)


        if rotation_done==0
            disp('Find the rotation of the segment')
        end
        if rotation_done==1

            im_sel_afterRot=im_afterRot(x_sel_small,y_sel_small); % after rotation, center the segment
            im_sel_afterRot_fluo=im_afterRot_fluo(x_sel_small,y_sel_small);

            im_sel_afterRot_vec=mean(im_sel_afterRot,1); % average along the segment
            im_sel_afterRot_vec_fluo=mean(im_sel_afterRot_fluo,1)/100;

            figure(1232); colormap('gray');
            imagesc(im_sel_afterRot_fluo)
            
            if final_selection_done==0
                disp('Center the segment')
            end
            if final_selection_done==1
                
                % peform an 1D Gaussian fit to estiamte signal
                x=1:length(im_sel_afterRot_vec); x=x-mean(x);
                
                im_sel_afterRot_vec_fluo=im_sel_afterRot_vec_fluo-min(im_sel_afterRot_vec_fluo(:));
                im_sel_afterRot_vec_fluo=im_sel_afterRot_vec_fluo/100; % normalize so that the fluorescence and scattering have similar magnitudes
                %
                gaussEqn = 'a*exp(-(x-b)^2/c^2 )+d'; % segment and local background
                [f1,a1] = fit(x(:),im_sel_afterRot_vec(:),gaussEqn,'Start', [min(im_sel_afterRot_vec(:)),0,2,0], 'Lower', [min(im_sel_afterRot_vec(:))*2,-4,1,-0.1], 'Upper', [0,4,10,0.1]);
                [fluo1,a1] = fit(x(:),im_sel_afterRot_vec_fluo(:),gaussEqn,'Start', [mean(im_sel_afterRot_vec_fluo(:)),0,2,0], 'Lower', [min(im_sel_afterRot_vec_fluo(:)),-4,1,0], 'Upper', [max(im_sel_afterRot_vec_fluo(:))*1.2,4,10,max(im_sel_afterRot_vec_fluo(:))/2]);


                figure(1); plot(x,im_sel_afterRot_vec); hold on
                plot(x,im_sel_afterRot_vec_fluo-min(im_sel_afterRot_vec_fluo));

                figure(2111); 
                if jj==frame_vec(1)
                    clf
                end
                hold on
                plot(x,f1(x))
                plot(x,fluo1(x))

                hold off

                scatter_integral_tmp=One_micron_in_pixels*abs(f1.a)*(f1.c)*sqrt(pi);
                fluo_cell_structure_tmp=abs(fluo1.a)*(fluo1.c);

                % find the image with the highest particle signal for C-SCAT
                if scat_signal(1)<scatter_integral_tmp
                   scat_signal(1)= scatter_integral_tmp;

                   con_pol=confint(f1);
                end

                % find the image with the highest particle signal for
                % fluoresence
                if fluo_signal(1)<fluo_cell_structure_tmp
                    fluo_signal(1)=fluo_cell_structure_tmp;
                    con_fluo=confint(fluo1);
                end
            end
        end
    end
end

% convert signals and estimate uncertainty 
pol_cell_structure(1)=scat_signal(1)/signal_beads*pol_bead/form_factor_segment; % 
fluo_cell_structure(1)=fluo_signal(1);

pol_cell_uncert(1)=One_micron_in_pixels*sqrt(pi)/signal_beads*pol_bead/0.8*(-min(abs(con_pol(:,1)))*min(con_pol(:,3)) + max(abs(con_pol(:,1)))*max(con_pol(:,3)))/2;
fluo_cell_uncert(1)=(-min(abs(con_fluo(:,1)))*min(con_fluo(:,3)) + max(abs(con_fluo(:,1)))*max(con_fluo(:,3)))/2;

%% the segments included in the manuscript
% first we need to quantify the dark current signal in the image by finding
% the lowest average value in the z-stack
bg_dark_use=100; % just a temporary value
for kk=1:size(scattering_mat,3)
    frame_scat=scattering_mat(:,:,kk);
    avg_bg=mean(frame_scat(:));
    if avg_bg<bg_dark_use
        bg_dark_use=avg_bg;
    end
end
%



power_measurement=31.6; % illumination power for the specific measurement
power_beads=75; % illumination power for the reference beads measurement
signal_beads=19; % signal for the reference bead
pol_bead=3*(4*pi/3*55^3)*(1.58^2-1.33^2)/(1.58^2+2*1.33^2); % theoretical polarizability for the reference beads
form_factor_segment=0.8; % theoretically estimated correction, come from the width of the illumication beam

% a loop for all segements manually selected
partInd_vec=1:7;
pol_cell_structure=zeros(1,length(partInd_vec));
pol_cell_uncert=zeros(1,length(partInd_vec));
fluo_cell_structure=zeros(1,length(partInd_vec));
fluo_cell_uncert=zeros(1,length(partInd_vec));
scat_signal=zeros(1,length(partInd_vec));
fluo_signal=zeros(1,length(partInd_vec));

for part_ind=partInd_vec

    if part_ind==1
        frame_start=111;
        x_sel=340:510;
        y_sel=470:640;
        rotation=30;
        x_sel_small=50:200;
        y_sel_small=85:125;

        sel_box_coord_x=[477,578,616,525,477];
        sel_box_coord_y=[473,321,344,504,473];
    elseif part_ind==2
        frame_start=111;
        x_sel=470:660;
        y_sel=464:646;
        rotation=24;
        x_sel_small=40:160;
        y_sel_small=110:135;

        sel_box_coord_x=[518,581,605,550,518];
        sel_box_coord_y=[598,479,493,616,598];
    elseif part_ind==3
        frame_start=109;
        x_sel=740:915;
        y_sel=470:620;
        rotation=25;
        x_sel_small=40:180;
        y_sel_small=90:120;

        sel_box_coord_x=[504,600,618,535,504];
        sel_box_coord_y=[840,670,716,898,840];
    elseif part_ind==4
        frame_start=112;
        x_sel=360:500;
        y_sel=620:740;
        rotation=24;
        x_sel_small=30:150;
        y_sel_small=70:95;

        sel_box_coord_x=[646,699,724,671,646];
        sel_box_coord_y=[451,354,373,476,451];
    elseif part_ind==5
        frame_start=111;
        x_sel=400:600;
        y_sel=750:900;
        rotation=23;
        x_sel_small=79:200;
        y_sel_small=98:117;

        sel_box_coord_x=[777,838,857,803,777];
        sel_box_coord_y=[562,429,448,583,562];
    elseif part_ind==6
        frame_start=110;
        x_sel=860:980;
        y_sel=710:800;
        rotation=14;
        x_sel_small=20:110;
        y_sel_small=30:70;

        sel_box_coord_x=[717,753,783,752,717];
        sel_box_coord_y=[940,835,843,967,940];
    elseif part_ind==7
        frame_start=109;
        x_sel=720:880;
        y_sel=550:680;
        rotation=12;
        x_sel_small=50:150;
        y_sel_small=64:87;

        sel_box_coord_x=[565,600,635,607,565];
        sel_box_coord_y=[875,772,788,887,875];
    end

    frames_search=4; % range around the chosen focus to find the maximum signal
    frame_vec=(-frames_search:frames_search)+frame_start;
    for jj=frame_vec
        % for the fluoresence we use a maxium projection (less background
        % signal than for scattering)
        frame_avg_test_fluo(:,:,1)=double(fluoresence_mat(:,:,jj));
        frame_avg_test_fluo(:,:,2)=double(fluoresence_mat(:,:,jj+1));
        frame_avg_test_fluo(:,:,3)=double(fluoresence_mat(:,:,jj+2));
        frame_avg_test_fluo(:,:,4)=double(fluoresence_mat(:,:,jj+3));
        frame_avg_test_fluo(:,:,5)=double(fluoresence_mat(:,:,jj-1));
        frame_avg_test_fluo(:,:,6)=double(fluoresence_mat(:,:,jj-2));
        frame_avg_test_fluo(:,:,7)=double(fluoresence_mat(:,:,jj-3));
        fluo_image=max(frame_avg_test_fluo,[],3);

        % start with subtracting dark current
        frame_scat=double(scattering_mat(:,:,jj))-bg_dark_use;
        frame_scat(frame_scat<0)=0;
        frame_scat_sqrt=sqrt(frame_scat); % then take the square root to perform the C-SCAT processing


        frame_bg=imgaussfilt(frame_scat_sqrt,25); % wide Gaussian to estimate background signal
        Cscat_im_sub=frame_scat_sqrt-frame_bg; % subtract background
        Cscat_im_sub=imgaussfilt(Cscat_im_sub,0.3); % convolve with a Gaussian to supress noise
        Cscat_im_sub=Cscat_im_sub *2/sqrt(1000)*sqrt(power_beads/power_measurement); % normalize the data and take into consideration differences between cell and reference bead data

        % smaller image around segment
        im_small_field=Cscat_im_sub(x_sel,y_sel);
        im_small_fluo=fluo_image(x_sel,y_sel);
        
        % rotate segment to one axis
        im_afterRot=imrotate(im_small_field,rotation,'bilinear');
        im_afterRot_fluo=imrotate(im_small_fluo,rotation,'bilinear');

        % make so that the segment is in the center
        im_sel_afterRot=im_afterRot(x_sel_small,y_sel_small);
        im_sel_afterRot_fluo=im_afterRot_fluo(x_sel_small,y_sel_small);

        % average along the segment
        im_sel_afterRot_vec=mean(im_sel_afterRot,1);
        im_sel_afterRot_vec_fluo=mean(im_sel_afterRot_fluo,1)/100;

        x=1:length(im_sel_afterRot_vec); x=x-mean(x);
        im_sel_afterRot_vec_fluo=im_sel_afterRot_vec_fluo-min(im_sel_afterRot_vec_fluo(:));
        im_sel_afterRot_vec_fluo=im_sel_afterRot_vec_fluo/100;
        
        % 1D Gaussian fit to estiamte signal
        if max(im_sel_afterRot_vec(:))>-min(im_sel_afterRot_vec(:))
            im_sel_afterRot_vec=-im_sel_afterRot_vec;
        end
        gaussEqn = 'a*exp(-(x-b)^2/c^2 )+d';
        [f1,a1] = fit(x(:),im_sel_afterRot_vec(:),gaussEqn,'Start', [min(im_sel_afterRot_vec(:)),0,2,0], 'Lower', [min(im_sel_afterRot_vec(:))*2,-4,1,-0.1], 'Upper', [0,4,10,0.1]);
        [fluo1,a1] = fit(x(:),im_sel_afterRot_vec_fluo(:),gaussEqn,'Start', [mean(im_sel_afterRot_vec_fluo(:)),0,2,0], 'Lower', [min(im_sel_afterRot_vec_fluo(:)),-4,1,0], 'Upper', [max(im_sel_afterRot_vec_fluo(:))*1.2,4,10,max(im_sel_afterRot_vec_fluo(:))/2]);


        scatter_integral_tmp=One_micron_in_pixels*abs(f1.a)*(f1.c)*sqrt(pi);% 
        fluo_cell_structure_tmp=abs(fluo1.a)*(fluo1.c);
        
        % save maximum signal per segment
        if scat_signal(part_ind)<scatter_integral_tmp
            scat_signal(part_ind)= scatter_integral_tmp;

            con_pol=confint(f1);
        end

        if fluo_signal(part_ind)<fluo_cell_structure_tmp
            fluo_signal(part_ind)=fluo_cell_structure_tmp;
            con_fluo=confint(fluo1);
        end

    end

    % convert to polarizailbity, estimate fluoresence signal, and fit
    % uncertainty
    pol_cell_structure(part_ind)=scat_signal(part_ind)/signal_beads*pol_bead/form_factor_segment;
    fluo_cell_structure(part_ind)=fluo_signal(part_ind);

    pol_cell_uncert(part_ind)=One_micron_in_pixels*sqrt(pi)/signal_beads*pol_bead/form_factor_segment*(-min(abs(con_pol(:,1)))*min(con_pol(:,3)) + max(abs(con_pol(:,1)))*max(con_pol(:,3)))/2;
    fluo_cell_uncert(part_ind)=(-min(abs(con_fluo(:,1)))*min(con_fluo(:,3)) + max(abs(con_fluo(:,1)))*max(con_fluo(:,3)))/2;
end

%% plot the data (as in the manuscript)

RI_cyto=1.36; % refractive index of the cytoplams
dia_actin=0:0.1:73; % diameter range

figure(2342); clf; hold on
eff_RI=1.57; % refractive index of actin
Actin_pol_est=3*(dia_actin/2).^2*pi*1000.*(eff_RI.^2-RI_cyto^2)./((eff_RI.^2+2*RI_cyto^2)); % theoretical polarizailbity for one micron segment

plot(dia_actin,Actin_pol_est/(10^4),'linewidth',2) % plot theoretical curve

ylabel('Polarizability per \mum [10^4 nm]','Interpreter','tex')
xlabel('Diameter [nm]','Interpreter','tex')
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
box on
set(gca,'linewidth',2)
set(gca,'FontSize', 15)
set(gcf, 'PaperPositionMode', 'auto');
grid off
box on
set(gca,'LineWidth',2)
ylim([0,120])

% plot data
for kk=1:length(pol_cell_structure)
    if kk==1
        plot([0,max(dia_actin)],[1,1]*pol_cell_structure(kk)/10^4,'--','LineWidth',2,'Color',[0.93,0.69,0.13])
    else
        plot([0,max(dia_actin)],[1,1]*pol_cell_structure(kk)/10^4,'--','LineWidth',2,'Color',[0.93,0.69,0.13],'HandleVisibility','off')
        
    end

    [a,b]=min(abs(Actin_pol_est-pol_cell_structure(kk)));
    single_dia_vec(kk)=dia_actin(b);
end

legend('Refractive index=1.57','Measured segments','Location','northwest')

yyaxis right
% convert signal to number of actin filamenets
dia_actin=7;
ER_pol_single=3*(dia_actin/2).^2*pi*1000.*(eff_RI.^2-RI_cyto^2)./((eff_RI.^2+2*RI_cyto^2));
ylim([0,120/(ER_pol_single/10000)])
ylabel('Number actin filaments','Interpreter','tex')

number_filaments=pol_cell_structure/ER_pol_single;

xlim([0,72])
% saveas(gcf, '../results/actin_pol_size_plot.png');

%% plot polarizability-fluoresence 
figure(120312); 
clf
hold on
errorbar(pol_cell_structure/(10^4),fluo_cell_structure*100,fluo_cell_uncert*100,fluo_cell_uncert*100,pol_cell_uncert/(10^4),pol_cell_uncert/(10^4),'*','LineWidth',2)
ylim([0,60])
xlim([0,122])
ylabel('Integrated fluoresence signal')
xlabel('Polarizability per \mum [10^4 nm]')

LinEq = 'a*x';
[f1,a1] = fit(pol_cell_structure(:)/(10^4),fluo_cell_structure(:)*100,LinEq,'Start', 1, 'Lower', 0, 'Upper',10);
plot(0:180,(0:180)*f1.a,'--',LineWidth=2)
hold off

set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
box on
set(gca,'linewidth',2)
set(gca,'FontSize', 15)
set(gcf, 'PaperPositionMode', 'auto');
grid off
box on
set(gca,'LineWidth',2)
legend('Actin','Linear fit','Location','northwest')
pause(0.001)

% saveas(gcf, '../results/actin_pol_fluo_plot.png');