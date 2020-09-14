clear; close all; clc

% % % hard-coding tempvec
temp_lims = [100 300]; 

% directories where to get stuff.
% datadir = '/Users/molly/Downloads/laura pilot data/';
% datadir = '/Users/molly/Downloads/TEMPITAPS_april/ece pilot/';
% stimdir = '/Users/molly/Downloads/TEMPITAPS_april/stimorders/';

datadir = '/Users/nicoleahuizinga/Desktop/TEMPITAPS_april/Ece/';
stimdir = '/Users/nicoleahuizinga/Desktop/TEMPITAPS_april/stimorders/';
%stimdir = '/Users/molly/Downloads/';

subnr = 2; % populate this

trial{1} = [1 2 3];
trial{2} = [1 2 3];%1:40; % can use this structure to throw out some trials in case they are trash

for s = subnr
    [trial_thresh,num_outliers,time_outliers] = deal([]);
    
    stimorder = [zeros(length(trial{s}),1), 2*ones(length(trial{s}),1)]; % TEMPORARY
    load(strcat(stimdir,'BKKmirrorss.mat'))
%     load(strcat(stimdir,'Rstimorder_',num2str(s),'.mat'))
%     load(strcat(stimdir,'Rstimorder_3.mat'))

      for kk = trial{s}
        file = strcat(datadir,'main_pt',num2str(s),'_trial',num2str(kk),'.wav');
%         file = strcat(datadir,'main_pt42_trial',num2str(kk),'.wav');
        [data,sf] = audioread(file);
        t  = 1/sf:1/sf:length(data)/sf;         
        
        % get the stimulus ------------------------------------------------
        y = BKKmirrors{kk};
        y = y(1:end-sf); % this chops off exactly 1 s, which is the noise
        tstim = 1/sf:1/sf:length(y)/sf;
%         figure,  plot(tstim,y)
        
        % % get rid of babies (noise floor) -----------------------------------------
        bthresh = .01; % this threshold is manual for now
        y(y < bthresh) = 0; % set everything below threshold to 0
%             figure, plot(tstim,y) % plot the nice clean signal
        
        % first derivative --------------------------------------------------------
        deriv = diff(y); % get the first derivative
        Td = tstim(1:end-1); % get a new time vector
%         figure, plot(Td,deriv) % plot the derivative            
        
        % find onsets based on derivative % --------------------------------------
        Thresh = bthresh; % why not just keep the threshold chosen above?
        all_onsets = find(deriv > Thresh); % find all onsets (this will be wayyyyyy too many)
        tones = []; % = taps
        distance = diff(all_onsets); % what is the distance between each "onset"?
        for ii = 1:length(all_onsets)-1 % loop through all of them
            if ii == 1 % always keep the first one
                tones = [tones all_onsets(ii)];
            else
                if distance(ii-1) > .05*sf % this is another manual threshold, here I don't take onsets that are closer than 100 ms to each other (people can't tap that fast)
                    tones = [tones all_onsets(ii)];
                end
            end
        end
        
      end
      
end
% %             figure, plot(Td,deriv), hold on % plot the derivative
% %             plot(Td(tones),0,'or') % with the onsets on top
        
        % make 1s & 0s out of the red circles (Td(tones))
        Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
        Tdlgth(tones) = 1;
%         figure, stem(Td,Tdlgth)

        rhythm_onsets = Td(tones);
%         aud_delay = rhythm_onsets(1);
%         rhythm_onsets_zeroed = rhythm_onsets - aud_delay;
        rhythm_onsets_zeroes = rhythm_onsets - rhythm_onsets(1);

        nints = length(diff(rhythm_onsets));
%         if stimorder(kk,2) == 1
%             tempvec = temp_lims(1) * (temp_lims(2)/temp_lims(1)).^((0:nints-1)/(nints-1)); % can make only for individual rhythm, need to know stimulus identity.
%         elseif stimorder(kk,2) == 2
%             tempvec = temp_lims(2) * (temp_lims(1)/temp_lims(2)).^((0:nints-1)/(nints-1)); % can make only for individual rhythm, need to know stimulus identity.
%         end  

        
        % get the tapping data (this plot if for inspection / thresholding)
        x = data(:,1);
%         figure,  plot(t,x)
        
        % this is a new little interactive module for you to set the threshold
        % on a trial-by-trial basis. Not fully automated, but neither is my
        % brain right now.
        ok = 0;
        while ok == 0
            figure,  plot(t,x)

            % % get rid of babies (noise floor) -----------------------------------------
            bthresh = input('Enter threshold to try:  '); % this threshold is manual for now
            x(x < bthresh) = 0; % set everything below threshold to 0
            %     figure, plot(t,x) % plot the nice clean signal
            
            % first derivative --------------------------------------------------------
            deriv = diff(x); % get the first derivative
            Td = t(1:end-1); % get a new time vector
            figure, plot(Td,deriv) % plot the derivative
            
            check = input('Is it ok (y=1|n=0)?  ');
            close all
            ok = check;
        end
        trial_thresh(kk) = bthresh; % save the single-trial threshold 
        
        % find onsets based on derivative % --------------------------------------
        Thresh = bthresh; % why not just keep the threshold chosen above?
        all_onsets = find(deriv > Thresh); % find all onsets (this will be wayyyyyy too many)
        tones = []; % = taps
        distance = diff(all_onsets); % what is the distance between each "onset"?
        for ii = 1:length(all_onsets)-1 % loop through all of them
            if ii == 1 % always keep the first one
                tones = [tones all_onsets(ii)];
            else
                if distance(ii-1) > .1*sf % this is another manual threshold, here I don't take onsets that are closer than 100 ms to each other (people can't tap that fast)
                    tones = [tones all_onsets(ii)];
                end
            end
        end
        %     figure, plot(Td,deriv), hold on % plot the derivative
        %     plot(Td(tones),0,'or') % with the onsets on top
        
        % make 1s & 0s out of the red circles (Td(tones))
        Tdlgth = zeros(size(Td)); %a matrix of zeros the same size as Td - a vector here
        Tdlgth(tones) = 1;
        %     figure, stem(Td,Tdlgth) % plot the TAPS ONLY
        
        %% now how to analyze this?
        taptimes_raw = Td(tones);
%         taptimes_raw = taptimes_raw - aud_delay; % we no longer know this
%         taptimes_zeroed = taptimes_raw - taptimes_raw(1);
        taptimes_zeroed = taptimes_raw - rhythm_onsets(1); % now have to assume rhythm play is perfect until new info
        
        % need the specific rhythm info here (will need stimorder specific to
        % subject to decode)
        % get stimulus times raw and temp vec per subject
        
        % get the intervals between taps
        ints = diff(taptimes_zeroed);
%         figure, subplot(1,2,1), plot(ints,'ok-'), hold on
        figure, plot(taptimes_zeroed(1:end-1),ints,'ok-'), hold on

%         scatter(1:length(ints),ints,'k')%,lsline
        
        % plot the tempvec with the produced tap series.
        if stimorder(kk,2) == 1
            oversampled_tempvec = temp_lims(1) * (temp_lims(2)/temp_lims(1)).^((0:10000-1)/(10000-1)); % can make only for individual rhythm, need to know stimulus identity.
        elseif stimorder(kk,2) == 2
            oversampled_tempvec = temp_lims(2) * (temp_lims(1)/temp_lims(2)).^((0:10000-1)/(10000-1)); % can make only for individual rhythm, need to know stimulus identity.
        end
        oversampled_tempvec = oversampled_tempvec / 1000;
        oversampled_timevec = linspace(0,rhythm_onsets(end),10000);
        quart_temp = oversampled_tempvec;
        half_temp  = 2*oversampled_tempvec;
        whole_temp = 4*oversampled_tempvec;

        plot(oversampled_timevec,quart_temp,'m')
        plot(oversampled_timevec,half_temp,'c')
        plot(oversampled_timevec,whole_temp,'g')

        
        tap_temp = [];
        %for tt = 1:length(taptimes_raw)
            %tap_temp = [tap_temp oversampled_tempvec(nearest(oversampled_timevec,taptimes_raw(tt)))];
        %   end
        
        % in this case tt is 113,
        for tt = 1:length(taptimes_raw)
            tap_temp = [tap_temp oversampled_tempvec(closeby(oversampled_timevec,taptimes_raw(tt)))];
        end
        
        % get the temp vecs for each metrical level
        quart_tap = tap_temp(1:end-1);
        half_tap  = 2*tap_temp(1:end-1);
        whole_tap = 4*tap_temp(1:end-1);
        
        % plot all the metrical levels on the produced interval time series
%         plot(1:length(ints),quart_temp,'m')
%         plot(1:length(ints),half_temp,'c')
%         plot(1:length(ints),whole_temp,'g')
        %         figure, scatter(tap_temp(1:end-1),ints)
        
        % determine which metrical level the intervals belong to, and clean
        % the ones that are then considered outliers
        out_perc = .5; % percentage of intended interval to be within, if outside of this, outlier
        % this is oversampled for plotting
        quart_range = [quart_temp + (quart_temp * out_perc); quart_temp - (quart_temp * out_perc)];
        half_range  = [half_temp + (half_temp * out_perc); half_temp - (half_temp * out_perc)];
        whole_range = [whole_temp + (whole_temp * out_perc); whole_temp - (whole_temp * out_perc)];
        
        % and this is for individual taps for calculations. 
        quart_tap_range = [quart_tap + (quart_tap * out_perc); quart_tap - (quart_tap * out_perc)];
        half_tap_range  = [half_tap + (half_tap * out_perc); half_tap - (half_tap * out_perc)];
        whole_tap_range = [whole_tap + (whole_tap * out_perc); whole_tap - (whole_tap * out_perc)];

        plot(oversampled_timevec,quart_range(1,:),'m--',oversampled_timevec,quart_range(2,:),'m--')
        plot(oversampled_timevec,half_range(1,:),'c--',oversampled_timevec,half_range(2,:),'c--')
        plot(oversampled_timevec,whole_range(1,:),'g--',oversampled_timevec,whole_range(2,:),'g--')
        
        % now remove outliers ---------------------------------------------
        % this is to remove some at the beginning in case they aren't
        % stable yet
        ok = 0;
        while ok ~= 1
            figure, plot(taptimes_zeroed(1:end-1),ints,'ok-'), hold on
            plot(taptimes_zeroed(1:end-1),quart_tap,'m')
            plot(taptimes_zeroed(1:end-1),half_tap,'c')
            plot(taptimes_zeroed(1:end-1),whole_tap,'g')
            plot(taptimes_zeroed(1:end-1),quart_tap_range(1,:),'m--',taptimes_zeroed(1:end-1),quart_tap_range(2,:),'m--')
            plot(taptimes_zeroed(1:end-1),half_tap_range(1,:),'c--',taptimes_zeroed(1:end-1),half_tap_range(2,:),'c--')
            plot(taptimes_zeroed(1:end-1),whole_tap_range(1,:),'g--',taptimes_zeroed(1:end-1),whole_tap_range(2,:),'g--')
            how_many_unstable = input('How many unstable taps to remove?  ');
            plot(1:how_many_unstable,ints(1:how_many_unstable),'rx','MarkerSize',10)
            check = input('Is this ok?  ');
            close all
            ok = check;
        end
        ints(1:how_many_unstable) = nan;
        
        % first need to correct the metrical level assignment in case there
        % is an accidental one (then it will be compared to the wrong
        % outlier criteria                
        % first choose which level
        met_lvl = [];
        for mla = 1:length(ints)
            d = abs([ints(mla)-quart_tap(mla) ints(mla)-half_tap(mla) ints(mla)-whole_tap(mla)]);
            met_lvl = [met_lvl find(d == min(d))];
        end

        met_thresh = 2; % this many in a run that fall outside the neighbors will be discareded
        nd = diff(met_lvl); % level change direction from tap to tap
        md_idx = [];
        for ll = 1:length(nd)-met_thresh
            if nd(ll) ~= 0
                tmp = nd(ll:ll+met_thresh);
                
            
        
        
        
        
        [outliers, clean_trialdata] = deal([]);
        for ay = 1:length(ints)
            if met_lvl(ay) == 1
                if ints(ay) <= quart_tap_range(1,ay) && ints(ay) >= quart_tap_range(2,ay)
                    clean_trialdata(ay) = ints(ay);
                else
                    clean_trialdata(ay) = nan;
                    outliers = [outliers, [ay; ints(ay)]];
                end
            elseif met_lvl(ay) == 2
                if ints(ay) <= half_tap_range(1,ay) && ints(ay) >= half_tap_range(2,ay)
                    clean_trialdata(ay) = ints(ay);
                else
                    clean_trialdata(ay) = nan;
                    outliers = [outliers, [ay; ints(ay)]];
                end
            elseif met_lvl(ay) == 3
                if ints(ay) <= whole_tap_range(1,ay) && ints(ay) >= whole_tap_range(2,ay)
                    clean_trialdata(ay) = ints(ay);
                else
                    clean_trialdata(ay) = nan;
                    outliers = [outliers, [ay; ints(ay)]];
                end
            end
        end
        figure
        subplot(1,2,1), plot(taptimes_zeroed(1:end-1),ints,'ok-'), hold on
        plot(taptimes_zeroed(outliers(1,:)),outliers(2,:),'xr','MarkerSize',10)
        plot(taptimes_zeroed(1:end-1),quart_tap,'m')
        plot(taptimes_zeroed(1:end-1),half_tap,'c')
        plot(taptimes_zeroed(1:end-1),whole_tap,'g')
        plot(taptimes_zeroed(1:end-1),quart_tap_range(1,:),'m--',taptimes_zeroed(1:end-1),quart_tap_range(2,:),'m--')
        plot(taptimes_zeroed(1:end-1),half_tap_range(1,:),'c--',taptimes_zeroed(1:end-1),half_tap_range(2,:),'c--')
        plot(taptimes_zeroed(1:end-1),whole_tap_range(1,:),'g--',taptimes_zeroed(1:end-1),whole_tap_range(2,:),'g--')

        subplot(1,2,2), plot(taptimes_zeroed(1:end-1),clean_trialdata,'ok-'), hold on
        plot(taptimes_zeroed(1:end-1),quart_tap,'m')
        plot(taptimes_zeroed(1:end-1),half_tap,'c')
        plot(taptimes_zeroed(1:end-1),whole_tap,'g')
        plot(taptimes_zeroed(1:end-1),quart_tap_range(1,:),'m--',taptimes_zeroed(1:end-1),quart_tap_range(2,:),'m--')
        plot(taptimes_zeroed(1:end-1),half_tap_range(1,:),'c--',taptimes_zeroed(1:end-1),half_tap_range(2,:),'c--')
        plot(taptimes_zeroed(1:end-1),whole_tap_range(1,:),'g--',taptimes_zeroed(1:end-1),whole_tap_range(2,:),'g--')

        pause
        close all

        
        
        
        
        
        % how to tell if there's a switch?
        sw = [];
        switch_thresh = .2; % yet another threshold
        for jj = 1:length(ints)-1
            %         tmp(jj) = ints(jj+1) - ints(jj);
            if abs(ints(jj+1) - ints(jj)) > switch_thresh
                sw = [sw jj+1];
            end
        end
        plot(sw,(ints(sw)),'+r','MarkerSize',25)
        
        
        
        
        
        
    end
    
end
    
    %% time point at which participant started tapping & this variability
    
    %% identifying when tapping behavior changes.....
    
    %% mean and std of the intervals
    
    m = mean(ints);
    sd = std(ints);
    
    %% landfill - extra code from pre-online version
    
    
    
    
