%LMTF DATA SANITY CHECK this function pulls down filenames from the
%database that do not have certain fields populated already, and
%repopulates those fields with the appropriate information. This is to
%ensure that all items in the LMTF table in the database have the same
%information gathered in the same fashion.
%
%AS OF 11/29/16, the data collected/input is as follows:
%   stim_dur - stimulus duration, rounded
%   displayLum - brightness of the display as measured by the
%   spectrophotometer, rounded to whole number
%Commented out information is not saved to the DB.
%EG
function LMTFDataSanityCheck()
%init
num_DP = 0;
conn = database('Nex_Paradigm_Sort','','','Vendor','MySql','Server','128.95.153.12'); %connect to sorted files database
flist = fetch(conn, 'SELECT fileID FROM LMTF WHERE subjID = ''G'' AND notes LIKE (''eye%'')');
    % Looking at color directions tested in the files
    for i = 1:length(flist())
        stro = nex2stro(findfile(flist{i},'N:\NexFiles'));
        % % stimulus durations
        %     stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
        %     stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));
        %     dur = stimoff_t-stimon_t;
        %     if sum(isnan(dur))>0
        %         dur(isnan(dur)) = [];
        %     end
        %     stim_dur = round(mean(dur)*1000); %in units of milliseconds
        %     % % get the framerate of the monitor at time of data collection
        %     framerate = stro.sum.exptParams.framerate;
        %     %nframes = dur*framerate;
        %      % % Getting monitor background luminance
        %     storedColorDirs = stro.trial(:,[11 12]);
        %     correctedColorDirs = mkbasis(storedColorDirs')';
        %     unique(correctedColorDirs);
        %     bkg = stro.sum.exptParams.bkgndrgb; %already been through the gamma table - if around .5 then that's halfway through monitor intensity
        %     ms = stro.sum.exptParams.mon_spd;
        %     if ~isnan(ms)
        %         try
        %             if size(ms,1)==303
        %                 mon_spd = reshape(ms, 101, 3);
        %                 wavelength_step = 4;
        %             else
        %                 mon_spd = reshape(ms, 201, 3);
        %                 wavelength_step = 2;
        %             end
        %         catch
        %             keyboard;
        %         end
        %         load('T_xyz1931')
        %         vlambda = T_xyz1931(2,:);
        %         [mon_spd] = SplineSpd([380:wavelength_step:780]', mon_spd, [380:5:780]');
        %         lum = 680*vlambda*mon_spd;
        %         displayLum = lum*bkg; %in units of cdm^-2
        %     else
        %         displayLum = 0;
        %     end
        %   CHECK NUMBER OF DATA POINTS PER FILE
        file_dp = unique(stro.trial(:,15)');
        num_DP = num_DP + file_dp(end); %save just num_dp(end)
        %   FILE TFS TESTED
        %   tfs_tested = unique(stro.trial(:,13)'); %maybe save all?
        %     updateDB_query = sprintf('UPDATE LMTF SET dispLum = %d, stimDur = %d WHERE fileID = ''%s''', round(displayLum,0), stim_dur, flist{i});
        %     updateDB = exec(conn, updateDB_query);
        %     if updateDB.Message
        %         keyboard;
        %     else
        %         continue;
        %     end
    end
    disp(num_DP);
close(conn);
end