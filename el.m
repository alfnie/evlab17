function varargout=el(option,varargin)
%
% PREPROCESSING:
% el preprocess <subjectID> [, functional_runs , preprocessing_pipeline_file])
% el preprocess <subjectID> <modelID>
%
% FIRST-LEVEL ANALYSES:
% el firstlevel <subjectID> <modelID>
%


persistent rootfolder;

evlab17 init silent;
fileout=[];
varargout=cell(1,nargout);
ROOTFOLDER=el('utils_rootfolder');; 

switch(option)
    case {'rootfolder','utils_rootfolder'}
            if nargin>1&&~isempty(varargin{1}), rootfolder=varargin{1}; end
            if isempty(rootfolder), varargout={fullfile(fileparts(which(mfilename)),'ROOTFOLDER')};
            else varargout={rootfolder};
            end
            
    case 'preprocessing'
        % adapted from msieg preprocess_PL2017
        % el('preprocessing',subject_id [, functional_runs , preprocessing_pipeline_file])
        % e.g. el('preprocessing','408_FED_20160617a_3T2',[ 3 7 9 11 13 15 17 21 23])
        
        subject=char(varargin(1));
        fd = find(subject=='_');
        
        if strmatch('FED',subject(fd(1)+1:fd(2)-1))>0;
            sub = subject(5:fd(3)+3);
        elseif strmatch('KAN',subject(fd(1)+1:fd(2)-1))>0;
            sub = subject(5:fd(3)+2);
        elseif strmatch('POLY',subject(fd(1)+1:fd(2)-1))>0;
            sub = subject(5:fd(3)+5);
        else
            'subject name problem';
        end
        
        
        %go to subject
        pathSubject=fullfile(rootfolder,'SUBJECTS',subject);
        
        %Create run_summary in python
        cd(pathSubject);
        commandStr = ['python ' fullfile(rootfolder,'Useful_Scripts/summarize_experiments_preproc.py') ' -s ' subject ' -o run_summary.txt'];
        [status, commandOut] = system(commandStr)
        
        %check arguments
        if length(varargin) == 1;
            %Get functional runs from run summary if not entered manually
            rsID = fopen(fullfile(pathSubject,'run_summary.txt'),'rt');
            text = textscan(rsID,'%s');
            text=text{:};
            
            func_inds1 = strmatch('RUNS',text);
            func_inds2 = strmatch('EXPERIMENT',text);
            func_runs = cellfun(@str2num,text(func_inds1+1:func_inds2-1))';
            
            % check for mprage
            matches_mprage = strfind(text,'T1_MPRAGE_1iso');
            any_mprage = any(vertcat(matches_mprage{:}));
            
            if any_mprage;
                mprage_ind = strmatch('T1_MPRAGE_1iso',text);
                mprage_ind=mprage_ind(1);
                struct_run1=char(text(mprage_ind-2));
                struct_run=str2num(struct_run1(1:end-1))
                
                %define preprocessing config file if not entered manually
                preproc_config_file = 'pipeline_preproc_DefaultMNI_PlusStructural.cfg';
            elseif ~any_mprage;
                preproc_config_file = 'pipeline_preproc_DefaultMNI.cfg';
            end
            
            
            
        elseif length(varargin) == 2;
            
            %manual func runs
            func_runs=varargin(2);
            func_runs=func_runs{:};
            
            %define preprocessing config file if not entered manually
            preproc_config_file = 'pipeline_preproc_DefaultMNI_PlusStructural.cfg';
            
        elseif length(varargin) == 3;
            
            %manual func runs
            func_runs=varargin(2);
            func_runs=func_runs{:};
            
            %manual preproc config
            preproc_config_file = char(varargin(3));
        end
        
        
        
        %Open and write data.cfg
        space = ' ';
        header1 = '#dicoms';
        header2 = '#functionals';
        header3 = '#structurals';
        header4 = '#RT';
        
        fid=fopen('data.cfg','w');
        
        if any_mprage;
            all_runs=[struct_run func_runs];
        elseif ~any_mprage;
            all_runs = func_runs;
        end
        
        fprintf(fid, [ header1 '\n']);
        for i = 1:length(all_runs);
            fprintf(fid, fullfile(rootfolder,'SUBJECTS',subject,[sub '_dicoms'],['*-' num2str(all_runs(i)) '-1.dcm']));
            fprintf(fid, [ space '\n']);
        end;
        
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ header2 '\n']);
        fprintf(fid, [num2str(func_runs) '\n']);
        
        if any_mprage;
            fprintf(fid, [ space '\n']);
            fprintf(fid, [ space '\n']);
            fprintf(fid, [ header3 '\n']);
            fprintf(fid, [num2str(struct_run) '\n']);
        end
        
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ header4 '\n']);
        fprintf(fid, ['2' '\n']);
        
        fclose(fid);
        
        %run preproc
        cd ../../ANALYSIS
        %addpath('/software/evlab17');
        %addpath('/software/conn/');
        
        mkdir(fullfile(rootfolder,'SUBJECTS',subject,'nii'));
        cd(fullfile(rootfolder,'SUBJECTS',subject,'nii'));
        
        %keyboard
        evlab17_run_preproc(fullfile(rootfolder,'SUBJECTS',subject,'data.cfg'),preproc_config_file);

    
    case 'firstlevel'
        % adapted from msieg firstlevel_PL2017
        %E.G. el('firstlevel','408_FED_20160617a_3T2','langlocSN')
        subject = char(varargin(1));
        expt = char(varargin(2));
        
        %Get contrasts
        fid = fopen(fullfile(rootfolder,'ANALYSIS/contrasts_by_expt.txt'));
        cons = textscan(fid, '%s', 'Delimiter', '\n');
        cons = cons{:};
        f_expt = strmatch(expt,cons);
        expt_cons = cons(f_expt(1)+1:f_expt(2)-1);
        
        %If npmod
        if length(varargin) > 2 & strmatch(varargin(3),'npmod');
            num_events = cell2mat(varargin(4));
            len_con = length(expt_cons);
            for i = 1:len_con;
                for j = 1:num_events;
                    %modcon
                    spl_con = strsplit(char(expt_cons(i)));
                    con = char(spl_con(1));
                    spl_con(1) = {[char(spl_con(1)) '_' sprintf('%02d',j)]};
                    f_mod = find(mod(1:length(spl_con),2)==0);
                    for k = 1:length(f_mod);
                        spl_con(f_mod(k)) = {[char(spl_con(f_mod(k))) '_EVENT' sprintf('%02d',j)]};
                    end;
                    join_con = strjoin(spl_con);
                    mod_cons(j+num_events*(i-1),1) = {join_con};
                    
                    %expt_con1
                    expt_cons1(i,1) = {con};
                    expt_cons1(i,j+1) = {strjoin(spl_con(2:end))};
                end;
                expt_cons2(i,1) = {strjoin(expt_cons1(i,:))};
            end
            expt_cons=[expt_cons2; mod_cons];
        end
        
        %Get cat file name
        cd(fullfile(rootfolder,'SUBJECTS',subject));
        cat_file = ls(['*' expt '.cat']);
        %cat_file = conn_dir(fullfile(rootfolder,'SUBJECTS',subject,['*' expt '.cat']),'-ls');
        %assert(numel(cat_file)==1,'either multiple or no matches to %s',fullfile(rootfolder,'SUBJECTS',subject,['*' expt '.cat']));

        %Get preprocessing .mat output name
        cd nii
        p = pwd;
        pp_files = dir([p '/evlab*mat']);
        pp_files = {pp_files.name};
        pp_file = char(pp_files(end));
        cd ..
        
        %Open and write modelfiles*.cfg
        header1 = '#dataset';
        space = ' ';
        header2 = '#design';
        header3 = '#model_name';
        header4 = '#contrasts';
        fid=fopen(strcat('modelfiles_',expt,'.cfg'),'w');
        
        fprintf(fid, [ header1 '\n']);
        fprintf(fid, fullfile(rootfolder,'SUBJECTS',subject,'nii',pp_file));
        
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ header2 '\n']);
        fprintf(fid, fullfile(rootfolder,'SUBJECTS',subject,cat_file));
        
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ header3 '\n']);
        fprintf(fid, [expt '\n']);
        
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ header4]);
        for i = 1:length(expt_cons);
            fprintf(fid, [ space '\n']);
            fprintf(fid, char(expt_cons(i,:)));
        end;
        
        fclose(fid);
        
        %run firstlevel
        cd /software/evlab17/
        evlab17_run_model(fullfile(rootfolder,'SUBJECTS',subject,['modelfiles_' expt '.cfg']),'pipeline_model_Default.cfg');
        
        cd(fullfile(rootfolder,'ANALYSIS'));



    otherwise
        [varargout{1:nargout}]=evlab17(option,varargin{:});
end

end
