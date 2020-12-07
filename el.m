function varargout=el(option,varargin)
% EVLAB tools
%
% PREPROCESSING SYNTAX:
%
%   el('preprocessing',subjectID [, functional_runs , preprocessing_pipeline]) imports from DICOM and preprocesses the functional&structural data for subject 'subjectID'
%      subjectID                : subject folder (e.g. '408_FED_20160617a_3T2')
%      functional_runs          : (optional) array of DICOM session numbers associated with functional data (e.g. [1,2,3]). If unspecified, this information will be automatically derived from the DICOM file headers
%      preprocessing_pipeline   : (optional) .cfg file defining preprocessing steps (e.g. pipeline_preprocessing_DefaultMNI.cfg). If unspecified, the pipeline_preproc_DefaultMNI_Denoise_PlusStructural.cfg pipeline will be used
%                                 (in the absence of structural data the pipeline_preproc_DefaultMNI_Denoise.cfg pipeline will be used instead)
%
%      e.g.  el('preprocessing', '408_FED_20160617a_3T2');
%
%
%   el('preprocessing.append',subjectID, preprocessing_pipeline) runs additional preprocessing steps to the (already previously preprocessed) functional&anatomical data for subject 'subjectID'
%
%      e.g.  el('preprocessing.append', '408_FED_20160617a_3T2', '/data/preproc_addsteps.cfg');
%
%
%   submitID = el('submit','preprocessing',subjectID [, functional_runs , preprocessing_pipeline]) runs import&preprocessing steps on remote node
%   submitID = el('submit','preprocessing.append',subjectID, preprocessing_pipeline) runs additional preprocessing steps on remote node
%
%      e.g.  el('submit', 'preprocessing', '408_FED_20160617a_3T2');
%      e.g.  el('submit', 'preprocessing.append', '408_FED_20160617a_3T2', '/data/preproc_addsteps.cfg');
%
%   dataID = el('preprocessing',...) returns an evlab17*.mat file identifying one imported&preprocessed dataset (for cases when several different import steps are run, e.g. each dataset focuses on different sessions which may require different preprocessing steps)
%   dataID = el('preprocessing.append',dataID,...) runs additional preprocessing steps on selected dataset
%
%
% CONFIGURATION OPTIONS:
%
%   el('default',fieldname) returns value of configuration setting fieldname
%   el('default',fieldname,fieldvalue) changes value of configuration setting fieldname to fieldvalue
%   
%   el('default','folder_subjects',foldername)
%      foldername               : default root directory where subject folders -e.g. 408_FED_20160617a_3T2- may be found. By default this is defined as evlab17/ROOTFOLDER
%   el('default','folder_dicoms',foldername)
%      foldername               : subdirectory name within folder_subjects where DICOM files may be found. Default '*dicoms'
%   el('default','dicom_disregard_functional',list)
%      list                     : list of keywords defining SeriesDescription values of DICOM sessions to be excluded from the list of valid functional runs
%   el('default','dicom_isstructural',list)
%      list                     : list of keywords defining SeriesDescription values of DICOM sessions to be interpreted as anatomical/structural scans
% 
% FIRST-LEVEL ANALYSES:
% el firstlevel <subjectID> <modelID>
%

persistent defaults;
if isempty(defaults), 
    defaults=struct(...
        'folder_subjects',fullfile(fileparts(which(mfilename)),'ROOTFOLDER') ,...
        'folder_dicoms','*dicoms' ,...
        'dicom_isstructural',{{'^T1_MPRAGE_1iso'}} ,...
        'dicom_disregard_functional',{{'^localizer','^AAScout','^AAHScout','^MoCoSeries','^T1_MPRAGE_1iso','^DIFFUSION_HighRes'}} );
end


evlab17 init silent;
fileout=[];
varargout=cell(1,nargout);

switch(lower(option))
    case 'default'
        if isempty(varargin), varargout={defaults}; 
        elseif isfield(defaults,varargin{1})
            if numel(varargin)>1, defaults.(varargin{1})=varargin{2};
            else varargout={defaults.(varargin{1})};
            end
        else
            if numel(varargin)>1, evlab17('default',varargin{:});
            else varargout={evlab17('default',varargin{1})};
            end
        end
    case 'submit'
        if ~nargout, conn('submit',@el,varargin{:}); % e.g. el submit preprocessing 408_FED_20160617a_3T2
        else [varargout{1:nargout}]=conn('submit',@el,varargin{:});
        end
            
    case 'preprocessing'
        % adapted from msieg preprocess_PL2017
        % el('preprocessing',subject_id [, functional_runs , preprocessing_pipeline_file])
        % e.g. el preprocessing 408_FED_20160617a_3T2
        
        pwd0=pwd;
        subject=char(varargin{1}); % subject id
        subject_path=fullfile(defaults.folder_subjects,subject);
        subject_path_dicoms = fullfile(subject_path,defaults.folder_dicoms);
        if any(subject_path_dicoms=='*'), subject_path_dicoms = conn_dir(subject_path_dicoms ,'-dir','-R','-cell','-sort'); end
        assert(~isempty(subject_path_dicoms),'unable to find dicom folder %s\n',fullfile(subject_path,defaults.folder_dicoms)); 
        if iscell(subject_path_dicoms), subject_path_dicoms=subject_path_dicoms{1}; end
        %subject_keys=regexp(subject,'_','split');
        %subject_path_dicoms = fullfile(subject_path,strjoin([subject_keys(2:4),{'dicoms'}],'_')); % 268_FED_20170929a_3T2_PL2017_unsmoothed/FED_20170929a_3T2_dicoms, 268_KAN_EvDB_20150317a_PL2017/KAN_EvDB_20150317a_dicoms, 230_KAN_prodsemphon_12_PL2017/KAN_prodsemphon_12_dicoms, 183_POLY01_20160420_3T1/POLY01_20160420_3T1_dicoms

        func_runs=[];
        struct_run=[];
        if numel(varargin)<2||isempty(varargin{2})
            Series=conn_dcmdir(fullfile(subject_path_dicoms,'*-1.dcm'),false);
            idx=find(cellfun('length',{Series.SeriesDescription})>0&cellfun('length',{Series.SeriesNumber})>0);
            SeriesDescription={Series(idx).SeriesDescription};
            SeriesNumber=[Series(idx).SeriesNumber];
            struct_run = SeriesNumber(find(cellfun(@(x)~isempty(regexp(char(x),strjoin(defaults.dicom_isstructural,'|'))),SeriesDescription)>0,1,'first')); % keep first structural
            func_runs = setdiff(SeriesNumber(find(cellfun(@(x)isempty(regexp(char(x),strjoin(defaults.dicom_disregard_functional,'|'))),SeriesDescription)>0)),struct_run); % keep all functionals
            assert(~isempty(func_runs),'unable to find any functional runs in %s\n',subject_path_dicoms);
            try, 
                fid=fopen(fullfile(subject_path_dicoms,'runs.csv'),'wt');
                fprintf(fid,'SeriesNumber,SeriesDescription\n');
                for n=1:numel(SeriesNumber), fprintf(fid,'%s,%s\n',num2str(SeriesNumber(n)),SeriesDescription{n}); end
                fclose(fid);
                fprintf('DICOM series information stored in %s\n',fullfile(subject_path_dicoms,'runs.csv'))
            end
        else % run numbers
            func_runs=varargin{2};
            if ischar(func_runs), func_runs=str2num(func_runs); end
        end
            
        if numel(varargin)<3||isempty(varargin{3})
            if isempty(struct_run), preproc_config_file = 'pipeline_preproc_DefaultMNI_Denoise.cfg';
            else preproc_config_file = 'pipeline_preproc_DefaultMNI_Denoise_PlusStructural.cfg';
            end
        else preproc_config_file = varargin{3};
        end
            
        all_runs=[struct_run func_runs];
        fid=fopen(fullfile(subject_path,'data.cfg'),'wt');
        fprintf(fid,'\n#dicoms\n');
        for n=1:numel(all_runs), fprintf(fid,'%s\n',fullfile(subject_path_dicoms,['*-' num2str(all_runs(n)) '-1.dcm'])); end
        fprintf(fid,'\n#functionals\n');
        fprintf(fid,'%s\n',num2str(func_runs));
        if ~isempty(struct_run)
            fprintf(fid,'\n#structurals\n');
            fprintf(fid,'%s\n',num2str(struct_run));
        end
        fprintf(fid,'\n#RT nan\n'); % note: will read RT from .json files
        fclose(fid);
        
        %run preproc
        [ok,msg]=mkdir(fullfile(subject_path,'nii'));
        cd(fullfile(subject_path,'nii'));
        fileout=evlab17_run_preproc(fullfile(subject_path,'data.cfg'),preproc_config_file,[],varargin(4:end));
        cd(pwd0);
        varargout{1}=fileout; 
    
    case 'preprocessing.append'
        % el('preprocessing.append',subject_id, preprocessing_pipeline_file)
        
        pwd0=pwd;
        subject=char(varargin{1}); % subject id
        if isempty(regexp(subject,'\.mat$'))
            subject_path=fullfile(defaults.folder_subjects,subject,'nii');
            files=conn_dir(fullfile(subject_path,'evlab17*.mat'),'-ls');
            assert(~isempty(files), 'unable to find any evlab17*.mat files in %s\n',subject_path);
            subject=files{end};
        else subject=conn_fullfile(subject); 
        end
        evlab17_module('load',subject);
        preproc_config_file = varargin{2}; % preprocessing pipeline
        if isempty(preproc_config_file), preproc_config_file={[]};
        else preproc_config_file={preproc_config_file, []};
        end
        %run preproc
        cd(fileparts(subject));
        fileout=evlab17_run_preproc(preproc_config_file{:},'dataset',subject, varargin(3:end));
        cd(pwd0);
        varargout{1}=fileout; 
    
    case 'firstlevel'
        return
        % adapted from msieg firstlevel_PL2017
        %E.G. el('firstlevel','408_FED_20160617a_3T2','langlocSN')
        subject = char(varargin(1));
        expt = char(varargin(2));
        
        %Get contrasts
        fid = fopen(fullfile(defaults.rootfolder,'ANALYSIS/contrasts_by_expt.txt'));
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
        cd(fullfile(defaults.rootfolder,defaults.subjectfolder,subject));
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
        fprintf(fid, fullfile(defaults.rootfolder,defaults.subjectfolder,subject,'nii',pp_file));
        
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ space '\n']);
        fprintf(fid, [ header2 '\n']);
        fprintf(fid, fullfile(defaults.rootfolder,defaults.subjectfolder,subject,cat_file));
        
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
        evlab17_run_model(fullfile(defaults.rootfolder,defaults.subjectfolder,subject,['modelfiles_' expt '.cfg']),'pipeline_model_Default.cfg');
        
        cd(fullfile(defaults.rootfolder,'ANALYSIS'));



    otherwise
        [varargout{1:nargout}]=evlab17(option,varargin{:});
end

end
