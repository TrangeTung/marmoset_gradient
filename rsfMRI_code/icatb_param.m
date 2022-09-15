function sesInfo = icatb_param(input_folder,outputDir,filename,ICs,prefix)
load ICA_ica_parameter_info.mat;
sesInfo.outputDir = outputDir;
sesInfo.numComp = ICs;
sesInfo.ICA_Options = {'epsilon',0.000100000000000000,'maxNumIterations',1000,'maxFinetune',5,'sampleSize',1,'numOfIC',40,'verbose','on','approach','defl','finetune','off','stabilization','off','g','pow3','only','all'};
sesInfo.userInput.ICA_Options = sesInfo.ICA_Options;
sesInfo.ICA_Options{10} = sesInfo.numComp;
for i = 1:length(input_folder)
    cd(input_folder{i});
    dir_f = dir(filename);
    name = dir_f.name;
    sesInfo.inputFiles(i).name = '';
    all_func = MY_select_file_for_SPM(input_folder,['^',name],[1 Inf]);
    for j = 1:length(all_func)
        sesInfo.inputFiles(i).name(j,:) = all_func{j};
    end
    sesInfo.icaOutputFiles(2).ses(i).name = [];
    sesInfo.icaOutputFiles(2).ses(i).name = strcat(prefix,'_sub01_component_ica_s',num2str(i),'_.nii');
end
sesInfo.numOfScans = length(dir_f);
sesInfo.diffTimePoints = repmat(sesInfo.numOfScans,[1,length(input_folder)]);
sesInfo.numOfSub = length(input_folder);
sesInfo.numOfSess = 1;
sesInfo.numOfDataSets = length(input_folder);

data_info = spm_vol(strcat(dir_f(1).folder,'\',dir_f(1).name));
sesInfo.HInfo.DIM = data_info.dim;
sesInfo.HInfo.VOX = abs(diag(data_info(1).mat(:,1:3)));
sesInfo.HInfo.V = data_info;

sesInfo.userInput.pwd = sesInfo.outputDir;
sesInfo.userInput.param_file = [sesInfo.userInput.pwd '\' prefix '_ica_parameter_info.mat'];
sesInfo.userInput.prefix = prefix;
sesInfo.userInput.files = sesInfo.inputFiles;
sesInfo.userInput.diffTimePoints = sesInfo.diffTimePoints;
sesInfo.userInput.numOfSess = sesInfo.numOfSess;
sesInfo.userInput.numOfSub = length(input_folder);
sesInfo.userInput.numOfGroups1 = sesInfo.numOfDataSets;
sesInfo.userInput.numComp = sesInfo.numComp;
sesInfo.userInput.HInfo = sesInfo.HInfo;
sesInfo.userInput.ICA_Options = sesInfo.ICA_Options;
sesInfo.userInput.perfType = 'maximize performance';
sesInfo.userInput.numOfPC1 = ICs;
sesInfo.userInput.backReconType = 'gica';

sesInfo.param_file = strcat(prefix,'_ica_parameter_info');
sesInfo.data_reduction_mat_file = strcat(prefix,'_pca_r');
sesInfo.ica_mat_file = strcat(prefix,'_ica');
sesInfo.back_reconstruction_mat_file = strcat(prefix,'_ica_br');
sesInfo.calibrate_components_mat_file = strcat(prefix,'_ica_c');
sesInfo.aggregate_components_an3_file = strcat(prefix,'_agg__component_ica_');
sesInfo.backReconType = 'gica';
sesInfo.userInput.algorithm = 2;

sesInfo.reduction(1).numOfGroupsBeforeCAT = length(input_folder);
sesInfo.reduction(1).numOfGroupsAfterCAT = sesInfo.reduction.numOfGroupsBeforeCAT;
sesInfo.reduction(1).numOfPCBeforeCAT = sesInfo.diffTimePoints;
sesInfo.reduction(1).numOfPCAfterReduction = ICs;
sesInfo.reduction(1).numOfPrevGroupsInEachNewGroupAfterCAT = ones(sesInfo.numOfScans,1);
sesInfo.reduction(1).numOfPCInEachGroupAfterCAT = sesInfo.diffTimePoints;

sesInfo.reduction(2).numOfGroupsBeforeCAT = length(input_folder);
sesInfo.reduction(2).numOfGroupsAfterCAT = -1;
sesInfo.reduction(2).numOfPCBeforeCAT = ICs;
sesInfo.reduction(2).numOfPCAfterReduction = 0;
sesInfo.reduction(2).numOfPrevGroupsInEachNewGroupAfterCAT = [];
sesInfo.reduction(2).numOfPCInEachGroupAfterCAT = [];

sesInfo.reduction(3).numOfGroupsBeforeCAT = -1;
sesInfo.reduction(3).numOfGroupsAfterCAT = -1;
sesInfo.reduction(3).numOfPCBeforeCAT = 0;
sesInfo.reduction(3).numOfPCAfterReduction = 0;
sesInfo.reduction(3).numOfPrevGroupsInEachNewGroupAfterCAT = [];
sesInfo.reduction(3).numOfPCInEachGroupAfterCAT = [];


end

