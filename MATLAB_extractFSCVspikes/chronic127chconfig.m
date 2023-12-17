sessionnum=127;

ncschannels={'p1','p3','p1-p3','pl3',...
    'p1-pl3','p3-pl3',...
    'cl1','cl3','cl1-cl3','cl4','cl1-cl4','cl3-cl4',...
    'cl5','cl4-cl5','cl1-cl5','s5','s4','s5-s4','s3','s4-s3',...
    'eyex','eyed','lickx','pulse'};    %all working channels

letterdrive='Y:';   
fscvdir='patra_fscv';   %FSCV folder
ephysdir='patra_ephys'; %Corresponding Ephys folder
paths{1}=fullfile(letterdrive,'data_MIT',fscvdir,'patra_chronic127_08142018','1dr','cvtotxt',filesep);
paths{2}=fullfile(letterdrive,'data_MIT',ephysdir,'2018-08-14_09-52-42',filesep);

