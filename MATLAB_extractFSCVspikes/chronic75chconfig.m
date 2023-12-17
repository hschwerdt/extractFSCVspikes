sessionnum=75;

ncschannels={'p1','p5','p1-p5','pl1','pl1-p5',...
    'cl1','cl1-cl4','cl3','cl4','cl3-cl4'...
    'cl6','cl4-cl6','cl3-cl6','s5','s4','s5-s4','s3','s4-s3',...
    'eyex','eyey','eyed','lickx','pulse'};    %

letterdrive='Y:';   %COPYHERE
fscvdir='patra_fscv'; %COPYHERE
ephysdir='patra_ephys';
paths{1}=fullfile(letterdrive,'data_MIT',fscvdir,'patra_chronic75_05232018','1dr','cvtotxt',filesep);
paths{2}=fullfile(letterdrive,'data_MIT',ephysdir,'2018-05-23_08-53-18',filesep);
