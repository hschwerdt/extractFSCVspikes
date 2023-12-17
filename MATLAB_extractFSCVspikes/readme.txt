extractSpikesFSCV(subject,sessnum,varargin) is the main function 
subject = 'patra' for our example data
sessnum = session number (integer)

%%Example Calls to extractSpikesFSCV

extractSpikesFSCV('patra',127,'sites',{'p1','p3','pl3','cl1','cl3'},'nlx')
extractSpikesFSCV('patra',83,'sites',{'p1','p5','pl1','cl3','cl6'},'nlx')