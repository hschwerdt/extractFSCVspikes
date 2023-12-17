function [lfpchs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map,varargin)
%get csc ids based on channel names that we used to identify
lfpchs=[];
bipolarchs={};
otherchs=[];
sumchs=[];
zid=[];
bcount=0;
lcount=0;
ocount=0;
argnum=1;
cflag=0;
while argnum <= length(varargin)
    switch varargin{argnum}       
        case 'sumchannels'
            %channels to sum (ie. lick channels if not default 39:41)
            argnum=argnum+1;
            sumchs=varargin{argnum};
        case 'cleolick37'
            %cleo lick ch only 37, no sum ch
            cflag=1;
        case ''
            %ignore
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end


for id=1:length(ncschannels)
    xid=find(ismember(csc_map,ncschannels{id})==1);
    xlab=xid-1;     %label with csc values
    cscnum=str2num(csc_map{xlab});
    %if # is > 100, must be bipolar config where csc num separated by 0
    if cscnum>100
        zids=strfind(csc_map{xlab},'0');
        %if more than 1 zero in label
        %if next to each other then the final one that is not end of label
        if length(zids)>1
            if zids(2)==zids(1)+1
                %if next to each other, then 2nd one
                zid=zids(2);
            end
            if zids(2)==length(csc_map{xlab})
                %if at end, then must be first one
                zid=zids(1);
            end
        else
            %just single zero
            zid=zids;
        end
        bcount=bcount+1;
        bipolarchs{bcount}=[str2num(csc_map{xlab}(1:zid-1)),...
            str2num(csc_map{xlab}(zid+1:end))];   
    else 
        if cscnum<33
                %then recorded from 1st headstage, lfp ch
                lcount=lcount+1;
                lfpchs(lcount)=cscnum;
        end
        if cscnum>=33
            %then physiological signal from 2nd headstage
            ocount=ocount+1;
            otherchs(ocount)=cscnum;
        end
    end
end

if cflag==1
    %cleo lick on ch37
    sumchs=37;
end
if isempty(sumchs) && (ismember('lick', ncschannels) || ismember('lickx',ncschannels))
    %default lick channels except for early sesions need to supply varargin
    sumchs=[39:41];    
    %make sure these are also included in otherchs
    otherchstemp=otherchs;
    nonsum=find(~ismember(otherchs,sumchs)==1);
    sumid=find(ismember(otherchs,sumchs)==1);
    otherchs(sumid:sumid+length(sumchs)-1)=sumchs;
    otherchs=[otherchs otherchstemp(sumid+1:end)];
end

end
        
