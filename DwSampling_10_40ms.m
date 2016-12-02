function [TimeData,TimeDataDwsampling] = DwSampling_10_40ms(TimeData12sf,DelayT,SF_num)

before32 = zeros(32,1);                                   % 提前采样的前32个点
after32data = load('after32data.mat');                    % 延时采样的后32个点
after32 = after32data.after32data;
clear after32data;

if (SF_num == 40)
    TimeDataDwsampling = [];
    length12 = size(TimeData12sf,1);
    TimeDataAll = [before32;TimeData12sf;after32 ];
    TimeData12sf = TimeDataAll(33+DelayT:length12+32+DelayT,1);  
    TimeData = TimeData12sf;    
    for i = 0:SF_num/10-1
        TimeData2 = TimeData(307200*i+1:307200*i+30720*10,:); 
        length10 = size(TimeData2,1);
        
        TimeDataTmp = zeros(length10,6);
        filter_dataTmp = zeros(length10,5);
        TimeDataDwsamplingTmp = zeros(length10,5);
        
        TimeDataTmp(:,1) = TimeData2;
        [TimeDataDwsampling_tmp] = DwSampling(TimeDataTmp,length10,10);
        TimeDataDwsampling = [TimeDataDwsampling;TimeDataDwsampling_tmp];
    end;
    clear length12 TimeDataAll TimeData12sf TimeData2 i length10 TimeDataTmp filter_dataTmp TimeDataDwsamplingTmp TimeDataDwsampling_tmp before32 after32;
    
else
    length12 = size(TimeData12sf,1);
    TimeDataAll = [before32;TimeData12sf;after32 ];
    TimeData = TimeDataAll(33+DelayT:length12+32+DelayT,1);   clear TimeDataAll;      % 加入延时的实际采样数据，
    
    TimeDataTmp = zeros(length12,6);
    filter_dataTmp = zeros(length12,5);
    TimeDataDwsamplingTmp = zeros(length12,5);
    
    TimeDataTmp(:,1) = TimeData;    
    [TimeDataDwsampling] = DwSampling(TimeDataTmp,length12,SF_num);
     clear TimeData12sf length12 TimeDataAll  TimeDataTmp filter_dataTmp TimeDataDwsamplingTmp before32 after32; 
end;