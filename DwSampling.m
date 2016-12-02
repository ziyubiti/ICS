% input:  LenData;
% output: SF_num Dwsampling Data; 30720/cs_para.Ratio*SF_num


function [TimeDataDwsampling] = DwSampling(TimeDataTmp,LenData,SF_num)

global cs_para;

length12 = LenData;


if strcmp(cs_para.filter_method, 'DSP_HB' )
    
    coee = cs_para.coee65;
    Nstage = log2(cs_para.Ratio);
    
    %     Hhb = dfilt.df2t(coee);
    %     Hcas = dfilt.cascade(Hhb,Hhb,Hhb,Hhb,Hhb);
    %     hcas = coeffs(Hcas);
    %     freqz(hcas);
    %     He = convert(Hcas,'df2t');
    
    if(length(coee) == 33)
        filterconvTmp = zeros(30752,5);
        for k = 1:5                                                % 5级滤波，每次2倍降采样
            filterconvTmp(1:30720/(2.^(k-1))+32,k) = conv(coee,TimeDataTmp(1:30720/(2.^(k-1)),k));            % 30720+33-1
            filter_dataTmp(1:30720/(2.^(k-1)),k) = filterconvTmp(17:30720/(2.^(k-1))+16,k);                   % 每次conv 延时（33-1）/2=16
            TimeDataDwsamplingTmp(1:30720/(2.^k),k) = filter_dataTmp(1:2:30720/(2.^(k-1)),k);
            TimeDataTmp(1:30720/(2.^k),k+1) = TimeDataDwsamplingTmp(1:30720/(2.^k),k);
        end;
    else       % 65阶
        filterconvTmp = zeros(length12+64,5);
        for k = 1:Nstage                                                % 5级滤波，每次2倍降采样
            filterconvTmp(1:length12/(2.^(k-1))+64,k) = conv(coee,TimeDataTmp(1:length12/(2.^(k-1)),k));            % 30720+65-1
            filter_dataTmp(1:length12/(2.^(k-1)),k) = filterconvTmp(33:length12/(2.^(k-1))+32,k);                   % 每次conv 延时（65-1）/2=32;   不考虑延时，取最前面的data长数；
            TimeDataDwsamplingTmp(1:length12/(2.^k),k) = filter_dataTmp(1:2:length12/(2.^(k-1)),k);
            TimeDataTmp(1:length12/(2.^k),k+1) = TimeDataDwsamplingTmp(1:length12/(2.^k),k);
        end;
    end;
    
    
    
    %         %filter_data = filter(coee,1,TimeData);                    % 30720
    %         filterconv = conv(coee,TimeData);                         % 30720+127-1
    %         filter_data = filterconv(17:30736,:);                     % 每次conv 延时（33-1）/2=16
    TimeDataDwsampling(:,1) = TimeDataTmp(1:30720/cs_para.Ratio*SF_num,Nstage+1);           % 1920;   960*10 =9600 ,960*12 =11520
    
    
elseif strcmp(cs_para.filter_method, '1FIR_R32' )
    coee = cs_para.coee737;
    filterconvTmp = zeros(length12+size(coee)-1,1);
    filterconvTmp(1:length12+size(coee)-1,1) = conv(coee,TimeDataTmp(1:length12,1));
    filter_dataTmp(1:length12,1) = filterconvTmp((size(coee)-1)/2+1:length12+(size(coee)-1)/2,1);
    
    TimeDataDwsampling(:,1) = filter_dataTmp(1:cs_para.Ratio:30720*SF_num,1);
    
elseif strcmp(cs_para.filter_method, '4HB_1FIR' )
    coee1 = cs_para.coee65;
    coee2 = cs_para.coee49;
    for k = 1:4                                                % 4级HB滤波，每次2倍降采样
        filterconvTmp(1:length12/(2.^(k-1))+64,k) = conv(coee1,TimeDataTmp(1:length12/(2.^(k-1)),k));            % 30720+65-1
        filter_dataTmp(1:length12/(2.^(k-1)),k) = filterconvTmp(33:length12/(2.^(k-1))+32,k);                   % 每次conv 延时（65-1）/2=32
        TimeDataDwsamplingTmp(1:length12/(2.^k),k) = filter_dataTmp(1:2:length12/(2.^(k-1)),k);
        TimeDataTmp(1:length12/(2.^k),k+1) = TimeDataDwsamplingTmp(1:length12/(2.^k),k);
    end;
    k = k + 1;
    filterconvTmp(1:length12/(2.^(k-1))+size(coee2)-1,k) = conv(coee2,TimeDataTmp(1:length12/(2.^(k-1)),k));
    filter_dataTmp(1:length12/(2.^(k-1)),k) = filterconvTmp((size(coee2)-1)/2+1:length12/(2.^(k-1))+(size(coee2)-1)/2,k);
    
    TimeDataDwsampling(1:length12/(2.^k),1) = filter_dataTmp(1:2:length12/(2.^(k-1)),k);
    
    
    
end;