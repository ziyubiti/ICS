


Temp2 = TimeDataDwsampling(960*3+1:960*4,1);
TimeDataDwsampling2 = TimeDataDwsampling;%output;
 figure;plot(real(Temp2)/norm(Temp2)*sqrt(960));grid on;
figure;plot(real(TimeDataDwsampling2)/norm(Temp2)*sqrt(960));grid on;

TimeDataDwsampling2 = TimeDataDwsampling(960+1:960*2,1);

  for i = 1:3
 xcor(i,:) = circle_corr(TimeDataDwsampling2.',time_seq(i,:));
      xcor_num22(i,:) = xcor(i,:);
      Pxcor_num22(i,:) = abs(xcor_num22(i,:)).*abs(xcor_num22(i,:));      % power
        P2xcor_num22 = Pxcor_num22;                   % plot and search
        

       
         if (strcmp(cs_para.pss_noise_method,'PART'))
             [peak_pss(1,i), pos(1,i)] = max(Pxcor_num22(i,:));% 1st peak
             Pxcor_num22(i,pos(1,i)) = 0;
             [peak_pss(2,i), pos(2,i)] = max(Pxcor_num22(i,:));% 2st peak
             Pxcor_num22(i,pos(2,i)) = 0;
             noise_power_pss(i,1) = mean(Pxcor_num22(i,:));
             Ppss_th(i,1) = cs_para.pss_threshold*noise_power_pss(i,1);
         elseif(strcmp(cs_para.pss_noise_method,'ALL'))
                 noise_power_pss(i,1) = mean(Pxcor_num22(i,:));
                 Ppss_th(i,1) = cs_para.pss_threshold*noise_power_pss(i,1);
         end;
         
%          markpeak = 1103;
%          P2xcor_num22(:,[1:markpeak-8 markpeak+8:end]) = 0;
         
         
         figure(5*i);
         stem(P2xcor_num22(i,:));      grid on; hold on;
         plot([0 size(P2xcor_num22,2)],[Ppss_th(i,1) Ppss_th(i,1)],'r',[0 size(P2xcor_num22,2)],[noise_power_pss(i,1) noise_power_pss(i,1)],'c');
         ylabel('Pxcor combine22');
 
 
 
 
 
  end;
  
