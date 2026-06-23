% compute bleedthrough and cross-talk
close all
clear all

path='C:\Users\Pr Simon Haziza\Desktop\GitHub\FluorescenceSpectra\SpectraProperties';

% Noise
% fiber_af_488
% fiber_af_561
% tissue_af_488
% tissue_af_561
% HbO
% % Signal
% gevi1_fluo
% ref_fluo
% gevi2_fluo
%
% % Detector
% responsivity_apd
%%
% figure
% for iFile=1:6
% plot(dev3626.demods(1).sample{1, iFile}.frequency,10*log10(dev3626.demods(1).sample{1, iFile}.auxin1pwr))
% hold on
% end
% hold off
%% uSMAART configuration > Filter + ASAP3

folder='Proteins';name='GFP';brightness_gfp=33.5;
% folder='Proteins';name='mNeonGreen';brightness_gfp=92.8;
[g,lambda]=importSpectrum(path,folder,name);

folder='Proteins';name='cyOFP';brightness_cyofp=30.4;
ref=importSpectrum(path,folder,name);

folder='Proteins';name='mRuby3';brightness_mruby2=42.9;
% folder='Proteins';name='mRuby3';brightness_mruby3=57.6;
r=importSpectrum(path,folder,name);

folder='Dichroics';name='Di01-R488_561'; % uSMAART
% folder='Dichroics';name='FF493_574-Di01_v2';
% folder='Dichroics';name='Di03-R488_561';
% folder='Dichroics';name='FF493_574-Di01';
dic1=importSpectrum(path,folder,name);
 
folder='Dichroics';name='FF562-Di03'; % uSMAART
dic2=importSpectrum(path,folder,name);

% folder='Dichroics';name='FF564-Di01';
% dic2=importSpectrum(path,folder,name);

% folder='Filters';name='FF01-525_30';
% folder='Filters';name='FF02-529_24';
folder='Filters';name='FF02-520_28'; % uSMAART
% folder='Filters';name='ET_537_29';
% folder='Filters';name='FF01-531_40';
filtG=importSpectrum(path,folder,name);

% folder='Filters';name='FF01-609_62';
folder='Filters';name='FF01_630_92'; % uSMAART
filtR=importSpectrum(path,folder,name);

% output_greenCh=[g(:,2) r(:,2) ref(:,2) ]'.*(dic1).*(1-dic2).*filtG;
% output_redCh=[g(:,2) r(:,2) ref(:,2) ]'.*(dic1).*(dic2).*filtR;
%
% figure;
% plot(lambda,([output_greenCh+output_redCh]),'linewidth',2)
% legend('green','red','ref')

WL1=488;WL2=561;
output_greenCh_488=[brightness_gfp*g(lambda==WL1,1)*g(:,2)...
                    brightness_cyofp*ref(lambda==WL1,1)*ref(:,2)...
                    brightness_mruby3*r(lambda==WL1,1)*r(:,2)]'.*dic1.*(1-dic2).*filtG;
output_redCh_488=[brightness_gfp*g(lambda==WL1,1)*g(:,2)...
                  brightness_cyofp*ref(lambda==WL1,1)*ref(:,2)...
                  brightness_mruby3*r(lambda==WL1,1)*r(:,2)]'.*dic1.*dic2.*filtR;
output_redCh_561=[brightness_gfp*g(lambda==WL2,1)*g(:,2)...
                  brightness_cyofp*ref(lambda==WL2,1)*ref(:,2)...
                  brightness_mruby3*r(lambda==WL2,1)*r(:,2)]'.*dic1.*dic2.*filtR;

% replace all NaN values by 0
output_greenCh_488(isnan(output_greenCh_488))=0;
output_redCh_488(isnan(output_redCh_488))=0;
output_redCh_561(isnan(output_redCh_561))=0;

figure;
subplot(311)
plot(lambda,([output_greenCh_488]),'linewidth',2)
legend('green','ref','red')
title(['Green Channel @' num2str(WL1)])
% ylim([0 1])
xlim([500 700])

subplot(312)
plot(lambda,([output_redCh_488]),'linewidth',2)
legend('green','ref','red')
title(['Red Channel @' num2str(WL1)])
% ylim([0 1])
xlim([500 700])

subplot(313)
plot(lambda,([output_redCh_561]),'linewidth',2)
legend('green','ref','red')
title(['Red Channel @' num2str(WL2)])
% ylim([0 1])
xlim([500 700])

% normalize per channel assuming 1:1 fluorescent molecule (brightness-aware)
temp = trapz(output_greenCh_488,2);
Raw(:,1)=temp;%Ch_norm(:,1)=temp./max(temp);
temp = trapz(output_redCh_488,2); %Xtalk channel
Raw(:,2)=temp;%Ch_norm(:,2)=temp./max(temp);
temp = trapz(output_redCh_561,2);
Raw(:,3)=temp;%Ch_norm(:,3)=temp./max(temp);

Ch_norm=round(100.*Raw./max(Raw,[],1),1)
Fluor_norm=round(100.*Raw./max(Raw,[],2),1)

% % normalize per GEVI
% temp = trapz(output_greenCh_488(:,500-400+1:700-400+1),2);
% q(:,1,1)=temp;%q(:,1,2)=temp./max(temp);
% temp = trapz(output_redCh_488(:,500-400+1:700-400+1),2);
% q(:,2,1)=temp;%q(:,2,2)=temp./max(temp);
% temp = trapz(output_redCh_561(:,500-400+1:700-400+1),2);
% q(:,3,1)=temp;%q(:,3,2)=temp./max(temp);

% disp('Bleedthrough GFP into Red for dual CT - cyOFP channel contamination (%)')
% 100*Raw(1,2)/Raw(1,1)
% 
% disp('Bleedthrough GFP into Red for single CT - mRuby2 or Varnam2 channel contamination (%)')
% 100*Raw(1,3)/Raw(1,1)
% 
% disp('Bleedthrough cyOFP into Green for dual CT (%)')
% 100*Raw(2,1)/Raw(2,2)
% 
% disp('Bleedthrough mRuby into Green for single CT (%)')
% 100*Raw(3,1)/Raw(3,3)

% idx=@(L) L-400+1;
% greenCh_488BR=log10(100*0.05.*dic1(idx(488)).*(1-dic2(idx(488))).*filtG(idx(488)));
% greenCh_561BR=log10(100*0.05.*dic1(idx(561)).*(1-dic2(idx(561))).*filtG(idx(561)));
% redCh_488BR=log10(100*0.05.*dic1(idx(488)).*(dic2(idx(488))).*filtR(idx(488)));
% redCh_561BR=log10(100*0.05.*dic1(idx(561)).*(dic2(idx(561))).*filtR(idx(561)));
% 
% % disp([q(1,1,1) q(3,1,2) ])
% disp([q(:,:,2)])
% disp(round([greenCh_488BR greenCh_561BR redCh_488BR redCh_561BR],1));

% disp(q)
%% TEMPO configuration > Filter + Ace1

folder='Proteins';name='mNeonGreen';brightness_mneongreen=92.8;
[g,lambda]=importSpectrum(path,folder,name);

folder='Proteins';name='mRuby3';
r=importSpectrum(path,folder,name);

folder='Proteins';name='cyOFP';
ref=importSpectrum(path,folder,name);

folder='Dichroics';name='FF493_574-Di01';
dic1=importSpectrum(path,folder,name);

folder='Dichroics';name='FF564-Di01';
dic2=importSpectrum(path,folder,name);

folder='Filters';name='ET_537_29';
filtG=importSpectrum(path,folder,name);

folder='Filters';name='FF01_630_92';
filtR=importSpectrum(path,folder,name);

% output_greenCh=[g(:,2) r(:,2) ref(:,2) ]'.*(dic1).*(1-dic2).*filtG;
% output_redCh=[g(:,2) r(:,2) ref(:,2) ]'.*(dic1).*(dic2).*filtR;

% figure;
% plot(lambda,log10([output_greenCh+output_redCh]),'linewidth',2)
% legend('green','red','ref')

output_greenCh_488=[g(lambda==488,1)*g(:,2) r(lambda==488,1)*r(:,2) ref(lambda==488,1)*ref(:,2) ]'.*(dic1).*(1-dic2).*filtG;
output_redCh_488=[g(lambda==488,1)*g(:,2) r(lambda==488,1)*r(:,2) ref(lambda==488,1)*ref(:,2) ]'.*(dic1).*(dic2).*filtR;
output_redCh_561=[g(lambda==561,1)*g(:,2) r(lambda==561,1)*r(:,2) ref(lambda==561,1)*ref(:,2) ]'.*(dic1).*(dic2).*filtR;

figure;
subplot(311)
plot(lambda,([output_greenCh_488]),'linewidth',2)
legend('green','red','ref')
title('Green Channel @488')
ylim([0 1])
xlim([500 700])

subplot(312)
plot(lambda,([output_redCh_488]),'linewidth',2)
legend('green','red','ref')
title('Red Channel @488')
ylim([0 1])
xlim([500 700])

subplot(313)
plot(lambda,([output_redCh_561]),'linewidth',2)
legend('green','red','ref')
title('Red Channel @561')
ylim([0 1])
xlim([500 700])

temp = trapz(output_greenCh_488(:,500-400+1:700-400+1),2);
q(:,1)=temp./max(temp);
temp = trapz(output_redCh_488(:,500-400+1:700-400+1),2);
q(:,2)=temp./max(temp);
temp = trapz(output_redCh_561(:,500-400+1:700-400+1),2);
q(:,3)=temp./max(temp);
q









