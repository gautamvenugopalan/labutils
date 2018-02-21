%%Wrapper function to plot all the files in a given directory
%%Can handle spectra or TFs taken with SR785 or AG analyzers
%%Calls importdata40m, outputs a structure with handles to plots

function [h,data,filenames] = plotData40m(dirname,analyzer,dataType,nPts)

plotFlag=0;

filenames = dir(strcat(dirname,'*.txt'));

data = zeros(length(filenames),nPts,3);
for i=1:length(filenames)
[a,b,c] = importdata40m(filenames(i).name,analyzer,dataType,plotFlag);
for j=1:length(a)
data(i,j,1) = a(j);
data(i,j,2) = b(j);
    if strcmp(dataType,'TF') == 1
        data(i,j,3) = c(j);
    end
end
end

%At this point, data should be loaded into the array 'data' - first index
%denotes file, second index denotes datapoint #, third index denotes
%f/mag/phase

%Now actually plot the data

if strcmp(dataType,'spec') == 1
    h = zeros(length(filenames),1);
    figure
    hold on
    for i=1:length(filenames)
        f = data(i,:,1);
        mag = data(i,:,2);
        h(i)=loglog(f,mag,'LineWidth',3);
    end
    grid on
    xlabel('Frequency (Hz)','FontSize',16)
    ylabel('ASD (V/rtHz)','FontSize',16)
    set(gca,'FontSize',16)
    set(gca,'XScale','log')
    set(gca,'YScale','log') 
    axis tight
    hold off
    
elseif strcmp(dataType,'TF') == 1
    h = zeros(length(filenames),2);
    figure
    hold on
    subplot(2,1,1)
    hold on
    subplot(2,1,2)
    hold on
    for i=1:length(filenames)
        f = data(i,:,1);
        mag = data(i,:,2);
        phase = data(i,:,3);
        subplot(2,1,1)
        h(i,1)=semilogx(f,mag,'LineWidth',3);
        subplot(2,1,2)
        h(i,2)=semilogx(f,phase,'LineWidth',3);
    end
    subplot(2,1,1)
    grid on
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Magnitude (dB)','FontSize',14)
    set(gca,'FontSize',14)
    axis tight
    set(gca,'XScale','log')
    hold off
    
    subplot(2,1,2)
    hold on
    grid on
    xlabel('Frequency (Hz)','FontSize',14)
    ylabel('Phase (degrees)','FontSize',14)
    set(gca,'FontSize',14)
    set(gca,'YTick',[-180:45:180]);
    set(gca,'XScale','log')
    axis tight
    hold off
%     legend(h(:,2),filenames.name)
end


end
