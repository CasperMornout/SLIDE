function [] = micronbar1(length, bgAlpha, bgColor, FontSize, FontColor)

% micron bar settings
f=gcm;
mp=getappdata(f.currentAxes,'mapPlot');
mp.micronBar.scanUnit='um'; % change scan unit
mp.micronBar.length=length; % change bar length
mp.micronBar.backgroundAlpha=bgAlpha; % change transparency
mp.micronBar.backgroundColor=bgColor; % change background color
mp.micronBar.lineColor='w'; % change bar color
mp.micronBar.txt.FontSize=FontSize; % change text font size
mp.micronBar.txt.Color=FontColor; % change text font color
mp.micronBar.txt.VerticalAlignment= 'cap'; % change text location
mp.micronBar.txt.String = ['\rm{\textbf{\boldmath{' mp.micronBar.txt.String '}}}'];
end