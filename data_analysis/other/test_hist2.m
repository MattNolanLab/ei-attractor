 events = 1000000;
  x1 = sqrt(0.05)*randn(events,1)-0.5; x2 = sqrt(0.05)*randn(events,1)+0.5;
  y1 = sqrt(0.05)*randn(events,1)+0.5; y2 = sqrt(0.05)*randn(events,1)-0.5;
  x= [x1;x2]; y = [y1;y2];
 
 %For linearly spaced edges:
  xedges = linspace(-1,1,128); yedges = linspace(-1,1,128);
  histmat = hist2(x, y, xedges, yedges);
  figure; pcolor(xedges,yedges,histmat'); colorbar ; axis square tight; 