function dataReceiverPlotter(src,event)
    persistent tempData;
    persistent tempTime;
    global data 
    global times
    if(isempty(tempData))
         tempData = [];
    end
    if(isempty(tempTime))
         tempTime = [];
    end
    figure(209387);
    plot(event.TimeStamps, event.Data)
    tempData = [tempData;event.Data];
    tempTime = [tempTime;event.TimeStamps];
    data = tempData;
    times = tempTime;
end