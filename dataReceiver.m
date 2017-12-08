function dataReceiver(src,event)
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
    tempData = [tempData;event.Data];
    tempTime = [tempTime;event.TimeStamps];
    data = tempData;
    times = tempTime;
end