%% General Parameters
function [videoData,videoTime,current,voltage,times2] = LC_run(runtime,DAQ,Camera,RPi,show)
global data
global times

if show
    clear dataReceiverPlotter
else
    clear dataReceiver
end
% create our clean up object
cleanupObj = onCleanup(@cleanMeUp);
 
%RPi = false;
%DAQ = true;
%Camera = false;
%show = true; %if true, displays some of the data at the end. There will aways be a real-time show.
%runtime = 10;

% the below parameters assume a 1khz driving signal. 
ml_sigma = 20000; %wavelet width. Bigger means more precision, but less speed.
ml_len = 1000; %wavelet length. Should likely be set from above, but it's here
ml_period = 10; 
v_cal = 0.062725918672321;
c_cal = 0.035365538339312;

laser_pulse_rate = 1; % in seconds


%% Setup functions and feature-specific parameters
if RPi
    rpi_ip = '169.254.54.174';
    rpi_user = 'admin';
    rpi_passw = 'RPi144';
    rpi_pin_cmd = 'gpio write 29 ';
    rpi_pin_setup = 'gpio mode 29 out';
    fprintf('Setting up rpi GPIO for laser\n');
end

if Camera
    fprintf('Starting video\n');
    vid = videoinput('pointgrey');%, 1, 'MONO16_1392x1040');
    xPos =0;%249;
    yPos = 0;%327;
    width = 2592;%900;%518;400
    height = 1944;%   
    framerate = 10; %10 fps, as set in the flycap2 software. Can be as high as 13 at this
                    % resolution, but 10 is a nice number. If you change it
                    % in the flycap software, change this here!
    totalframes = runtime*framerate;
    vid.ROIPosition = [xPos yPos width height];
    vid.FramesPerTrigger = totalframes;
end

if DAQ
    x = 1:ml_len;
    ml = exp(-(x-(ml_len/2)).^2/(2*ml_sigma)).*exp(1i*2*pi*x/ml_period);
    fprintf('If needed, setting up the daq board\n');
    if ~exist('s','var')
        DAQ_setup();
    end
    fprintf('Starting daq\n')
    s.DurationInSeconds = runtime;
    if show
        lh = addlistener(s,'DataAvailable',@dataReceiverPlotter);
    else
        lh = addlistener(s,'DataAvailable',@dataReceiver); 
    end
end


%% Starting all the requested devices:
if RPi
    ssh2_simple_command(rpi_ip,rpi_user,rpi_passw,rpi_pin_setup)
    toggle = '1';
end

if Camera
    start(vid);
end

if DAQ
    startBackground(s);
end

%% Gathering Data
tic;

if RPi
    tdot = [];
    numpulse = ceil(laser_pulse_rate/runtime);
    i=0;
    while i < numpulse
    ssh2_simple_command(rpi_ip,rpi_user,rpi_passw,[rpi_pin_cmd toggle]);
    t = toc;
    tdot = [tdot;t];
    if strcmp(toggle,'1')
        toggle='0';
    else
        toggle='1';
    end
    i=i+1;
    pause(laser_pulse_rate);
    end
end

if DAQ
    while(~s.IsDone)
    %fprintf('.')
    t = toc;
    if t>(1.5*runtime)
        break;
    end
    pause(1)
    end
end
fprintf('\n')
close(gcf);

if Camera
    while toc<(runtime)
        pause(1);
    end
    [videoData,videoTime] = getdata(vid);
else
    videoData = 0;
    videoTime = 0;
end
%figure(24378);
%plot(times,data); %  plot global data
%[times,i] = sort(times);
%data = data(i,:);
%hold on;
%plot(tdot,zeros(length(tdot)),'r*');
if DAQ
    voltage = v_cal*abs(conv(data(:,1),ml,'valid'));
    current = c_cal*abs(conv(data(:,2),ml,'valid'));
    %time = 0:0.0001:(sec_len-0.0001);
    vraw = conv(data(:,1),ml,'valid');
    craw = conv(data(:,2),ml,'valid');
    vphase = atand(real(vraw)./imag(vraw));
    cphase = atand(real(craw)./imag(craw));
    pdif = cphase-vphase;
else
    voltage=0;
    current=0;
    pdif = 0;
end
if show
    figure(78354);
    clf();
    hold on;
    plot(times(500:(length(times)-500)),abs(current).*abs(voltage),'bo-','MarkerSize',1);
    %plot(times(500:(length(times)-500)),abs(voltage),'ro-','MarkerSize',1); 
    %plot(times(500:(length(times)-500)),abs(cphase-vphase),'go','MarkerSize',1);
    %plot(time,scaled_switch,'y');
    title(['Current: ' num2str(mean(current)) 'u-amps Voltage: ' num2str(mean(voltage)) 'V']);
    xlabel([num2str(mean(abs(vphase-cphase))) 'deg phase difference'])
    %plot(tdot,mean(abs(current).*abs(voltage))*ones(length(tdot)),'r*');
end
times2 = times;

if Camera && show
    implay(videoData,10)
end

function cleanMeUp()

end

end 