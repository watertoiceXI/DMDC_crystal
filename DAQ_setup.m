% setup DAQ system
s = daq.createSession('ni');

addAnalogInputChannel(s,'Dev1', 'ai1', 'Voltage'); 
addAnalogInputChannel(s,'Dev1', 'ai2', 'Voltage'); 

s.Rate = 10000;