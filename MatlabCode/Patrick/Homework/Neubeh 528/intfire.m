% Basic integrate-and-fire neuron 
% R Rao 2007

%clear
% input current
I = 1 % nA
fmax = 1000

% Variables
Ifinal = 10
deltaI = .01
I_trace = []; % current trace for plotting
FR = []
%freqs = FR(FR~=0)

% Cycle through input currents
for f = freq %1:fmax
    
    % capacitance and leak resistance
    C = 1 % nF
    R = 40 % M ohms
    
    % I & F implementation dV/dt = - V/RC + I/C
    % Using h = 1 ms step size, Euler method
    
    V = 0;
    tstop = 200;
    abs_ref = 5; % absolute refractory period
    ref = 0; % absolute refractory period counter
    V_trace = []; % voltage trace for plotting
    V_th = 10; % spike threshold
    
    for t = 1:tstop
        
        I = sin(t)*f;

        if ~ref
            V = V - (V/(R*C)) + (I/C);
        else
            ref = ref - 1;
            V = 0.2*V_th; % reset voltage
        end
        
        if (V > V_th)
            V = 50;  % emit spike
            ref = abs_ref; % set refractory counter
        end
        
        V_trace = [V_trace V];
        
        
    end
    FR = [FR sum(V_trace==50)];
    I_trace = [I_trace I];
    %I = I+deltaI;
figure; hold on;
plot(sin(1:tstop).*f)
plot(V_trace)
end


figure(1)
plot(V_trace)
xlabel('Time (ms)')
ylabel('Voltage (mV)')

% figure(2)
% plot(1:fmax,FR)
% title('Input Frequency vs Firing Rate')
% xlabel('Input Frequency')
% ylabel('Firing Rate (sp/s)')

