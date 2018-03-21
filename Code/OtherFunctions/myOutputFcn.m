function status = myOutputFcn(t,y,flag, Kp, Ki, T0, r)
% OutputFcn sample
persistent u
switch flag
    case 'init'
        u = plantControllerPI(r, y, Kp, Ki, T0)'; %mass increases linearly with time
    case ''
        for i = 1:numel(t)
            u = [u; plantControllerPI(r, y(:,i), Kp, Ki, T0)'];
        end
    case 'done' % when it's done
        assignin('base','uODE',u); % get the data to the workspace.
end
status = 0;