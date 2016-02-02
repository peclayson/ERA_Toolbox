function ra_start

% Construct a questdlg with three options
questdlg('What would you like to do?', ...
	'Reliability Analyses', ...
	'Process new data','View results','Exit','Exit');
% Handle response
switch task
    case 'Process new data'
        ra_startproc;
    case 'View results'
        ra_startview;
    case 'Exit'
        return;
end

end
