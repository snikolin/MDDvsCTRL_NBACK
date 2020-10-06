function comps = MDDvsCTRL_ICAcheck(comps)
% Identify which components are to be removed
% Done through automated process with follow-up visual confirmation
% Automated process to suggest components for removal
% At the prompt, list the components to reject
% Example: [3,6,7,8,9];

%% Automated ICA inspection
[blinks, eye_movements, muscle, gen_disc, suggested_comps,all_art] = MDDvsCTRL_ICAauto(comps);

%% Visual ICA inspection
MDDvsCTRL_ICAplot(comps,suggested_comps);

%% List components for rejection
% Suggestion of components to remove
disp('Automated ICA reject components are:')
disp(['Blinks: ', num2str(blinks)]);
disp(['Lateral eye movements: ', num2str(eye_movements)]);
disp(['Generic discontinuities: ', num2str(gen_disc)]);
disp(['Muscle: ', num2str(muscle)]);
disp(['All artifacts: ', num2str(all_art)]);
disp('.');
disp(['Suggested: ', num2str(suggested_comps)]);

% Creates prompt for user to include which components need to be rejected
prompt          = 'list components to reject. Approximately 10-15% of total:';
x               = input(prompt);
comps.rejected  = x;
close
