function cfg = config()

cfg.INPUT_DIR   = 'input_sounds';   % Raw sound files from Task 2
cfg.PHASE1_DIR  = 'output_sounds';  % Phase 1 outputs (*_mono16k.wav & Phase 1 plots)
cfg.PHASE2_DIR  = 'output_phase2';  % Phase 2 outputs (Phase 2 plots & phase2_results.mat)

cfg.TARGET_FS   = 16000;            % Required sampling rate, Hz (Task 3.6)
cfg.PLOT_ALL    = true;             % true = plot every file; false = plot only PLOT_FILE
cfg.PLOT_FILE   = '';               % '' = first file found, or e.g. 'male_quiet.wav'

% Phase 1 Settings 
cfg.COS_FREQ    = 1000;            
cfg.PLAY_AUDIO  = true;

% Phase 2 Setting
cfg.N_CHANNELS  = 8;               
cfg.F_MIN       = 100;             
cfg.F_MAX       = 7950;            
cfg.SPACING     = 'log';           
cfg.BP_ORDER    = 4;               
cfg.LPF_CUTOFF  = 400;             
cfg.LPF_ORDER   = 4;               
cfg.SAVE_MAT    = true;            