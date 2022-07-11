function [options] =  getOptionsdefault()
    options.plot_flag       = false;
    options.PCAflag         = true;
    options.nPC_thresh      = 97;
    options.l               = logspace(-5, +5, 200);
    options.Nperm           = 10000;
    options.save_workspace  = true;
    options.zscoreY         = true;
    options.thr_sbj_overlap = 0;
    options.thr_les_size    = [0 100];
end
