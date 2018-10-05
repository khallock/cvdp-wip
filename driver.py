import os
import datetime
import sys
#import ncl_scripts.functions as ncl
import subprocess
import glob
import shutil
import time

# ============================================================================================
#
# CVDP driver script. To run the CVDP at the command line type: python driver.py
# To run the CVDP at the command line, put it in background mode, and write the terminal output
# to a file named file.out, type: python driver.py >&! file.out &
#
# ============================================================================================

class CVDPDriver:
    # ======== BEGIN USER MODIFICATIONS ===========================================================

    outdir            = './cvdp_output/'   # location of output files   (must end in a '/')
    # It is recommended that a new or empty directory be pointed to here
    # as existing files in outdir can get removed.

    namelists_only       = False               # Set to True to only create the variable namelists. Useful
    # upon running the package for the first time to verify that the correct
    # files are being selected by the package. (See files in namelist_byvar/ directory)                                               # Set to False to run the entire package.

    obs                  = True                # True = analyze and plot observations (specified in namelist_obs), False = do not
    scale_timeseries     = False               # True = scale timeseries so that x-axis length is comparable across timeseries, False = do not
    output_data          = True                # True = output selected calculated data to a netCDF file. Make sure .nc files from previous CVDP
    #        runs are not in outdir or they will get added to or modified.
    compute_modes_mon    = True                # True = compute DJF, MAM, JJA, SON, Annual and Monthly Atmospheric Modes of Variability
    # False = do not compute the Monthly Atmospheric Modes of Variability  (saves computation time)
    # - - - - - - - - - - - - - - - - - -
    opt_climo         = 'Full'                 # Full  = remove climatology based on full record of each simulation,
    # Custom = set climatological period using climo_syear (climatological start year) and climo_eyear (climatological end year)

    # use these values if opt_climo == 'Custom'.   When climo_syear and climo_eyear are positive, remove the climatology/annual cycle based on these years.
    climo_syear    = -30                  #   Both settings should be within the range of years of all specified model runs and observational datasets.
    climo_eyear    = 0                    # When climo_syear is negative, remove the climatology/annual cycle relative to the end of each model run
                                          #   being removed from the last 26 years of each model run and observations.
    # - - - - - - - - - - - - - - - - - -
    colormap          = 0               # 0 = default colormaps, 1 = colormaps better for color blindness

    output_type       = 'png'           # png = create png files, ps = create postscript files as well as png files (for web viewing).

    png_scale         = 1.5             # Set the output .png size. Value between .1->5.  Any value > 1 (< 1) increases (decreases) png size.
    # When output_type = 'png' a value of 1 will result in a png sized 1500 (H) x 1500 (W) before automatic cropping of white space
    # When output_type = 'ps'  a value of 1 will result in a png density setting of 144 before automatic cropping of white space

    webpage_title     = 'Title goes here'          # Set webpage title

    tar_output        = False                 # True = tar up all output in outdir and remove individual files, False = do not
    # Note: ALL files in outdir will be tarred up and then removed from the outdir directory.

    # ---Advanced Options----------------------------------------------------------------------
    zp = './py_scripts/'    # directory path of CVDP python scripts. (must end in a '/')
    # Examples: 'py_scripts/' if all code is local, or on CGD or CISL systems: '~asphilli/CESM-diagnostics/CVDP/Release/v4.1.0/py_scripts/'
    # Regardless of this setting the following files should be in one directory: namelist, driver.ncl, and namelist_obs.
    # If pointing to code in ~asphilli make sure the driver script version #s match between this script and the script in ~asphilli.

    ncl_exec = 'ncl'       # This can be changed to a different path if a different version of NCL needs to be used, such as '/different/path/to/bin/ncl'

    run_style = 'parallel'   # parallel = allow simple python-based parallelization to occur. X number of CVDP NCL scripts will be called at once.
    #            X is set via max_num_tasks. Terminal output will be harder to follow.
    # serial = call CVDP NCL scripts serially. (Default)

    modular = 'True'      # True = Run only those CVDP scripts specified in modular_list.
    # False = Run all CVDP scripts (Default)

    modular_list = ['tas.trends_timeseries']  # When modular = 'True' list the CVDP scripts that will be run.
    # Example: modular_list = 'amoc,amo,pr.trends_timeseries'
    # For a list of available scripts see complete_list at line 72.

    machine_casesen = 'True'  # True = Your filesystem is case sensitive  (Default)
    # False = Your filesystem is case insensitive

    #
    # user-settings associated with task-parallel mode of running
    #
    poll_interval = 1.  # seconds  SET TO 1 FOR DEBUGGING   POLL_INTERVAL = 15.  # seconds

    # max number of CVDP NCL scripts can be run in parallel is CVDP configured for task-parallelism
    max_concurrent = 3  # if in doubt, try 3
    # ========END  MODIFICATIONS===========================================================

    version = '5.0.0'

    complete_list = ['psl.nam_nao', 'psl.pna_npo', 'tas.trends_timeseries', 'snd.trends', 'psl.trends', 'amo', 'pdo',
                     'sst.indices', 'pr.trends_timeseries', 'psl.sam_psa', 'sst.mean_stddev', 'psl.mean_stddev',
                     'pr.mean_stddev', 'sst.trends_timeseries', 'amoc', 'tas.mean_stddev', 'snd.mean_stddev',
                     'aice.mean_stddev', 'aice.trends_timeseries,ipo']

    outfiles = ('ts', 'trefht', 'psl', 'prect', 'snowdp', 'moc', 'maxnum', 'aice_nh', 'aice_sh')


    def run(self):
        print(f'Starting: Climate Variability Diagnostics Package ({datetime.datetime.now().strftime("%c")})')

        for gg in self.outfiles:                 # line 77 NCL script
            path = 'obs_' + gg
            if os.path.exists(path):
                os.remove(path)

        path = self.outdir + 'metrics_orig.txt'
        if os.path.exists(path):        # remove metrics_orig.txt file if present
            os.remove(path)

        if self.opt_climo == 'Custom':       # line 88 NCL script
            if self.climo_syear >= self.climo_eyear:
                print('Specified custom climatology start year (climo_syear) cannot be greater than or equal to the specified end year (climo_eyear), exiting CVDP.')
                quit()
        else:
            self.climo_syear = -999
            self.climo_eyear = -999

        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)

        myEnvVars = {'OUTDIR' : self.outdir, 'OBS' : str(self.obs), 'SCALE_TIMESERIES' : str(self.scale_timeseries),
                     'OUTPUT_DATA' : str(self.output_data), 'VERSION' : self.version, 'PNG_SCALE' : str(self.png_scale),
                     'OPT_CLIMO' : self.opt_climo, 'CLIMO_SYEAR' : str(self.climo_syear), 'CLIMO_EYEAR' : str(self.climo_eyear),
                     'COMPUTE_MODES_MON' : str(self.compute_modes_mon), 'OUTPUT_TYPE' : self.output_type,
                     ### CORRECT ONE    'MACHINE' : str(self.machine_casesen), 'COLORMAP' : str(self.colormap), 'CVDP_SCRIPTS' : self.zp
                     'MACHINE': str(self.machine_casesen), 'COLORMAP': str(self.colormap), 'CVDP_SCRIPTS': './ncl_scripts/'
                     }

        # theEnv = os.environ
        # theEnv.update(myEnvVars)

        # inject environment variables into the running process; subprocesses will inherit these from the parent
        for envVar,envVal in myEnvVars.items():
            os.putenv(envVar, envVal)

        # get the name of the python interpreter -- this will preserve any conda environment paths
        exec_str = sys.executable

        p = subprocess.run(f'{exec_str}  {self.zp}/namelist.py', shell=True)     # create variable namelists

        if self.namelists_only:
            print('Variable namelists have been created. Examine files in namelist_byvar/ directory to verify CVDP file selection.')
            print('Finished: Climate Variability Diagnostics Package ({datetime.datetime.now().strftime("%c")})')
            for gg in self.outfiles:
                path = 'obs_' + gg
                if os.path.exists(path):
                    os.remove(path)

        if self.modular:       # line 120 NCL script
            modular_list = [ml.strip() for ml in self.modular_list]
            # ### THIS IS THE RIGHT ONE     modular_list = [self.zp + ml + '.ncl' for ml in self.modular_list]
            self.modular_list = ['./ncl_scripts/' + ml + '.ncl' for ml in self.modular_list]      #THIS ONE FOR TESTING
            if self.run_style == 'parallel':
                self.manageTasks(self.modular_list)
            else:
                for script in self.modular_list:
                    subprocess.run([self.ncl_exec, script], stdin=subprocess.DEVNULL)
        else:
            self.complete_list = [self.zp + cl + '.ncl' for cl in self.complete_list]
            if self.run_style == 'parallel':
                self.manageTasks(self.complete_list)
            else:
                for script in self.complete_list:
                    subprocess.run([self.ncl_exec, script], stdin=subprocess.DEVNULL)

        if self.output_data:       # line 145 NCL script
            print('skipping ncfiles.append -- requires NetCDF operators')
            # ##p = subprocess.run([ncl_exec, './ncl_scripts/'+'ncfiles.append.ncl']) ####subprocess.run([ncl_exec, zp+'ncfiles.append.ncl'], \
                 ## stdin=subprocess.DEVNULL )

        if self.output_type == 'png':
            ofiles = glob.glob(self.outdir + '*.png')
            for gg in ofiles:
                p = subprocess.run(['convert','-trim','+repage','-border','8','-bordercolor','white',gg,gg])
                print(p.args)  # ##RLB
        else:
            ofilesS = glob.glob(self.outdir + '*.' + self.output_type)
            for gg in ofilesS:
                filesize = os.path.size(gg)
                if filesize < 10000:
                    print('Removing: ' + gg)
                    os.remove(gg)

            ofiles = glob.glob(self.outdir+'*.'+self.output_type)
            ofiles_png = ofiles.replace('.'+self.output_type,'.png')
            d_opt = 144*self.png_scale
            print(f'Converting {self.output_type} files to .png')
            for gg in ofiles:
                p = subprocess.run(['convert','-density',d_opt,'-trim','+repage','-border','8',
                                    'bordercolor','white','-flatten',gg,gg])
            print(f'Done with {self.output_type}->png conversion')

        shutil.copy('./ncl_scripts/' + 'cas-cvdp.png', self.outdir) # ##shutil.copyfile(zp+'cas-cvdp.png', outdir)
        byvarFiles = glob.glob('namelist_byvar/*')
        for bv in byvarFiles:
            shutil.copy(bv, self.outdir)
        shutil.copy('namelist', self.outdir)
        if self.obs:
            shutil.copy('namelist_obs', self.outdir)

        # #########TEMPORARY FOR TESTING#########
        self.zp = './ncl_scripts/'

        met_files = glob.glob(self.outdir+'metrics.*.txt')      # line 183 NCL script
        if len(met_files) == 9:         # all 9 metrics text files are present, create metrics table(s)
            p = subprocess.run([self.ncl_exec, self.zp+'metrics.ncl'])

        p = subprocess.run([self.ncl_exec, 'webtitle="'+self.webpage_title+'"', self.zp+'webpage.ncl'],
                           stdin=subprocess.DEVNULL)

        # -------------------------------
        if self.tar_output:       # line 194 NCL script
            if os.path.exists(self.outdir+'cvdp.tar'):
                os.remove(+self.outdir+'cvdp.tar')
            p = subprocess.run('cd '+self.outdir+'# tar -cf cvdp.tar *')
            p = subprocess.run('cd '+self.outdir+'# rm *.png *.ps *.txt *.html *.nc namelist*')

        # -------------------------------
        # Cleanup
        for gg in self.outfiles:
            obsfile = 'obs' + gg
            if os.path.exists(obsfile):
                os.remove(obsfile)

        print('Finished: Climate Variability Diagnostics Package ({datetime.datetime.now().strftime("%c")})')


    # ----------------------- task parallelism functions ---------------------------

    def launchTask(self, script):
        # Note -- under python 3.6.x, we have found that for whatever reason, without the
        # stdin-subprocess.DEVNULL, NCL scripts will run to completion and then block
        # indefinitely in the lexer section waiting for input on stdin (???)
        task = subprocess.Popen(self.ncl_exec+ " " + script, shell=True, stdin=subprocess.DEVNULL)

        return task

    def manageTasks(self, taskList):

        # fire off up-to MAX_CONCURRENT subprocesses...
        tasks = list()
        for i, task in enumerate(taskList):
            if i >= self.max_concurrent:
                break
            tasks.append(self.launchTask(task))

        taskList = taskList[len(tasks):]  # remove those scripts we've just launched...

        while len(tasks) > 0:
            finishedList = []
            for task in tasks:
                retCode = task.poll()
                if retCode is not None:
                    #             print "Task status ", task.pid, ": ", task.poll()
                    finishedList.append(task)

                    # more scripts to be run?
                    if len(taskList) > 0:
                        tasks.append(self.launchTask(taskList[0]))
                        del taskList[0]

            for task in finishedList:
                tasks.remove(task)

            time.sleep(self.poll_interval)
        #    print "."      # Feedback to show the script is doing something; not necessary

        print("runTasks.py: Done with CVDP calculation scripts")

if __name__ == '__main__':
    cvdp = CVDPDriver()
    cvdp.run()
