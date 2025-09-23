The module structure is always the same:

**setup.sh**: this is executed each time the 
application starts, it should download the 
necessary data for the module and annotate
it. Any pre processing or set up of sort
should be handle there.
It takes no args.

**enrich.sh**: this is one is called each time
the core pipeline runs. 
args:
USER_DIR= the job specifc dir
OUPUT_DIR= output dir
BACKGROUND_NAME= which bg to use
P_VALUE= the p va;
CORRECTION= correctino method

**multisample.sh**: runs the final part of 
the pipeline 3. 
args:
USER_DIR= same as enrich.sh
METHOD= this is not used for now, but in the
future you might give the user an option on
what method to use (fisher / ttest /whatever)

all other files are module specific.

Note for the future: a good idea would be to
move the graph logic inside each module, for 
now we support only 3 plot types, which one
is going to be used depends on the file name.
check App/Utils/view.py for more. 