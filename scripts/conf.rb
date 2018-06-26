
#0 PARALLEL, 1 ITERATIVE
PARALLEL_OR_ITERATIVE = 1
PARALLEL_TIME_SLEEP = 60

PATH=ENV["AGILE"] + "/"
PATHMODEL = PATH + "/model/scientific_analysis/data/"
#0 no remote shell, 1 remote shell
#A file token contain
#host, token
#if token = 0 transfer only process
#if token > 1 transfer data and process
REMOTE_SHELL_COMMAND_MULTI = 0

SRCLOCCONFLEVEL = 5.9914659
SQRTTS_CUTTONEXTSTEP = 2.0

#BASEDIR_ARCHIVE = "/scratch1/users/bulgarelli/ARCHIVE/"
#BASEDIR_ARCHIVE = "/ARCHIVE/"
BASEDIR_ARCHIVE = "/AGILE_PROC3/"

TYPE_MATRIX = "I0025"

#ARCHIVE_ID = 0 BUILD17, 1 = BUILD15 or 16
ARCHIVE_ID = 0


load PATH + "scripts/Fits.rb"
load PATH + "scripts/DataUtils.rb"
load PATH + "scripts/AgileFOV.rb"
load PATH + "scripts/MultiOutput.rb"
#load PATH + "scripts/MultiOutput6.rb"
load PATH + "scripts/AlikeUtils.rb"
load PATH + "scripts/DataConversion.rb"
load PATH + "scripts/Parameters.rb"

load "date.rb"
