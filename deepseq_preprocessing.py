#!/usr/bin/env python
""" This is mostly a wrapper around fastq/fasta preprocessing utilities (fastx_trimmer, cutadapt, fastx_collapser), that runs the programs on your input file and keeps track of read-counts and other metadata, saving it in a *_info.txt file.

Use the -T, -A, -C options to provide full option lists to each wrapped program - the option lists will be passed to the program with no parsing or changes except for specifying input/output files and possibly verbosity levels (the full command-lines, as ran, are printed to stdout).  I may add separate options for commonly used options of the wrapped programs in the future.

Output - there are actually two outfiles generated (with names based on the -o option - calling the value X here):
    - sequence file - <X>.fq or .fa, depending on input format etc - contains the sequences
    - metadata file - <X>_info.txt - contains general information: 
        exact command used, date/time, system, user, etc
        read-length distribution before/after each step (may be optional), 
        summary output from the wrapped programs (normally printed to stdout), 
        sequence counts from each stage, with the total count/percentage of discarded sequences
    - temporary/intermediate files - <X>_tmp* - contain intermediate sequence data, or output from wrapped programs.
        Removed after program finishes running with no errors, unless option -k is used.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, Oct 2011

USAGE: deepseq_preprocessing.py [options] infile """

import sys, os, subprocess
from general_utilities import write_header_data
from read_length_distribution import read_length_distribution


def check_readlengths(infile, OUTFILE=None, printing=True, description=None):
    """ Make sure infile exists and is fa/fq (sys.exit if not); print read length distribution to OUTFILE and/or stdout.
    Return total number of reads. """
    if not os.path.exists(infile):    
        sys.exit("No %s file found! Exiting program. You should examine and remove the intermediate files."%description)
        # MAYBE-TODO should the intermediate files still be removed if exiting?
    read_length_data = read_length_distribution([infile],include_zeros=False,verbose=1,OUTPUT=None)
    if description is not None:     header = "# %s file (%s) data:"%(description, infile)
    else:                           header = "# %s file data:"%infile
    output = "%s:\n%s"%(header, ''.join(read_length_data))
    if OUTFILE is not None:     OUTFILE.write(output+'\n')
    if printing:                print output
    read_count = int(read_length_data[-1].split(' ')[1])
    return read_count


def run_command_and_print_info(command, LOGFILE=None, printing=True, shell=True, program_name=None):
    """ Run command using subprocess.call; first print a line describing that to LOGFILE and/or stdout.
    The shell arg to subprocess.call is given by shell; LOGFILE should be an open file object; 
    program_name is only used for printing, and the first word of the command will be used by default. """
    if program_name is None:
        program_name = command.split(' ')[0]
    output = "### Running %s: %s"%(program_name, command)
    if LOGFILE is not None:     LOGFILE.write(output+'\n')
    if printing:                print output
    subprocess.call([command], shell=shell)


def print_text_from_file(infile, OUTFILE=None, printing=True, add_newlines=0):
    """ Write all text from infile to OUTFILE (if not None), also print to stdout if printing is set. 
    Return line counts.  Infile should be a filename; OUTFILE should be an open file object. """
    line_count = 0
    for line in open(infile):
        if OUTFILE is not None:     OUTFILE.write(line)
        if printing:                print line,
        line_count += 1
    if add_newlines:
        if OUTFILE is not None:     OUTFILE.write('\n'*add_newlines)
        if printing:                print '\n'*add_newlines,
    return line_count


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    from optparse import OptionParser
    parser = OptionParser(__doc__)
    parser.add_option('-T','--full_fastx_trimmer_options', metavar='"text"', default='-f 5',  
                      help="Full set of options to pass to fastx_trimmer, as a quoted string. Set to NONE to not run it."
                      +' For help on available options, run "fastx_trimmer -h" on the command-line. Default "%default".')
    parser.add_option('-A','--full_cutadapt_options', metavar='"text"', default='-a GTTGGAACCAAT -e 0.1 -n 1 -m 20 -M 21', 
                      help="Full set of options to pass to cutadapt, as a quoted string. Set to NONE to not run it."
                      +' For help on available options (long), run "cutadapt -h" on the command-line. Default "%default".')
    # TODO should I be using the -O LENGTH, --overlap=LENGTH cutadapt option?  Test!  Put in defaults if so.
    parser.add_option('-C','--collapse_to_unique', action="store_true", default=False, 
                      help="Run output through fastx_collapser to collapse all identical-sequence reads to unique ones. "
                      +"(converts fastq to fasta) (default %default).")
    parser.add_option('-o','--outfile_basename', metavar='X', default='test', 
                      help="Basename for the two outfiles (X.fa or X.fq, and X_info.txt). Default %default.")
    # MAYBE-TODO add options for commonly-used fastx_trimmer and cutadapt options (fastx_collapser doesn't have any); if I do that, make sure that for each wrapped program, either specific options OR the general option (-T/-A) is provided and used, NOT BOTH, and print the information for the user!
    parser.add_option('-k', '--keep_tmpfiles', action="store_true", default=False, 
                      help="Don't delete the temporary/intermediate files (all will have the same prefix as the outfile).")
    parser.set_defaults(verbosity=1)
    parser.add_option('-q', '--quiet', action="store_const", const=1, dest='verbosity',
                      help="Only print commands and the output summary to STDOUT.")
    parser.add_option('-Q', '--very_quiet', action="store_const", const=0, dest='verbosity',
                      help="Don't print anything to STDOUT.")
    parser.add_option('-v', '--verbose', action="store_const", const=2, dest='verbosity',
                      help="Print read-length counts at each step, and each program output, to STDOUT.")
    try:
        (options, [infile]) = parser.parse_args()
    except ValueError:
        parser.print_help()
        sys.exit("Error: exactly one infile required!")
    # MAYBE-TODO implement option with multiple infiles? Need to make sure they're the same fa/fq type etc...

    ### outfile and tmpfile names
    infile_suffix = os.path.splitext(infile)[1]
    outfile_suffix = '.fa' if options.collapse_to_unique else infile_suffix
    outfile = options.outfile_basename + outfile_suffix
    infofile = options.outfile_basename + '_info.txt'
    tmpfile1 = tmpfile1_original = options.outfile_basename + '_tmpfile1' + infile_suffix
    tmpfile2 = tmpfile2_original = options.outfile_basename + '_tmpfile2' + infile_suffix
    tmp_cutadapt_info = options.outfile_basename + '_tmp_cutadapt_info.txt'
    tmp_collapser_info = options.outfile_basename + '_tmp_collapser_info.txt'
    
    with open(infofile,'w') as INFOFILE:

        ### write header data
        write_header_data(INFOFILE,options)
        INFOFILE.write('\n')

        ### 1. look at the infile; make sure it's readable, etc
        #       (read_length_distribution uses HTSeq and autodetects fa/fq format)
        starting_readcount = check_readlengths(infile, INFOFILE, bool(options.verbosity>1), "original input")

        ### 2. run fastx_trimmer
        if options.full_fastx_trimmer_options == 'NONE':
            if options.verbosity>0:   print "# Not running fastx_trimmer, since NONE was passed to -T."
            tmpfile1 = infile
            trimmed_readcount = starting_readcount
        else:
            command = "fastx_trimmer %s -i %s -o %s"%(options.full_fastx_trimmer_options, infile, tmpfile1)
            run_command_and_print_info(command, INFOFILE, bool(options.verbosity>0), shell=True)
            trimmed_readcount = check_readlengths(tmpfile1, INFOFILE, bool(options.verbosity>1), "fastx_trimmer output")

        ### 3. run cutadapt
        if options.full_cutadapt_options == 'NONE':
            if options.verbosity>0:   print "# Not running cutadapt, since NONE was passed to -A."
            tmpfile2 = tmpfile1
            adapter_readcount = trimmed_readcount
        else:
            if not 'untrimmed' in options.full_cutadapt_options:
                options.full_cutadapt_options += ' --untrimmed /dev/null'
            command = "cutadapt %s -o %s %s > %s"%(options.full_cutadapt_options, tmpfile2, tmpfile1, tmp_cutadapt_info)
            run_command_and_print_info(command, INFOFILE, bool(options.verbosity>0), shell=True)
            print_text_from_file(tmp_cutadapt_info, INFOFILE, bool(options.verbosity>1))
            adapter_readcount = check_readlengths(tmpfile2, INFOFILE, bool(options.verbosity>1), "cutadapt output")

        ### 4. run fastx_collapser
        if not options.collapse_to_unique:
            if options.verbosity>0:   print "# Not running fastx_collapser, since -C option was not used."
            os.rename(tmpfile2,outfile)
            collapsed_readcount = adapter_readcount
            # Note for fastx_collapser, but also for the others - NONE is necessary here, can't just use '', because 
            #    fastx_collapser works fine with no options, so '' is a sensible input and can't be used to turn it off.
        else:
            command = "fastx_collapser -v -i %s -o %s > %s"%(tmpfile2, outfile, tmp_collapser_info)
            run_command_and_print_info(command, INFOFILE, bool(options.verbosity>0), shell=True)
            print_text_from_file(tmp_collapser_info, INFOFILE, bool(options.verbosity>1), add_newlines=1)
            collapsed_readcount = check_readlengths(outfile, INFOFILE, bool(options.verbosity>1), "fastx_collapser output")

        final_output = ["### Final read count info\n"]
        final_output.append("# starting read count:   %s\n"%starting_readcount)
        if not options.full_fastx_trimmer_options == 'NONE':
            final_output.append("# read count after trimming:   %s\n"%trimmed_readcount)
        if not options.full_cutadapt_options == 'NONE':
            final_output.append("# read count after adapter stripping:   %s\n"%adapter_readcount)
        discarded_readcount = starting_readcount-adapter_readcount
        discarded_percent = 100*float(discarded_readcount)/starting_readcount
        final_output.append("## reads removed:   %s  (%.0f%% of total)\n"%(discarded_readcount, discarded_percent))
        if options.collapse_to_unique:
            final_output.append("# read count after collapsing to unique sequences:   %s\n"%collapsed_readcount)
        for line in final_output:
            INFOFILE.write(line)
            if options.verbosity>0:
                print line,

        # MAYBE-TODO the read_length_distribution step etc should probably be optional, 
        #   since it's going to take extra time/cpu and the files may be large...  Make an option for that, implement.
        # To do things quickly, just use a single command with pipes, instead of multiple commands with intermediate files.
        # MAYBE-TODO also maybe have an intermediate level where the output doesn't include 
        #   the whole read-length distribution but just the total read count?

    ### Remove tmpfiles
    # need to use the tmpfile*_original names here because I do "tmpfile1 = infile" etc if skipping steps, 
    #   and I don't want to remove the infile!
    if not options.keep_tmpfiles:
        for tmpfile in [tmpfile1_original, tmpfile2_original, tmp_cutadapt_info, tmp_collapser_info]:
            if os.path.exists(tmpfile):     os.remove(tmpfile)

