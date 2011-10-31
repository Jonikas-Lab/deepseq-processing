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

USAGE: deepseq_preprocessing.py [options] infile -o outfile_basename """

# basic libraries
import sys, os, subprocess
# my modules
from general_utilities import write_header_data
from seq_count_and_lengths import seq_count_and_lengths


def check_readcount(infile, OUTFILE=None, printing=True, description=None, 
                    total_read_number_only=False, input_collapsed_to_unique=False):
    """ Make sure infile exists and is fa/fq (sys.exit if not); print read count and length distribution 
    to OUTFILE and/or stdout.  Return total number of reads. """
    if not os.path.exists(infile):    
        sys.exit("No %s file found! Exiting program. You should examine and/or remove the intermediate files."%infile)
    read_count_data = seq_count_and_lengths([infile], total_read_number_only, input_collapsed_to_unique, 
                                             include_zeros=False, verbosity=1, OUTPUT=None)
    if description is not None:     header = "# %s file (%s) data:"%(description, infile)
    else:                           header = "# %s file data:"%infile
    output = "%s:\n%s"%(header, ''.join(read_count_data))
    if OUTFILE is not None:     OUTFILE.write(output+'\n')
    if printing:                print output
    read_count = int(read_count_data[-1].split(' ')[1])
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


# MAYBE-TODO write unit-tests?


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    from optparse import OptionParser
    parser = OptionParser(__doc__)
    ### wrapped program options
    parser.add_option('-T','--full_fastx_trimmer_options', metavar='"text"', default='-f 5',  
                      help="Full set of options to pass to fastx_trimmer, as a quoted string. Set to NONE to not run it."
                      +' For help on available options, run "fastx_trimmer -h" on the command-line. Default "%default".')
    parser.add_option('-A','--full_cutadapt_options', metavar='"text"', 
                      default='-a GTTGGAACCAAT -e 0.1 -O 10 -n 1 -m 20 -M 21', 
                      help="Full set of options to pass to cutadapt, as a quoted string. Set to NONE to not run it."
                      +' For help on available options (long), run "cutadapt -h" on the command-line. Default "%default".')
    parser.add_option('-C','--collapse_to_unique', action="store_true", default=False, 
                      help="Run output through fastx_collapser to collapse all identical-sequence reads to unique ones. "
                      +"(converts fastq to fasta) (default %default).")
    # MAYBE-TODO add options for commonly-used fastx_trimmer and cutadapt options (fastx_collapser doesn't have any); if I do that, make sure that for each wrapped program, either specific options OR the general option (-T/-A) is provided and used, NOT BOTH, and print the information for the user!
    ### input/output options
    parser.add_option('-n','--total_read_number_only', action="store_true", default=False, 
                      help="Only check the total read number after each step, not read length counts (default %default).")
    parser.add_option('-u','--dont_check_uncollapsed_reads', action="store_true", default=False, 
                      help="Don't check whether the read counts after fastx_collapser (uncollapsed based on header info) match the ones before (by default this check is done and prints a warning if failed, but it takes a while on long files).")
    ### outfile options
    parser.add_option('-o','--outfile_basename', metavar='X', default='test', 
                      help="Basename for the two outfiles (X.fa or X.fq, and X_info.txt). Default %default.")
    parser.add_option('-k', '--keep_tmpfiles', action="store_true", default=False, 
                      help="Don't delete the temporary/intermediate files (all will have the same prefix as the outfile, "
                      +"plus _tmp*; for exact names see full commands printed to stdout and info file).")
    # MAYBE-TODO make option to generate and keep both the collapsed-to-unique file and the uncollapsed file?  That would be useful - may want it as a default.  Give them sensible names, like _full and _unique.
    # MAYBE-TODO the seq_count_and_lengths step etc could be optional, to speed things up - make option for that?  To do things quickly, use a single command with pipes, instead of multiple commands with intermediate files.
    ### printing options
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
        #       (check_readcount uses seq_count_and_lengths, which uses HTSeq and autodetects fa/fq format)
        starting_readcount = check_readcount(infile, INFOFILE, bool(options.verbosity>1), "original input", 
                                             options.total_read_number_only, False)

        ### 2. run fastx_trimmer
        if options.full_fastx_trimmer_options == 'NONE':
            if options.verbosity>0:   print "# Not running fastx_trimmer, since NONE was passed to -T option."
            tmpfile1 = infile
            trimmed_readcount = starting_readcount
        else:
            command = "fastx_trimmer %s -i %s -o %s"%(options.full_fastx_trimmer_options, infile, tmpfile1)
            run_command_and_print_info(command, INFOFILE, bool(options.verbosity>0), shell=True)
            trimmed_readcount = check_readcount(tmpfile1, INFOFILE, bool(options.verbosity>1), "fastx_trimmer output", 
                                                options.total_read_number_only, False)

        ### 3. run cutadapt
        if options.full_cutadapt_options == 'NONE':
            if options.verbosity>0:   print "# Not running cutadapt, since NONE was passed to -A option."
            tmpfile2 = tmpfile1
            adapter_readcount = trimmed_readcount
        else:
            if not 'untrimmed' in options.full_cutadapt_options:
                options.full_cutadapt_options += ' --untrimmed /dev/null'
            command = "cutadapt %s -o %s %s > %s"%(options.full_cutadapt_options, tmpfile2, tmpfile1, tmp_cutadapt_info)
            run_command_and_print_info(command, INFOFILE, bool(options.verbosity>0), shell=True)
            print_text_from_file(tmp_cutadapt_info, INFOFILE, bool(options.verbosity>1))
            adapter_readcount = check_readcount(tmpfile2, INFOFILE, bool(options.verbosity>1), "cutadapt output", 
                                                options.total_read_number_only, False)

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
            collapsed_readcount = check_readcount(outfile, INFOFILE, bool(options.verbosity>1), "fastx_collapser output", 
                                                  options.total_read_number_only, input_collapsed_to_unique=False)
            # make sure uncollapsed readcount is the same as before collapsing
            if not options.dont_check_uncollapsed_reads:
                uncollapsed_readcount = check_readcount(outfile, None, False, "", True, input_collapsed_to_unique=True)
                if not uncollapsed_readcount == adapter_readcount:
                    text = "ERROR: the uncollapsed read-count after fastx_collapser isn't the same as the before-collapser count!  Collapsing went wrong somehow, or the way fastx_collapser works changed since this program was written?\n"
                else:
                    text = "(checked that all the reads are still there if you uncollapse the numbers using header info)\n"
                print text
                INFOFILE.write(text+'\n')

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

    ### Remove tmpfiles
    # need to use the tmpfile*_original names here because I do "tmpfile1 = infile" etc if skipping steps, 
    #   and I don't want to remove the infile!
    if not options.keep_tmpfiles:
        for tmpfile in [tmpfile1_original, tmpfile2_original, tmp_cutadapt_info, tmp_collapser_info]:
            if os.path.exists(tmpfile):     os.remove(tmpfile)

