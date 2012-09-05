#!/usr/bin/env python
""" This is mostly a wrapper around fastq/fasta preprocessing utilities (cutadapt, fastx_collapser), that runs the programs on your input file and keeps track of read-counts and other metadata, saving it in a *_info.txt file.

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

USAGE: deepseq_preprocessing_wrapper.py [options] infile -o outfile_basename """

# basic libraries
from __future__ import division
import sys, os
# other packages
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
# my modules
from general_utilities import write_header_data, run_command_print_info_output, value_and_percentages
from seq_basic_utilities import write_fasta_line, name_seq_generator_from_fasta_fastq
import seq_count_and_lengths


FASTQ_ENCODINGS_FASTX_TOOLKIT = {'auto': '', 'sanger': '-Q33', 'illumina': '-Q64'}


def check_readcount(infile, OUTFILE=None, printing=True, description=None, 
                    total_read_number_only=False, input_collapsed_to_unique=False):
    """ Make sure infile exists and is fa/fq (sys.exit if not); print read count and length distribution 
    to OUTFILE and/or stdout.  Return total number of reads. """
    if not os.path.exists(infile):    
        sys.exit("No %s file found! Exiting program. You should examine and/or remove the intermediate files."%infile)
    read_count_data = seq_count_and_lengths.main([infile], total_read_number_only, input_collapsed_to_unique, 
                                                 include_zeros=False, verbosity=1, OUTPUT=None)[2]
    if description is not None:     header = "# %s file (%s) data:"%(description, infile)
    else:                           header = "# %s file data:"%infile
    output = "%s\n%s"%(header, ''.join(read_count_data))
    if OUTFILE is not None:     OUTFILE.write(output+'\n')
    if printing:                print output
    read_count = int(read_count_data[-1].split(' ')[1])
    return read_count


def _trim_prefix_single(seqname, seq, prefix_bases, TRIMMED_OUTFILE, WRONG_PREFIX_OUTFILE=None):
    """ If prefix_bases is a prefix of seq, trim it off, print to TRIMMED_OUTFILE, return 1; 
    otherwise print to WRONG_PREFIX_OUTFILE if not None, return 0. 
    """
    if seq.startswith(prefix_bases):
        seq_trimmed = seq[len(prefix_bases):]
        write_fasta_line(seqname, seq_trimmed, TRIMMED_OUTFILE)
        return 1
    else:
        if WRONG_PREFIX_OUTFILE is not None:
            write_fasta_line(seqname, seq, WRONG_PREFIX_OUTFILE)
        return 0
    # MAYBE-TODO add an option that specifies the number/percentage of allowed mismatches?


# MAYBE-TODO might be worth making this a separate program
def trim_prefix(prefix_bases, infile, trimmed_outfile, wrong_prefix_outfile=os.devnull, 
               INFOFILE=None, verbosity=1):
    """ Trim prefix_bases from seqs in infile and print to trimmed_outfile; print other seqs to wrong_prefix_outfile.

    Reads fasta or fastq files; outputs fasta files only.
    For each seq in infile, if seq starts with prefix_bases, trim them and print result to trimmed_outfile; 
     otherwise print full seq to wrong_prefix_outfile (if not None). 
    INFOFILE should be an open file handle to print summary info to, or None; 
     verbosity governs how much is printed to stdout
    """
    text = "### Trimming %s from start of each sequence in %s (output to %s, untrimmed to %s)\n"%(prefix_bases, infile, 
                                                                                trimmed_outfile, wrong_prefix_outfile)
    # MAYBE-TODO modify so it can output fastq too?
    if INFOFILE is not None:    INFOFILE.write(text+'\n')
    if verbosity>0:             print text
    N_trimmed, N_untrimmed = 0, 0
    with open(trimmed_outfile, 'w') as TRIMMED_OUTFILE:
        with open(wrong_prefix_outfile, 'w') as WRONG_PREFIX_OUTFILE:
            # MAYBE-TODO right now if wrong_prefix_outfile==None, /dev/null is used - it would be faster with a custom file-like object that doesn't do anything, see http://stackoverflow.com/questions/2929899/cross-platform-dev-null-in-python
            name_seq_generator = name_seq_generator_from_fasta_fastq(infile, verbosity>1)
            for name,seq in name_seq_generator:
                if_trimmed = _trim_prefix_single(name, seq, prefix_bases, TRIMMED_OUTFILE, WRONG_PREFIX_OUTFILE)
                if if_trimmed:  N_trimmed += 1
                else:           N_untrimmed += 1

    N_total = N_trimmed + N_untrimmed
    text = "Trimmed sequences: %s\nUntrimmed sequences: %s\n"%(value_and_percentages(N_trimmed, [N_total]), 
                                                               value_and_percentages(N_untrimmed, [N_total]))
    if INFOFILE is not None:    INFOFILE.write(text+'\n')
    if verbosity>0:             print text
    return N_trimmed, N_untrimmed


def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)
    ### wrapped program options
    parser.add_option('-T','--first_bases_to_trim', metavar='SEQ|NONE', default='ACTA',  
                      help="Sequence to trim from the beginning of each read - reads that don't start with that sequence "
                      +"will end up in another file. Set to NONE to not trim any bases.  Default %default.")
    parser.add_option('-A','--full_cutadapt_options', metavar='"text"', 
                      default='-a GTTGGAACCAAT -e 0.1 -O 10 -n 1 -m 20 -M 21', 
                      help="Full set of options to pass to cutadapt, as a quoted string. Set to NONE to not run it."
                      +' For help on available options (long), run "cutadapt -h" on the command-line. Default "%default".')
    parser.add_option('-C','--collapse_to_unique', action="store_true", default=False, 
                      help="Run output through fastx_collapser to collapse all identical-sequence reads to unique ones. "
                      +"(converts fastq to fasta) (default %default).")
    parser.add_option('-e','--fastq_encoding', type='choice', choices=FASTQ_ENCODINGS_FASTX_TOOLKIT.keys(), default='auto',
                      metavar='|'.join(FASTQ_ENCODINGS_FASTX_TOOLKIT.keys()), 
                      help="fastq quality encoding for fastx_toolkit (default %default).")
    # MAYBE-TODO add options for commonly-used cutadapt options (fastx_collapser doesn't have any); if I do that, make sure that either specific options OR the general option (-A) is provided and used, NOT BOTH, and print the information for the user!
    ### input/output options
    parser.add_option('-n','--total_read_number_only', action="store_true", default=False, 
                      help="Only check the total read number after each step, not read length counts (default %default).")
    parser.add_option('-u','--dont_check_uncollapsed_reads', action="store_true", default=False, 
                      help="Don't check whether the read counts after fastx_collapser (uncollapsed based on header info) "
                      +"match the ones before (by default this check is done and prints a warning if failed, but it "
                      +"takes a while on long files).")
    parser.add_option('-U','--never_check_readcounts_lengths', action="store_true", default=False, 
                      help="Don't check file readcounts/lengths at any stage (to save time) - note that this means "
                      +"no by-stage readcount summary at the end! (default: do check)")
    ### outfile options
    parser.add_option('-o','--outfile_basename', metavar='X', default='test_preprocessing_output', 
                      help="Basename for the two outfiles (which will be X.fa or X.fq, and X_info.txt). Default %default.")
    parser.add_option('-k', '--keep_tmpfiles', action="store_true", default=False, 
                      help="Don't delete the temporary/intermediate files (all will have the same prefix as the outfile, "
                      +"plus _tmp*; for exact names see full commands printed to stdout and info file).")
    # MAYBE-TODO make option to generate and keep both the collapsed-to-unique file and the uncollapsed file?  That would be useful - may want it as a default.  Give them sensible names, like _full and _unique.
    # MAYBE-TODO the seq_count_and_lengths step etc could be optional, to speed things up - make option for that?  To do things quickly without intermediate checks, use a single command with pipes, instead of multiple commands with intermediate files?  But that would make disentangling the stuff printed by each command complicated...
    ### printing options
    parser.set_defaults(verbosity=1)
    parser.add_option('-q', '--quiet', action="store_const", const=1, dest='verbosity',
                      help="Only print commands and the output summary to STDOUT.")
    parser.add_option('-Q', '--very_quiet', action="store_const", const=0, dest='verbosity',
                      help="Don't print anything to STDOUT.")
    parser.add_option('-v', '--verbose', action="store_const", const=2, dest='verbosity',
                      help="Print read-length counts at each step, and each program output, to STDOUT.")
    return parser


if __name__ == "__main__":
    """ Allows both running and importing of this file. """

    parser = define_option_parser()
    try:
        (options, [infile]) = parser.parse_args()
    except ValueError:
        parser.print_help()
        sys.exit("Error: exactly one infile required!")
    # MAYBE-TODO implement option with multiple infiles? Need to make sure they're the same fa/fq type etc...

    ### outfile and tmpfile names
    infile_suffix = os.path.splitext(infile)[1]
    outfile_suffix = '.fa'
    #outfile_suffix = '.fa' if options.collapse_to_unique else infile_suffix
    outfile = options.outfile_basename + outfile_suffix
    infofile = options.outfile_basename + '_info.txt'
    bad_prefix_file = options.outfile_basename + '_bad_prefix.fa'
    no_cassette_file = options.outfile_basename + '_no_cassette.fa'
    tmpfile1 = tmpfile1_original = options.outfile_basename + '_tmpfile1.fa'
    tmpfile2 = tmpfile2_original = options.outfile_basename + '_tmpfile2.fa'
    
    with open(infofile,'w') as INFOFILE:

        ### write header data
        write_header_data(INFOFILE,options)
        INFOFILE.write('\n')

        ### 1. look at the infile; make sure it's readable, etc
        #       (check_readcount uses seq_count_and_lengths, which uses HTSeq and autodetects fa/fq format)
        if not options.never_check_readcounts_lengths:
            starting_readcount = check_readcount(infile, INFOFILE, bool(options.verbosity>1), "original input", 
                                                 options.total_read_number_only, False)

        ### 2. Trim the first bases (from adapter)
        if options.first_bases_to_trim == 'NONE':
            if options.verbosity>0:   print "### Not trimming first bases, since NONE was passed to -T option.\n"
            tmpfile1 = infile
            trimmed_readcount = starting_readcount
        else:
            trim_prefix(options.first_bases_to_trim, infile, tmpfile1, bad_prefix_file, INFOFILE, options.verbosity)
            if not options.never_check_readcounts_lengths:
                trimmed_readcount = check_readcount(tmpfile1, INFOFILE, bool(options.verbosity>1), 
                                                    "first-base-trimming output", options.total_read_number_only, False)
            else: INFOFILE.write('\n')

        ### 3. run cutadapt to strip cassette sequence
        if options.full_cutadapt_options == 'NONE':
            if options.verbosity>0:   print "### Not running cutadapt, since NONE was passed to -A option.\n"
            tmpfile2 = tmpfile1
            cutadapt_readcount = trimmed_readcount
        else:
            for extra_seq_category in ('untrimmed', 'too-short', 'too-long'):
                if not extra_seq_category in options.full_cutadapt_options:
                    options.full_cutadapt_options += ' --%s-output %s'%(extra_seq_category, no_cassette_file)
            # NOTE: this currently requires my version of cutadapt, cutadapt_mod, to deal with too-long seqs correctly
            command = "cutadapt_mod %s -o %s %s"%(options.full_cutadapt_options, tmpfile2, tmpfile1)
            run_command_print_info_output(command, INFOFILE, bool(options.verbosity>0), shell=True)
            if not options.never_check_readcounts_lengths:
                cutadapt_readcount = check_readcount(tmpfile2, INFOFILE, bool(options.verbosity>1), "cutadapt output", 
                                                    options.total_read_number_only, False)

        ### 4. run fastx_collapser
        if not options.collapse_to_unique:
            if options.verbosity>0:   print "### Not running fastx_collapser, since -C option was not used.\n"
            os.rename(tmpfile2,outfile)
            collapsed_readcount = cutadapt_readcount
            # Note for fastx_collapser, but also for the others - NONE is necessary here, can't just use '', because 
            #    fastx_collapser works fine with no options, so '' is a sensible input and can't be used to turn it off.
        else:
            command = "fastx_collapser -v %s -i %s -o %s"%(FASTQ_ENCODINGS_FASTX_TOOLKIT[options.fastq_encoding], 
                                                           tmpfile2, outfile)
            run_command_print_info_output(command, INFOFILE, bool(options.verbosity>0), shell=True)
            INFOFILE.write('\n')
            if not options.never_check_readcounts_lengths:
                collapsed_readcount = check_readcount(outfile,INFOFILE,bool(options.verbosity>1),"fastx_collapser output", 
                                                      options.total_read_number_only, input_collapsed_to_unique=False)
            # make sure uncollapsed readcount is the same as before collapsing
            if not (options.never_check_readcounts_lengths or options.dont_check_uncollapsed_reads):
                uncollapsed_readcount = check_readcount(outfile, None, False, "", True, input_collapsed_to_unique=True)
                if not uncollapsed_readcount == cutadapt_readcount:
                    text = "ERROR: the uncollapsed read-count after fastx_collapser isn't the same as the before-collapser count!  Collapsing went wrong somehow, or the way fastx_collapser works changed since this program was written?\n"
                else:
                    text = "(checked that all the reads are still there if you uncollapse the numbers using header info)\n"
                print text
                INFOFILE.write(text+'\n')
            # also run fastx_collapser on bad_prefix_file and no_cassette_file
            if options.verbosity:
                print "### Running fastx_collapser on the additional output files. Not printing the output to info file."
            for extra_file in (bad_prefix_file, no_cassette_file):
                command = "fastx_collapser -v %s -i %s -o tmp.fa"%(FASTQ_ENCODINGS_FASTX_TOOLKIT[options.fastq_encoding], 
                                                                   extra_file)
                retcode = run_command_print_info_output(command, None, bool(options.verbosity>1), shell=True)
                if retcode in (0, None):
                    os.remove(extra_file)
                    os.rename('tmp.fa', extra_file)

        if not options.never_check_readcounts_lengths:
            final_output = ["### Final read count info\n"]
            final_output.append("# starting read count:   %s\n"%starting_readcount)
            if not options.first_bases_to_trim == 'NONE':
                final_output.append("# read count after trimming:   %s\n"%trimmed_readcount)
            if not options.full_cutadapt_options == 'NONE':
                final_output.append("# read count after cassette stripping:   %s\n"%cutadapt_readcount)
            discarded_readcount = starting_readcount-cutadapt_readcount
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
        for tmpfile in [tmpfile1_original, tmpfile2_original]:
            if os.path.exists(tmpfile):     os.remove(tmpfile)

# TODO write unit-tests or run-tests for all this

