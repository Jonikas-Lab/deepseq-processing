#!/usr/bin/env python
""" Preprocessing of insertional-mutant deepseq reads: strip first bases (from adapter) and cassette sequence, filter by length, collapse to unique (all stages are optional and customizable). 

1) Strip a given start sequence from each read; reads that didn't start with that sequence go in a separate file.
Uses my own trim_prefix function (used to be done with fastx_trimmer, but that didn't check what the sequence was)

2) Strip the cassette sequence from middle/end of read, and grab only the remaining reads that are within the expected length range; reads that don't contain the cassette or are too short/long go in a separate file. 
This is a wrapper around the cutadapt command-line program (or my modified version of it).

3) Collapse identical read sequences to single reads (encoding the original readcount in the header) to avoid aligning the same sequence multiple times. This discards the quality info (i.e. converts file from fastq to fasta), if it hasn't been done earlier. 
This is a wrapper around the fastx_collapser command-line program that does exactly that.

My program keeps track of read-counts and metadata and saves it all in a metadata file.

Output - there are actually four outfiles generated (with names based on the -o option - calling the value X here):
    - sequence file - <X>.fq or .fa, depending on input format etc - contains the final sequences
        after trimming the start sequence, the cassette sequence (cutadapt) and collapsing to unique (fastx_collapser)
    - metadata file - <X>_info.txt - contains general information: 
        exact command used, date/time, system, user, etc
        read-length distribution before/after each step (may be optional), 
        summary output from the wrapped programs (normally printed to stdout), 
        sequence counts from each stage, with the total count/percentage of discarded sequences
    - <X>_wrong-start.fa - sequences that didn't start with the start sequence given in the -F option;
        collapsed to unique if -C option.
    - <X>_no-cassette.fa - sequences that didn't contain the cassette sequence given in the -A option;
        collapsed to unique if -C option.
    - temporary/intermediate files - <X>_tmp* - contain intermediate sequence data (after start-trimming but
         before cassette-stripping, and after cassette-stripping but before collapsing to unique if -C option given).
        Removed after program finishes running with no errors, unless option -k is used.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011-2012

USAGE: deepseq_preprocessing_wrapper.py [options] infile -o outfile_basename """

# standard library
from __future__ import division
import sys, os
import unittest
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

    ### basic functionality options
    parser.add_option('-F','--first_bases_to_trim', metavar='SEQ|NONE', default='ACTA',  
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

    ### extra readcount/length check options
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

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    return parser


def main(args, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """
    try:
        [infile] = args
    except ValueError:
        parser = define_option_parser()
        parser.print_help()
        sys.exit("Error: exactly one infile required!")
    # MAYBE-TODO implement option with multiple infiles? Need to make sure they're the same fa/fq type etc...

    ### outfile and tmpfile names
    infile_suffix = os.path.splitext(infile)[1]
    outfile_suffix = '.fa'
    #outfile_suffix = '.fa' if options.collapse_to_unique else infile_suffix
    outfile = options.outfile_basename + outfile_suffix
    infofile = options.outfile_basename + '_info.txt'
    wrong_start_file = options.outfile_basename + '_wrong-start.fa'
    no_cassette_file = options.outfile_basename + '_no-cassette.fa'
    tmpfile1 = tmpfile1_original = options.outfile_basename + '_tmpfile1.fa'
    tmpfile2 = tmpfile2_original = options.outfile_basename + '_tmpfile2.fa'
    
    with open(infofile,'w') as INFOFILE:

        ### write header data
        write_header_data(INFOFILE,options)
        INFOFILE.write('\n')

        ### 0. look at the infile; make sure it's readable, etc
        #       (check_readcount uses seq_count_and_lengths, which uses HTSeq and autodetects fa/fq format)
        if not options.never_check_readcounts_lengths:
            starting_readcount = check_readcount(infile, INFOFILE, bool(options.verbosity>1), "original input", 
                                                 options.total_read_number_only, False)

        ### 1. Trim the first bases (from adapter)
        # MAYBE-TODO I could do this with cutadapt again, instead of with my own trim_prefix function... 
        #  Would that be faster, or better in any other way?
        # MAYBE-TODO could also do it with a multiplexing barcode-splitting tool (like fastx_barcode_splitter.pl), 
        #  since that's the eventual point of having those constant first bases there...
        if options.first_bases_to_trim == 'NONE':
            if options.verbosity>0:   print "### Not trimming first bases, since NONE was passed to -F option.\n"
            tmpfile1 = infile
            trimmed_readcount = starting_readcount
            untrimmed_readcount = 0
        else:
            trim_prefix(options.first_bases_to_trim, infile, tmpfile1, wrong_start_file, INFOFILE, options.verbosity)
            if not options.never_check_readcounts_lengths:
                trimmed_readcount = check_readcount(tmpfile1, INFOFILE, bool(options.verbosity>1), 
                                                    "first-base-trimming output", options.total_read_number_only, False)
                untrimmed_readcount = check_readcount(wrong_start_file, None, False, True, False)
                assert trimmed_readcount+untrimmed_readcount==starting_readcount,\
                        "Trimmed/untrimmed readcounts don't add up to starting readcount - check tmpfile!"\
                        +"(%s+%s != %s)"%(trimmed_readcount, untrimmed_readcount, starting_readcount)
            else: INFOFILE.write('\n')

        ### 2. run cutadapt to strip cassette sequence
        if options.full_cutadapt_options == 'NONE':
            if options.verbosity>0:   print "### Not running cutadapt, since NONE was passed to -A option.\n"
            tmpfile2 = tmpfile1
            cutadapt_readcount = trimmed_readcount
            no_cassette_readcount = 0
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
                no_cassette_readcount = check_readcount(no_cassette_file, None, False, True, False)
                assert cutadapt_readcount+no_cassette_readcount==trimmed_readcount,\
                        "Cassette/no-cassette readcounts don't add up to trimmed readcount - check tmpfile!"\
                        +"(%s+%s != %s)"%(cutadapt_readcount, no_cassette_readcount, trimmed_readcount)

        ### 3. run fastx_collapser to collapse the sequences to unique
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
                if options.verbosity: print text
                INFOFILE.write(text+'\n')
            # also run fastx_collapser on wrong_start_file and no_cassette_file
            if options.verbosity:
                print "### Running fastx_collapser on the additional output files. Not printing the output to info file."
            extra_collapsed_readcounts = {}    
            for extra_file in (wrong_start_file, no_cassette_file):
                command = "fastx_collapser -v %s -i %s -o tmp.fa"%(FASTQ_ENCODINGS_FASTX_TOOLKIT[options.fastq_encoding], 
                                                                   extra_file)
                retcode = run_command_print_info_output(command, None, bool(options.verbosity>1), shell=True)
                if retcode in (0, None):
                    os.remove(extra_file)
                    os.rename('tmp.fa', extra_file)
                extra_collapsed_readcounts[extra_file] = check_readcount(extra_file, None, False, "", True, 
                                                                         input_collapsed_to_unique=False)

        ### Final readcount check
        if not options.never_check_readcounts_lengths:
            final_output = ["### Final read count info for %s (main output file %s)\n"%(infile, outfile)]
            final_output.append("# starting total read count:\t%s\n"%starting_readcount)
            if not options.first_bases_to_trim == 'NONE':
                final_output.append('# "good" read count after start trimming (%% of total):\t%s\n'%
                                    value_and_percentages(trimmed_readcount, [starting_readcount]))
                final_output.append('#  "bad" read count (wrong-start) (%% of total):\t%s\n'%
                                    value_and_percentages(untrimmed_readcount, [starting_readcount]))
            if not options.full_cutadapt_options == 'NONE':
                final_output.append('# "good" read count after cassette stripping (%% of total, %% of trimmed):\t%s\n'%
                                    value_and_percentages(cutadapt_readcount, [starting_readcount, trimmed_readcount]))
                final_output.append('#  "bad" read count (no-cassette) (%% of total, %% of trimmed):\t%s\n'%
                                    value_and_percentages(no_cassette_readcount, [starting_readcount, trimmed_readcount]))
            final_output.append('## final "good" reads (in main output file) (%% of total):\t%s\n'%
                                value_and_percentages(cutadapt_readcount, [starting_readcount]))
            final_output.append('## final "bad" reads (in _wrong-start and/or _no-cassette files) (%% of total):\t%s\n'%
                                value_and_percentages(starting_readcount-cutadapt_readcount, [starting_readcount]))
            if options.collapse_to_unique:
                final_output.append('# "good" unique sequence count after collapsing reads to unique sequences '
                                    '(%% of read count):\t%s\n'
                                    %value_and_percentages(collapsed_readcount, [cutadapt_readcount]))
                if not options.first_bases_to_trim == 'NONE':
                    final_output.append('# wrong-start unique sequence count after collapsing (%% of read count):\t%s\n'
                            %value_and_percentages(extra_collapsed_readcounts[wrong_start_file], [untrimmed_readcount]))
                if not options.full_cutadapt_options == 'NONE':
                    final_output.append('# no-cassette unique sequence count after collapsing (%% of read count):\t%s\n'
                            %value_and_percentages(extra_collapsed_readcounts[no_cassette_file], [no_cassette_readcount]))
            for line in final_output:
                INFOFILE.write(line)
                if options.verbosity>0:  print line,

    ### Remove tmpfiles
    # need to use the tmpfile*_original names here because I do "tmpfile1 = infile" etc if skipping steps, 
    #   and I don't want to remove the infile!
    if not options.keep_tmpfiles:
        for tmpfile in [tmpfile1_original, tmpfile2_original]:
            if os.path.exists(tmpfile):     os.remove(tmpfile)


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    from testing_utilities import run_functional_tests
    test_folder = 'test_data'
    infile1 = test_folder + '/INPUT_raw_fastq.fq'
    # tests in (testname, test_description, arg_string, infiles) format
    test_runs = [ 
        ('not-unique', "full functionality without collapsing to unique", "", [infile1]),
        ('unique', "full functionality with collapsing to unique (for all output files)", "-C", [infile1]), 
        ]
    # convert tests into (testname, arg_and_infile_string) format, adding the options that are always used
    test_runs = [('pre__'+testname, descr, test_args+' -Q '+' '.join(infiles)) 
                  for testname,descr,test_args,infiles in test_runs]
    # argument_converter converts (parser,options,args) to the correct argument order for main
    argument_converter = lambda parser,options,args: (args, options)
    # use my custom function to run all the tests, auto-detect reference files, compare to output.
    return run_functional_tests(test_runs, define_option_parser(), main, test_folder, 
                                argument_converter=argument_converter, outfile_option='-o', 
                                append_to_outfilenames='') 


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # MAYBE-TODO add unit-tests


if __name__=='__main__':
    parser = define_option_parser()
    options,args = parser.parse_args()

    # if run with -t option, do unit tests and quit
    if options.test_functionality:
        print("*** You used the -t option - ignoring all other options/arguments, running the built-in test suite. ***")
        # to run tests for another file, have to use TextTestRunner, not unittest.main -  make a test suite with 
        #   autodetection of all tests (see http://docs.python.org/library/unittest.html#unittest.TestLoader)
        #print("\n * unit-tests for the ______ module")
        #test_suite_1 = unittest.defaultTestLoader.loadTestsFromModule(______
        #unittest.TextTestRunner(verbosity=1).run(test_suite_1)
        # to run tests for current module, just run unittest.main, passing it only the filename 
        #   (by default it takes all of sys.argv and complains about options/arguments it can't recognize)
        print("\n * unit-tests for this module (%s)"%sys.argv[0])
        unittest.main(argv=[sys.argv[0]])   # unittest.main automatically runs sys.exit()

    if options.test_run:
        print("*** You used the -T option - ignoring all other options and running the built-in example test runs. ***")
        test_result = do_test_run()
        sys.exit(test_result)

    # otherwise pass the arguments to the main function
    main(args, options)

