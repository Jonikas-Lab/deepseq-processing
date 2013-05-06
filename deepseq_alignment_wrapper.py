#!/usr/bin/env python
""" Align sequences against index (using bowtie), categorize, write to one/multiple outfiles, write metadata to info file.

Use the command-line options to specify bowtie options. For bowtie help, run "bowtie -h" on the command-line.  This program uses -v mode only.  
There is an option to specify two different indexes - one genome and one cassette.  The purpose of that is to always categorize genome-and-cassette multiple-alignments as cassette, without needing to see all of the genome alignments (since there can be thousands).  The processing is as follows:
 1) align the infile against the genome and cassette
 2) merge the two files together, taking only the better alignments for each read
 3) if there are equally good alignments to genome and cassette for one read, take only the cassette ones
 4) categorize the reads based on remaining alignments, into unaligned, cassette, genome-multiple, genome-unique.

The input should be a single fasta or fastq file (file format is determined based on extension), plus an optional preprocessing metadata file.

Output (all filenames are based on OUTFILE_BASENAME provided by command-line option): 
   - metadata file OUTFILE_BASENAME_info.txt - contains exact command used, date/time, system, user, etc, bowtie version, 
                        bowtie command(s) ran, summary bowtie output (normally printed to stdout), final category counts,
                        and all the preprocessing metadata, if a preprocessing metadata file was found.
 AND
   - single output file OUTFILE_BASENAME.sam - alignment file (normally SAM format), containing all filtered alignments
 OR
   - multiple output files:
     - OUTFILE_BASENAME_genomic-unique.sam - unique genomic alignments
     - OUTFILE_BASENAME_cassette.sam - cassette alignments, unique or multiple (multiple cause Warnings printed)
     - OUTFILE_BASENAME_multiple-genomic.* - format (sam or fasta) and how many alignments per read depends on options.
     - OUTFILE_BASENAME_unaligned.* - format (sam or fasta) depends on options.

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011-2012

USAGE: deepseq_alignment_wrapper.py [options] infile -o outfile_basename """

# standard library
from __future__ import division
import sys, os
import unittest
from collections import defaultdict
# other packages
import HTSeq
# my modules
from general_utilities import write_header_data, print_text_from_file, run_command_print_info_output, value_and_percentages
from deepseq_utilities import check_mutation_count_try_all_methods
from basic_seq_utilities import check_fasta_fastq_format, write_fasta_line, get_seq_count_from_collapsed_header
from mutant_analysis_classes import is_cassette_chromosome


def make_aln_dict_from_samfile(samfile, starting_dict=None):
    """ Read SAM-format file, return read:alignment_object_list dict (optionally start with starting_dict); uses HTSeq. """
    new_dict = defaultdict(lambda: [])
    if starting_dict is not None:
        new_dict.update(starting_dict)
    for aln in HTSeq.SAM_Reader(samfile):
        new_dict[aln.read.name] += [aln]
    return dict(new_dict)


def reduce_alignment_dict(readname_to_aln_list):
    """ Given a readname:alignment_list dict, remove all alignments with more than minimum errors per read. """
    for readname,aln_list in readname_to_aln_list.iteritems():
        # nothing to be done for single alignments per read
        if len(aln_list) == 1:
            continue
        # if there are multiple unaligned, just keep one, and go to the next read
        if all([aln.aligned==False for aln in aln_list]):
            readname_to_aln_list[readname] = aln_list[:1]
            continue
        # if there are aligned and unaligned reads, remove the unaligned, and keep going to the next step
        else:
            aln_list = [aln for aln in aln_list if aln.aligned] 
        # if there are aligned reads with different numbers of errors, remove the ones with more than min errors
        min_errors = min([check_mutation_count_try_all_methods(aln) for aln in aln_list])
        aln_list = [aln for aln in aln_list if check_mutation_count_try_all_methods(aln)==min_errors] 
        # finally save the filtered aln_list back to the dictionary
        readname_to_aln_list[readname] = aln_list


def prioritize_cassette_reads(readname_to_aln_list, if_cassette_function=is_cassette_chromosome):
    """ Given a readname:alignment_list dict, whenever there's a cassette alignment, remove all other alignments. """
    for readname,aln_list in readname_to_aln_list.iteritems():
        # nothing to be done for single alignments per read
        if len(aln_list) == 1:
            continue
        # if there are any cassette alignments, keep only the cassette alignments
        if any([if_cassette_function(aln.iv.chrom) for aln in aln_list]):
            aln_list = [aln for aln in aln_list if if_cassette_function(aln.iv.chrom)] 
            readname_to_aln_list[readname] = aln_list


def write_SAM_line_from_HTSeq_aln(htseq_aln, OUTFILE):
    """ Write SAM line from HTSeq alignment (get_sam_line() method); change to tab-separated, add newline to end. """
    OUTFILE.write(htseq_aln.get_sam_line().replace(' ','\t') + '\n')


def categorize_reads_print_to_files(readname_to_aln_list, UNALIGNED_FILE, CASSETTE_FILE, MULTIPLE_GENOMIC_FILE, 
                                    GENOMIC_UNIQUE_FILE, unaligned_as_fasta=True, multiple_to_write=-1, 
                                    input_collapsed_to_unique=False, no_warnings=False):
    """ Decide the proper category for each read, write to appropriate output file; return category counts. 
    
    Categories: unaligned, cassette (one or more cassette alignments - print warning if multiple), 
     genomic-unique (single non-cassette alignment), multiple-genomic (multiple non-cassette alignments. 
    If input_collapsed_to_unique, for the purpose of category counts each read will be counted as N reads, 
     with N determined from readname using the fastx-collapser encoding.
    In the output category counts, cassette-multiple is a special subcategory - anything in it is also counted in cassette.

    Each read is printed to the appropriate outfile (all outfiles should be open file handles); 
     for multiple-genomic, multiple_to_write lines will be written; if unaligned_as_fasta, unaligned reads
     will be written as fasta instead of SAM format (and so will multiple-genomic if multiple_to_write is 0).
    """
    category_readcounts = {'unaligned':0, 'cassette':0, 'multiple-genomic':0, 'genomic-unique':0, 'cassette-multiple':0}

    for readname,aln_list in sorted(readname_to_aln_list.items()):
        readcount = 1 if not input_collapsed_to_unique else get_seq_count_from_collapsed_header(readname)
        # if there's a single alignment, it's unaligned, cassette or genomic-unique
        if len(aln_list) == 1:
            aln = aln_list[0]
            if not aln.aligned:
                category_readcounts['unaligned'] += readcount
                if unaligned_as_fasta:  write_fasta_line(readname, aln.read.seq, UNALIGNED_FILE)
                else:                   write_SAM_line_from_HTSeq_aln(aln, UNALIGNED_FILE)
            elif is_cassette_chromosome(aln.iv.chrom):
                category_readcounts['cassette'] += readcount
                write_SAM_line_from_HTSeq_aln(aln, CASSETTE_FILE)
            else:
                category_readcounts['genomic-unique'] += readcount
                write_SAM_line_from_HTSeq_aln(aln, GENOMIC_UNIQUE_FILE)
        # if there are multiple alignments, it's cassette-multiple (weird!) or multiple-genomic
        else:
            assert all([aln.aligned for aln in aln_list]), "Shouldn't see multiple unaligned lines per read!"
            # multiple-cassette - shouldn't really happen, but write to CASSETTE_FILE
            # MAYBE-TODO come up with something better to do for multiple-cassette cases? If they ever happen.
            if any([is_cassette_chromosome(aln.iv.chrom) for aln in aln_list]):
                assert all([is_cassette_chromosome(aln.iv.chrom) for aln in aln_list]), "Mixed cassette/other!"
                category_readcounts['cassette'] += readcount
                if not no_warnings:
                    print "Warning: multiple cassette alignments! Printing all to cassette file.\n\t%s"%(aln_list)
                category_readcounts['cassette-multiple'] += readcount
                for aln in aln_list:
                    write_SAM_line_from_HTSeq_aln(aln, CASSETTE_FILE)
            # multiple genomic alignments - how many get written depends on multiple_to_write; 
            #  if it's 0, the outfile should be fasta, or else I guess it should be written as unaligned?
            #   (MAYBE-TODO writing single multiple as unaligned not implemented!)
            else:
                category_readcounts['multiple-genomic'] += readcount
                if multiple_to_write == 0:
                    if unaligned_as_fasta:
                        write_fasta_line(readname, aln_list[0].read.seq, MULTIPLE_GENOMIC_FILE)
                    else:
                        raise Exception("Writing 0 multiple alignments in SAM format NOT IMPLEMENTED!")
                else:
                    for aln in aln_list[:multiple_to_write]:
                        write_SAM_line_from_HTSeq_aln(aln, MULTIPLE_GENOMIC_FILE)
    return category_readcounts


def define_option_parser():
    """ Populates and returns an optparse option parser object, with __doc__ as the usage string."""
    from optparse import OptionParser
    parser = OptionParser(__doc__)

    ### test options
    parser.add_option('-t','--test_functionality', action='store_true', default=False, 
                      help="Run the built-in unit test suite (ignores all other options/arguments; default %default).")
    parser.add_option('-T','--test_run', action='store_true', default=False, 
                      help="Run on a test input file, check output against reference files. "
                          + "Ignores all other options/arguments. (default %default).")

    ### functionality options
    parser.add_option('-G','--genome_bowtie_index', metavar='BASENAME', default='Chlre4nm_chl-mit',  
                      help='Bowtie index file set to use for the genome (or everything combined). Default %default.')
    parser.add_option('-C','--cassette_bowtie_index', metavar='BASENAME', default='cassette-pMJ013b',  
                      help="Bowtie index file set to use for the insertion cassette; NONE to only use the index from -G. "
                      +"Default %default.")
    parser.add_option('-e','--allowed_errors', type='int', metavar='N', default=1,  
                      help='How many mismatches to allow in alignments (bowtie -v option); max 3. Default %default.')
    parser.add_option('-m','--multiple_to_show', type='int', metavar='N', default=10,  
                      help='How many multiple alignments to show (bowtie -k/-a option); -1 = no limit. Default %default.')
    parser.add_option('-B','--other_bowtie_options',default='--strata --best --tryhard -S --sam-nosq', metavar='"TEXT"',
                      help="Remaining options to pass to bowtie, as a quoted string (-m and -v will be stripped - "
                      + " specify them separately as -m and -e).  Default \"%default\".")

    ### infile/outfile options
    parser.add_option('-c','--input_collapsed_to_unique', action='store_true', default=False, 
                      help="Use to get correct original total read counts if the data was collapsed to unique sequences "
                          +"using fastx_collapser before alignment (default %default).")
    parser.add_option('-f','--input_metadata_file', metavar='FILE', default='AUTO', 
                      help="Metadata file for the infile. Can be a filename, AUTO to infer from infile name "
                          +"(<infile_basename>_info.txt), or NONE when there is no metadata file. "
                          +"Unless the value is NONE, a warning will be raised if not found. Default: %default")
    parser.add_option('-o','--outfile_basename', metavar='X', default='test', 
                      help="Basename for the outfiles (suffixes will vary). Default %default.")
    parser.add_option('-s','--dont_split_by_category', action='store_true', default=False, 
                      help="Don't split the output file into unaligned, multiple, cassette and unique-genome. "
                      +"Default %default.")
    parser.add_option('-q', '--quiet', action="store_true", default=False,
                      help="Don't print anything to STDOUT, including warnings (default %default).")
    # MAYBE-TODO more stdout verbosity levels?  Do I ever want the full bowtie output, or just the summary?  Summary seems fine...

    return parser

def main(args, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """

    try:
        [infile] = args
    except ValueError:
        parser.print_help()
        sys.exit("Error: exactly one infile required! %s infiles provided: %s"%(len(args), args))
        # MAYBE-TODO bowtie could take multiple infiles, but then I'd have to deal with multiple preprocessing metafiles...

    other_bowtie_options_split = options.other_bowtie_options.split(' ')
    if any([x in other_bowtie_options_split for x in ('-v -e --maqerr -n --seedmms -l --seedlen'.split(' '))]):
        raise Exception("Cannot include -v/-n/-e and related bowtie options in -B!  Use separate -e option for that; "
                        "note that this program allows -v bowtie mode only.")
    if any([x in other_bowtie_options_split for x in ('-m -k -a --all'.split(' '))]):
        raise Exception("Cannot include -m/-a bowtie options in -B!  Use separate -m option for that.")

    specific_bowtie_options = '-v %s'%options.allowed_errors
    if not any([x in options.other_bowtie_options for x in ('-f', '-q')]):
        infile_format = check_fasta_fastq_format(infile)
        if infile_format=='fasta':      specific_bowtie_options += ' -f'
        elif infile_format=='fastq':    specific_bowtie_options += ' -q'
        else:                           raise Exception("Cannot process auto-detected infile format %s!"%infile_format)

    # using a minimum of -k 2 (or -a) in order to make sure I can easily tell multiple from unique alignments
    if options.multiple_to_show == -1:  multiple_bowtie_option = '-a' 
    else:                               multiple_bowtie_option = '-k %s'%max(options.multiple_to_show, 2)

    # output file names: temporary for alignments, final (split or all), metadata info file. 
    outfile_suffix = '.sam' if any([x in options.other_bowtie_options for x in ['-S','--sam']]) else '.map'
    tmpfile_genome = options.outfile_basename + '_tmp_genome' + outfile_suffix
    if options.cassette_bowtie_index != 'NONE':
        tmpfile_cassette = options.outfile_basename + '_tmp_cassette' + outfile_suffix
    if options.dont_split_by_category:
        outfile_all = options.outfile_basename + outfile_suffix
    else:
        outfile_unaligned = options.outfile_basename + '_unaligned.fa'
        outfile_cassette = options.outfile_basename + '_cassette' + outfile_suffix
        outfile_multiple_genomic = options.outfile_basename + '_multiple-genomic'\
                                   + ('.fa' if options.multiple_to_show==0 else outfile_suffix)
        outfile_genomic_unique = options.outfile_basename + '_genomic-unique' + outfile_suffix
    infofile = options.outfile_basename + '_info.txt'

    with open(infofile,'w') as INFOFILE:

        ### write header data
        write_header_data(INFOFILE,options)

        ### run bowtie vs the main/genome index file
        # run 'bowtie --version' to get that data (print to INFOFILE but not stdout)
        INFOFILE.write('\n\n')
        run_command_print_info_output("bowtie --version", INFOFILE, printing_level=0, shell=True)
        # run the actual bowtie alignment command; always print output to stdout as well as INFOFILE
        #   (bowtie actually prints the summary to stderr, not stdout, so I need to print it to stdout in case there's 
        #    an error, so I can see the error message!  Or I could try to detect whether there was an error or not
        #    based on the output contents, but that seems like unnecessary work.)
        INFOFILE.write('\n\n')
        command = "bowtie %s %s %s %s %s %s"%(specific_bowtie_options, multiple_bowtie_option, 
                                      options.other_bowtie_options, options.genome_bowtie_index, infile, tmpfile_genome)
        run_command_print_info_output(command, INFOFILE, printing_level=(not options.quiet), shell=True)

        ### run bowtie vs the cassette index file if given
        if options.cassette_bowtie_index != 'NONE':
            INFOFILE.write('\n\n')
            command = "bowtie %s %s %s %s %s %s"%(specific_bowtie_options, '--all', options.other_bowtie_options, 
                                                  options.cassette_bowtie_index, infile, tmpfile_cassette)
            run_command_print_info_output(command, INFOFILE, printing_level=(not options.quiet), shell=True)

        ### Check that bowtie runs worked
        missing_alnfile_text = "Bowtie run against %s failed! See above or %s file for bowtie error message."
        if not os.access(tmpfile_genome, os.R_OK):
            sys.exit(missing_alnfile_text%(options.genome_bowtie_index, infofile))
        if options.cassette_bowtie_index != 'NONE' and not os.access(tmpfile_cassette, os.R_OK):
            sys.exit(missing_alnfile_text%(options.cassette_bowtie_index, infofile))
        # MAYBE-TODO make sure bowtie errors are printed to stdout even with -1?  Hard - bowtie is unfortunately ANNOYING 
        #  and uses stderr both for normal output and for errors, AND gives no returncode. 

        ### Parse the two alignment files, and merge them together (remove sub-optimal alignments,
        #    (and remove non-cassette ones if there are cassette ones with equal quality); remove alignment files.
        readname_to_aln_list = make_aln_dict_from_samfile(tmpfile_genome)
        if options.cassette_bowtie_index != 'NONE':
            readname_to_aln_list = make_aln_dict_from_samfile(tmpfile_cassette, starting_dict=readname_to_aln_list)
        # MAYBE-TODO right now I'm reading the entire files into memory before merging and processing them, 
        #  which takes a fair amount of memory - could instead write something that would read both alignment files
        #  in parallel and do the merging and output-writing read-by-read.  Do that if I start getting memory issues.
        reduce_alignment_dict(readname_to_aln_list)
        prioritize_cassette_reads(readname_to_aln_list, if_cassette_function=is_cassette_chromosome)
        # delete alignment tmpfiles now that they've been parsed
        os.remove(tmpfile_genome)
        if options.cassette_bowtie_index != 'NONE':
            os.remove(tmpfile_cassette)

        ### Decide the proper category for each read, and write the info to appropriate final output files
        if options.dont_split_by_category:
            with open(outfile_all,'w') as ALL_FILE:
                category_counts = categorize_reads_print_to_files(readname_to_aln_list, ALL_FILE, ALL_FILE, ALL_FILE, 
                                          ALL_FILE, unaligned_as_fasta=False, multiple_to_write=options.multiple_to_show, 
                                          input_collapsed_to_unique=options.input_collapsed_to_unique, 
                                          no_warnings=options.quiet)
        else:
            with open(outfile_unaligned, 'w') as UNALIGNED_FILE:
                with open(outfile_cassette, 'w') as CASSETTE_FILE:
                    with open(outfile_multiple_genomic, 'w') as MULTIPLE_GENOMIC_FILE:
                        with open(outfile_genomic_unique, 'w') as GENOMIC_UNIQUE_FILE:
                            category_counts = categorize_reads_print_to_files(readname_to_aln_list, UNALIGNED_FILE, 
                                                      CASSETTE_FILE, MULTIPLE_GENOMIC_FILE, GENOMIC_UNIQUE_FILE, 
                                                      unaligned_as_fasta=True, multiple_to_write=options.multiple_to_show, 
                                                      input_collapsed_to_unique=options.input_collapsed_to_unique, 
                                                      no_warnings=options.quiet)

        ### print category_readcounts to INFOFILE in a nice way
        text1 = "\n### FINAL ALIGNMENT CATEGORY COUNTS"
        cassette_multiple = category_counts.pop('cassette-multiple')
        total_reads = sum(category_counts.values())
        text2 = "# total reads:  %s"%total_reads
        if options.input_collapsed_to_unique: text2 +=" (uncollapsed readcounts)"
        lines = [text1, text2]
        for category,count in sorted(category_counts.items()):
            text = "# %s:  %s"%(category, value_and_percentages(count, [total_reads]))
            if category=='cassette' and cassette_multiple:  
                text += ' (Warning: %s multiple!!)'%cassette_multiple
            lines.append(text)
        INFOFILE.write('\n')
        for text in lines:
            INFOFILE.write(text + '\n')
            if not options.quiet: print text

        ### copy preprocessing metadata file to the bottom of the new metadata file
        INFOFILE.write("\n\n################## Metadata from input preprocessing ##################\n\n")
        if options.input_metadata_file == 'NONE':
            INFOFILE.write('Not looking for a metadata input file, as specified by options\n')
        else:
            if options.input_metadata_file == 'AUTO':
                # the correct info file for X.txt is X.fa, but for X_5prime.txt it can be either X_5prime.txt or X.txt, so try both.
                #  (in the new preprocessing version all files are X_*prime.txt and the info files are X_info.txt; 
                #   in the old version it was just X.txt and X_info.txt)
                # MAYBE-TODO add a test-case for this thing!  Probably too minor.
                metafile_basename = os.path.splitext(infile)[0] 
                options.input_metadata_file = metafile_basename + '_info.txt'
                if not os.path.exists(options.input_metadata_file):
                    if metafile_basename.endswith('_3prime') or metafile_basename.endswith('_5prime'):
                        options.input_metadata_file = metafile_basename[:-len('_3prime')] + '_info.txt'
                text = 'Automatically determining metadata input file name: %s\n'%options.input_metadata_file
                if not options.quiet:
                    print text,
            else:
                text = 'Metadata input file name provided in options: %s\n'%options.input_metadata_file
            INFOFILE.write(text+'\n')
            if os.path.exists(options.input_metadata_file):
                print_text_from_file(options.input_metadata_file, INFOFILE, printing=False)
            else:
                text = 'Metadata input file %s not found!\n'%options.input_metadata_file
                if not options.quiet:
                    print text,
                INFOFILE.write(text)

        # MAYBE-TODO could actually parse input metadata file and make sure the number of reads found by bowtie matches the one reported at the end of the file, etc...


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    from testing_utilities import run_functional_tests
    test_folder = "test_data"
    infile1 = "%s/INPUT_fasta_for_aln.fa"%test_folder
    # tests in (testname, [test_description,] arg_and_infile_string) format
    test_runs = [ 
        ('aln__basic-1-error',         '-e 1 -m 3     -G Chlre4nm_chl-mit -C cassette-pMJ013b %s -q'%infile1),
        ('aln__collapsed-input',       '-e 1 -m 3 -c  -G Chlre4nm_chl-mit -C cassette-pMJ013b %s -q'%infile1),
        ('aln__show-multiple-0-fasta', '-e 1 -m 0     -G Chlre4nm_chl-mit -C cassette-pMJ013b %s -q'%infile1),
        ('aln__dont-split',            '-e 1 -m 3  -s -G Chlre4nm_chl-mit -C cassette-pMJ013b %s -q'%infile1),
    ]
    # MAYBE-TODO add run-tests for more options/combinations?
    
    # argument_converter converts (parser,options,args) to the correct argument order for main
    argument_converter = lambda parser,options,args: (args, options)
    # use my custom function to run all the tests, auto-detect reference files, compare to output.
    return run_functional_tests(test_runs, define_option_parser(), main, test_folder, 
                                argument_converter=argument_converter, outfile_option='-o') 


class Testing(unittest.TestCase):
    """ Unit-tests this module. """

    def test__(self):
        sys.exit("NO UNIT-TESTS FOR THIS MODULE")
    # LATER-TODO add unit-tests?


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
