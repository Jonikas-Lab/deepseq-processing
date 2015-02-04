#! /usr/bin/env python2.7

"""
Take LEAP-Seq cassette-side read files (fastq), which should be <IB><cassette><genomic-flanking-seq>.  Look for the full or truncated cassette sequence after the expected IB length; output the IB and flanking seqs to separate fasta output files; output truncated cassette lengths to another output file, and a summary to stdout.
Output files will be <outfile_base>_IB.fa, <outfile_base>_flank.fa, <outfile_base>_unstripped.fa and <outfile_base>_cassette-info.txt (tab-separated).
 -- Weronika Patena, 2015
USAGE: deepseq_strip_cassette [options] infile outfile_base
"""

# standard library
from __future__ import division
import sys
import unittest
import collections
# other packages
# my modules
import basic_seq_utilities

class CassetteStrippingError(Exception): pass


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
    parser.add_option('-C', '--cassette_sequence', default=None, metavar='N',
                      help="Expected assette sequence (default %default)")
    parser.add_option('-c', '--min_cassette_length', type='int', default=5, metavar='N', 
                      help="Minimum cassette length to expect (since we're allowing truncations) (default %default).")
    parser.add_option('-l', '--allowed_IB_lengths', default='21,22,23', metavar='N,M', 
                      help="How many bases before the cassette to expect (comma-separated numbers) (default %default).")

    ### cosmetic options
    parser.add_option('-q','--quiet', action='store_true', default=False, 
                      help="Don't print anything to stdout (default %default).")
    return parser


def strip_cassette_seq(input_seq, cassette_seq, IB_lens, min_len=5):
    """ Split input_seq into sections before and after cassette_seq, allowing cassette truncations.

    The input_seq is expected to look like this: xxxxxxxxxCCCCCCCCCCCCCCCyyyyyyyy, 
    where CCCCCCCCCC is either all of cassette_seq or a prefix of it (of at least min_len), 
    and the length of xxxxxxxx is one of IB_lens. 

    For each possible length in IB_lens, find the longest possible prefix of the 
    cassette sequence at that position in input_seq, allowing no errors. 

    If the cassette seq is found for multiple lengths, raise Exception.
    If it's not found at all, return None.
    If it's found once, return a (seq_before_cassette, seq_after_cassette, cassette_length) tuple.
    """
    cassette_seq = cassette_seq.lower()
    input_seq = input_seq.lower()
    results = []
    for before_length in IB_lens:
        if_matched_bases = list(c==i for (c,i) in zip(cassette_seq, input_seq[before_length:]))
        # this would probably be faster if I encoded these as large integers and used bitwise operations
        try:                cassette_len = if_matched_bases.index(False)
        except ValueError:  cassette_len = len(if_matched_bases)
        if cassette_len >= min_len:
            results.append((before_length, cassette_len))
    if len(results) > 1:
        raise CassetteStrippingError("Cassette %s found in seq %s in multiple (pos,length) pairs! %s"%(
                                        cassette_seq, input_seq, results))
    if not len(results):
        return
    before_length, cassette_len = results[0]
    return input_seq[:before_length], input_seq[before_length+cassette_len:], cassette_len
    # TODO implement allowing mismatches - slightly complicated, need to define a trade-off between #mismatches and length
    # TODO implement allowing indels - more complicated, requires full alignment really - use Biopython?  Or even do an actual bowtie2 local alignment of each sequence to the expected cassette or something?


def strip_cassette_from_file(infile, outfile_base, cassette_sequence, allowed_IB_lengths, 
                             min_cassette_len=5, quiet=False):
    """ Strip cassette seq from middles of LEAP-Seq reads, outputting starts/ends to outfiles. 
    
    Keep track of stripped cassette lengths, and print summary to stdout.
    """
    if cassette_sequence is None:
        raise CassetteStrippingError("Must provide cassette sequence!")
    try:
        if len(cassette_sequence) < min_cassette_len:
            raise CassetteStrippingError("Cassette sequence must be longer than min_cassette_len! (%s, %s)"%(
                                                                            min_cassette_len, cassette_sequence))
    except TypeError:
        raise CassetteStrippingError("Cassette sequence must be a string! (%s)"%cassette_sequence)
    if type(allowed_IB_lengths) == str:
        allowed_IB_lengths = [int(x) for x in allowed_IB_lengths.split(',')]
    total_seqs = 0
    results_IB_lengths = collections.Counter()
    results_cassette_lengths = collections.Counter()
    resuts_N_unstripped = 0
    with open(outfile_base+'_cassette-info.txt', 'w') as OUTFILE_cassette:
        OUTFILE_cassette.write("seq_header\tcassette_length\n")
        with open(outfile_base+'_IB.fa', 'w') as OUTFILE_IB:
          with open(outfile_base+'_flank.fa', 'w') as OUTFILE_flank:
            with open(outfile_base+'_unstripped.fa', 'w') as OUTFILE_unstripped:
              for (name, seq) in basic_seq_utilities.name_seq_generator_from_fasta_fastq(infile):
                  total_seqs += 1
                  output = strip_cassette_seq(seq, cassette_sequence, allowed_IB_lengths, min_cassette_len)
                  if output is None:
                      resuts_N_unstripped += 1
                      OUTFILE_unstripped.write('>%s\n%s\n'%(name, seq))
                  else:
                      IB, flank, cassette_len = output
                      results_IB_lengths[len(IB)] += 1
                      results_cassette_lengths[cassette_len] += 1
                      OUTFILE_IB.write('>%s\n%s\n'%(name, IB))
                      OUTFILE_flank.write('>%s\n%s\n'%(name, flank))
                      OUTFILE_cassette.write('%s\t%s\n'%(name, cassette_len))
    if not quiet:
        print "Total sequences:       %s"%total_seqs
        print "Unstripped sequences:  %s (there may also be partially-stripped ones!)"%general_utilities.value_and_percentages(
                                                                                        resuts_N_unstripped, [total_seqs])
        print "IB lengths:"
        print '\n'.join("   %-7s: %s"%x for x in sorted(results_IB_lengths.items()))
        print "Stripped cassette lengths:"
        print '\n'.join("   %-7s: %s"%x for x in sorted(results_cassette_lengths.items()))


def main(args, options):
    """ Run the main functionality of the module (see module docstring for more information), excluding testing.
    The options argument should be generated by an optparse parser.
    """
    try:
        infile, outfile_base = args
    except ValueError:
        parser.print_help()
        sys.exit("\nError: exactly one infile and one outfile basename required!")
    strip_cassette_from_file(infile, outfile_base, options.cassette_sequence, 
                             options.allowed_IB_lengths, options.min_cassette_length, options.quiet)


def do_test_run():
    """ Test run: run script on test infile, compare output to reference file."""
    from testing_utilities import run_functional_tests
    test_folder = "test_data"
    infile1 = test_folder + '/INPUT_strip_cassette.fa'
    # tests in (testname, [test_description,] arg_and_infile_string) format
    test_runs = [
        ('strip__basic', 'basic cassette-stripping', '-C tgtgtgtg -c 3 -l 4,5 -q %s'%infile1), 
                ]
    # argument_converter converts (parser,options,args) to the correct argument order for main
    argument_converter = lambda parser,options,args: (args, options)
    # use my custom function to run all the tests, auto-detect reference files, compare to output.
    return run_functional_tests(test_runs, define_option_parser(), main, test_folder, 
                                argument_converter=argument_converter) 
    # LATER-TODO add run-tests!


class Testing(unittest.TestCase):
    """ Runs unit-tests for this module. """

    def test__strip_cassette_seq(self):
        # cassette found
        for cassette in 'CC CCCCCCC cccccccc ccaaggggg CCATG'.split():
            for lens in ([2], [1,2,3], [2,10]):
                self.assertEqual(strip_cassette_seq('aacctt', cassette, lens, 2), ('aa', 'tt', 2))
        self.assertEqual(strip_cassette_seq('aacctt', 'cct', [1,2,3], 2), ('aa', 't', 3))
        # cassette not found
        for cassette in 'ctct aaa gactttt'.split():
            self.assertIsNone(strip_cassette_seq('aacctt', cassette, [2], 2))
        # cassette found twice - error
        self.assertRaises(CassetteStrippingError, strip_cassette_seq, 'aatttt', 'ttt', [2,3], 2)


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
