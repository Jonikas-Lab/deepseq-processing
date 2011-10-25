#!/usr/bin/env python
""" This is mostly a wrapper around fastq/fasta preprocessing utilities (fastx_trimmer, cutadapt, fastx_collapser), that runs the programs on your input file and keeps track of metadata in a hopefully sensible manner.

I made separate options for commonly used options of the wrapped programs; if you want to do anything not covered by those, use the -T, -A, -C options to provide full option lists to each wrapped program - the option lists will be passed to the program with no parsing or changes.

Output - there are actually two outfiles generated (with names based on the -o option - calling the value X here):
    - sequence file - <X>.fq or .fa, depending on input format etc - contains the sequences
    - metadata file - <X>_info.txt - contains general information: 
        exact command used, date/time, system, user, 
        summary output from cutadapt (normally printed to stdout), 
        a count of discarded sequences from cutadapt and the other programs

 -- Weronika Patena, Jonikas Lab, Carnegie Institution, 2011

USAGE: deepseq_preprocessing.py [options] infile """

import sys, os, subprocess
from general_utilities import write_header_data
from read_length_distribution import read_length_distribution


def check_readlengths(description, infile, OUTFILE, quiet):
    read_length_data = read_length_distribution([infile],include_zeros=False,verbose=1,OUTPUT=None)
    output = "# %s (%s) data:\n%s"%(description, infile, ''.join(read_length_data))
    OUTFILE.write(output+'\n')
    if not quiet:   print output


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
    # MAYBE-TODO add options for commonly-used fastx_trimmer and cutadapt options (fastx_collapser doesn't have any)
    # MAYBE-TODO make sure that for each wrapped program, either specific options OR the general option (-T/-A) is provided and used, NOT BOTH, and print the information for the user!
    parser.add_option('-k', '--keep_tmpfiles', action="store_true", default=False, 
                      help="Don't delete the temporary/intermediate files (all will have the same prefix as the outfile).")
    parser.add_option('-q', '--quiet', action="store_true", default=False, help="Don't print anything to STDOUT.")
    try:
        (options, [infile]) = parser.parse_args()
    except ValueError:
        parser.print_help()
        sys.exit("Error: exactly one infile required!")
    # MAYBE-TODO implement option with multiple infiles? Need to make sure they're the same fa/fq type etc...

    # outfile and tmpfile names
    infile_suffix = os.path.splitext(infile)[1]
    outfile_suffix = '.fa' if options.collapse_to_unique else infile_suffix
    outfile = options.outfile_basename + outfile_suffix
    infofile = options.outfile_basename + '_info.txt'
    tmpfile1 = tmpfile1_original = options.outfile_basename + '_tmpfile1' + infile_suffix
    tmpfile2 = tmpfile2_original = options.outfile_basename + '_tmpfile2' + infile_suffix
    tmp_cutadapt_info = options.outfile_basename + '_tmp_cutadapt_info.txt'
    tmp_cutadapt_untrimmed = options.outfile_basename + '_tmp_cutadapt_untrimmed' + infile_suffix
    tmp_cutadapt_tooshort = options.outfile_basename + '_tmp_cutadapt_too-short' + infile_suffix
    tmp_collapser_info = options.outfile_basename + '_tmp_collapser_info.txt'
    
    with open(infofile,'w') as INFOFILE:

        # write header data
        write_header_data(INFOFILE,options)
        INFOFILE.write('\n')

        # 1. look at the infile; make sure it's readable, etc
        #       (read_length_distribution uses HTSeq and autodetects fa/fq format)
        check_readlengths("Original infile", infile, INFOFILE, options.quiet)

        # 2. run fastx_trimmer
        if options.full_fastx_trimmer_options == 'NONE':
            if not options.quiet:   print "Not running fastx_trimmer, since NONE was passed to -T."
            tmpfile1 = infile
        else:
            command = "fastx_trimmer %s -i %s -o %s"%(options.full_fastx_trimmer_options, infile, tmpfile1)
            output = "# Running fastx_trimmer: %s"%command
            INFOFILE.write(output+'\n')
            if not options.quiet:   print output
            subprocess.call([command], shell=True)
            check_readlengths("fastx_trimmer output", tmpfile1, INFOFILE, options.quiet)
            # TODO do something sensible if no reads are left!

        # 3. run cutadapt
        if options.full_cutadapt_options == 'NONE':
            if not options.quiet:   print "Not running cutadapt, since NONE was passed to -A."
            tmpfile2 = tmpfile1
        else:
            command = "cutadapt %s -o %s --too-short %s --untrimmed %s %s > %s"%(options.full_cutadapt_options, tmpfile2,\
                                     tmp_cutadapt_tooshort, tmp_cutadapt_untrimmed, tmpfile1, tmp_cutadapt_info)
            # TODO grab the text from tmp_cutadapt_info and put it in INFOFILE!
            # TODO count the sequences in tmp_cutadapt_tooshort and tmp_cutadapt_untrimmed and put info in INFOFILE!
            output = "# Running cutadapt: %s"%command
            INFOFILE.write(output+'\n')
            if not options.quiet:   print output
            subprocess.call([command], shell=True)
            check_readlengths("cutadapt output", tmpfile2, INFOFILE, options.quiet)
            # TODO do something sensible if no reads are left!

        # 4. run fastx_collapser
        if not options.collapse_to_unique:
            if not options.quiet:   print "Not running fastx_collapser, since -C option was not used."
            os.rename(tmpfile2,outfile)
            # Note for fastx_collapser, but also for the others - NONE is necessary here, can't just use '', because 
            #    fastx_collapser works fine with no options, so '' is a sensible input and can't be used to turn it off.
        else:
            command = "fastx_collapser -v -i %s -o %s 2> %s"%(tmpfile2, outfile, tmp_collapser_info)
            # TODO grab the text from tmp_collapser_outfile and put it in INFOFILE!
            output = "# Running fastx_collapser: %s"%command
            INFOFILE.write(output+'\n')
            if not options.quiet:   print output
            subprocess.call([command], shell=True)
            check_readlengths("fastx_collapser output", outfile, INFOFILE, options.quiet)
            # TODO do something sensible if no reads are left!

        # MAYBE-TODO the read_length_distribution step etc should probably be optional, 
        #   since it's going to take extra time/cpu and the files may be large...  Make an option for that, implement.
        # To do things quickly, just use a single command with pipes, instead of multiple commands with intermediate files.
        # MAYBE-TODO also maybe have an intermediate level where the output doesn't include 
        #   the whole read-length distribution but just the total read count?

    # remove tmpfiles
    # MAYBE-TODO make this optional?
    # need to use the tmpfile*_original names here because I do "tmpfile1 = infile" etc if skipping steps, 
    #   and I don't want to remove the infile!
    if not options.keep_tmpfiles:
        for tmpfile in [tmpfile1_original, tmpfile2_original, tmp_cutadapt_info, 
                        tmp_cutadapt_untrimmed, tmp_cutadapt_tooshort, tmp_collapser_info]:
            if os.path.exists(tmpfile):
                os.remove(tmpfile)



