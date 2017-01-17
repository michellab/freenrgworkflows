#!/usr/bin/env python
# -*- coding: utf-8 -*-

#!/usr/bin/env python

# This file is part of freenrgworkflows.
#
# Copyright 2016,2017 Julien Michel Lab, University of Edinburgh (UK)
#
# freenrgworkflows is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__author__ = "Antonia Mey"
__email__ = "antonia.mey@ed.ac.uk"


####################################################################################################
#
#   IMPORTS
#
####################################################################################################

from networkanalysis.networkanalysis import *
from networkanalysis.experiments import *
from networkanalysis.stats import *
import networkanalysis
from argparse import ArgumentParser, FileType
import numpy as np


####################################################################################################
#
#   MAIN PART
#
####################################################################################################

if '__main__' == __name__:

    ############################################################################
    #
    #   capture the command line arguments
    #
    ############################################################################

    parser = ArgumentParser()
    parser.add_argument(
        'files',
        help='networkx compatible csv file/files of the computed free energies with Sire',
        nargs='*',
        metavar='FILE'
    )
    parser.add_argument(
             "--target_compound",
             help="Name of the reference compound with respect to which the free energy should be computed",
             metavar='STRING'
    )
    parser.add_argument(
             "--intermed_ID",
             help="Name of the reference compound with respect to which the free energy should be computed",
             metavar='STRING'
    )
    parser.add_argument(
             '-o',
             '--network_output',
             help='File to write final free energies to based on network analysis',
             metavar='FILE'
    )
    parser.add_argument(
             '-e',
             '--experiments',
             help='File containing experimental data results',
             default = None,
             metavar='FILE'
    )
    parser.add_argument(
            "--stats",
            help="Print correclation statistics between computated and experimental data",
            action='store_true'
    )
    # parser.add_argument(
    #         "--maxiter",
    #         help="limit the number of fixed point iterations",
    #         type=int,
    #         default=100,
    #         metavar='INT'
    # )
    parser.add_argument(
            "--save_data",
            help="Saves network data output",
            action='store_true'
    )

    args = parser.parse_args()

    ############################################################################
    #
    #   check mandatory command line arguments
    #
    ############################################################################
    if 1 > len(args.files):
        print ("ERROR: you must give at least one networkanalysis networkx compatible file")
        exit(1)


    ############################################################################
    #
    #   write header
    #
    ############################################################################
    print ("\n\n################# NETWORKANALYSIS v. %s WITH NETWORKX ################################" %networkanalysis.__version__)
    print ("\n\n########################## Parameters ######################################")
    print ("filelist: \t\t\t\t%s" %args.files)
    print ("target compound: \t\t\t%s" %args.target_compound)
    print ("intermed_ID: \t\t\t\t%s" %args.intermed_ID)
    print ("Network computed free energies file: \t%s" %args.network_output)
    print ("IC50s datafile: \t\t\t%s" %args.experiments)
    print ("Correlation statistics:\t\t\t%s" %args.stats)
    print ("#############################################################################\n\n")


    #Do the network analysis
    pG = PerturbationGraph()
    pG.populate_pert_graph(args.files[0])
    if len(args.files) > 1:
        for f in args.files[1:]:
            pG.add_data_to_graph(f)
    pG.compute_weighted_avg_paths(args.target_compound)
    pG.format_free_energies(merge_BM=True, intermed_ID=args.intermed_ID, weighted = True)
    comp_DDG = pG.freeEnergyInKcal

    if args.save_data:
        pG.write_free_energies(comp_DDG, filename=args.network_output)

    #Read experimental data
    if args.experiments != None:
        ex = ExperimentalData()
        ex.compute_DDG_from_IC50s(args.experiments,reference=args.target_compound)
        exp_DDG = ex.freeEnergiesInKcal
        if args.stats:
            stats = freeEnergyStats()
            stats.generate_statistics(comp_DDG,exp_DDG,repeats=1000)

            print ("\n\n########################## Statistics ######################################")
            print (" R and error = %f ± %f" %(stats.R, stats.R_error))
            print (" R2 and error = %f ± %f" %(stats.R2, stats.R2_error))
            print (" tau and error = %f ± %f" %(stats.tau, stats.tau_error))
            print (" MUE and error = %f ± %f" %(stats.mue, stats.mue_error))
            print ("#############################################################################\n\n")

    #create plots 


    ############################################################################
    #
    #   say good bye
    #
    ############################################################################
    print ("\n#################################################################################################\n#")
    print ("#                 That's it, now it's time to put the kettle on ")
    print ("#                Thank you for using the network analysis package!")
    print ("#\n################################################################################################\n\n")

