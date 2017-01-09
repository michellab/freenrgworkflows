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
             '-o',
             '--network_output',
             help='File to write final free energies to based on network analysis',
             metavar='FILE'
    )
    # parser.add_argument(
    #         '-e',
    #         '--experiments',
    #         help='File containing experimental data results',
    #         metavar='FILE'
    # )
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
    print ("\n\n##########################NETWORKANALYSIS WITH NETWORKX ######################################")
    print ("#\n### PARAMETERS\n#")
    print (args.files)


    #Do the network analysis
    pG = PerturbationGraph()
    pG.populate_pert_graph(args.files[0])
    if len(args.files) > 1:
        for f in args.files[1:]:
            pG.add_data_to_graph(f)
    pG.compute_weighted_avg_paths(args.target_compound)
    comp_DDG = pG.freeEnergyInKcal

    if args.save_data:
        pG.save_average_paths(args.network_output, comp_DDG)

    #Read experimental data



    ############################################################################
    #
    #   say good bye
    #
    ############################################################################
    print ("#\n###################That's it, now it is time to put the kettle on ##############################\n#")
    print ("#                  Thank you for using the network analysis package!\n#")
    print ("#\n################################################################################################\n\n")

