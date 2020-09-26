import argparse
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from NetworkConstruction import network_class



parser = argparse.ArgumentParser(description='Creates a TFs-target genes network.')
parser.add_argument('-n',
                    dest='num_tfs',
                    help='Number of TFs to include in the network.',
                    required=True)
parser.add_argument('-d',
                    dest='direction',
                    help='Direction of perturbation. 0 for Knock out; 2 for Induced.',
                    required=True)
parser.add_argument('-e',
                    dest='edge_score_path',
                    help='Labeled edge score matrix (gene x TF).',
                    required=True)
parser.add_argument('-lr',
                    dest='rank_low',
                    help='Lowest edge rank before search stops. 0 for no limit',
                    required=False,
                    default='1250')
parser.add_argument('-mt',
                    dest='min_targets',
                    help='Minimum number of targets for a TF in the network.',
                    required=False,
                    default='2')
parser.add_argument('-tf',
                    dest='tf_list_path',
                    help='Path to list of TFs of interest. ',
                    required=False,
                    default=None)
parser.add_argument('-l',
                    dest='max_loops',
                    help='Maximum number of rounds to search for.',
                    required=False,
                    default='5')


args = parser.parse_args()

args.num_tfs = int(args.num_tfs)
args.direction = int(args.direction)
args.rank_low = int(args.rank_low)
args.min_targets = int(args.min_targets)
args.max_loops = int(args.max_loops)

network_class.Network(args).construct_network()




