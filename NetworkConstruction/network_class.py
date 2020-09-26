import pandas as pd
import numpy as np
import scipy.stats as sy
import os

from NetworkConstruction import utils


class Network:
    def __init__(self, args):
        self.args = args

    def get_matrices(self):
        cis_bp = utils.load_file('input/cisBP.csv', 0, None)

        edge_score_matrix = utils.load_file(self.args.edge_score_path, 0, 0)

        return cis_bp, edge_score_matrix

    def get_labels(self, gene_id, edge_score_matrix):

        """Returns gene and tf labels sorted in alphabetical order

        Optional input tf list and gene list
        """

        if self.args.tf_list_path is not None:
            tfs = utils.load_file(self.args.tf_list_path, None, None)
            tfs = sorted(tfs.values, key=str.lower)
        else:
            tfs = list(set(gene_id) & set(edge_score_matrix.columns.values))
            tfs = sorted(tfs, key=str.lower)


        genes = list(set([i for i in edge_score_matrix.index.values if i.startswith('Y')]) - set(tfs))
        genes = sorted(genes, key=str.lower)

        return tfs, genes

    def construct_network(self):

        cis_bp, edge_score_matrix = self.get_matrices()
        tfs, genes = self.get_labels(cis_bp['GeneID'], edge_score_matrix)

        # Get edge scores of sorted TF and gene labels
        tf_scores = abs(edge_score_matrix[tfs].loc[genes])
        if self.args.rank_low != 0:
            score_low = np.sort(tf_scores.to_numpy().flatten())[-1*self.args.rank_low]
            print('Edge score within rank', self.args.rank_low, ':', score_low)
        else:
            score_low = 0

        # start search
        network_tfs = []
        counter = 0
        max_loops = self.args.max_loops

        print("Starting network construction...")

        while len(network_tfs) < self.args.num_tfs and counter < max_loops:
            rank_high = self.find_tf_rank_threshold(tf_scores)
            score_high = np.sort(tf_scores.to_numpy().flatten())[int(-1*rank_high)]

            if score_high < score_low:
                print("Cannot find", self.args.num_tfs, "TFs with rank threshold of", self.args.rank_low)
                return
            print("Found", self.args.num_tfs, "TFs. Highest rank", int(rank_high), ", Score:", score_high)

            edge_matrix_high = self.binarize_scores(score_high, tf_scores)
            edge_matrix_low = self.binarize_scores(score_low, tf_scores)
            network_tf_index = np.where(utils.sum_col(edge_matrix_high.values) != 0.)[0]
            network_gene_index = np.where(utils.sum_col(np.transpose(edge_matrix_low.values)) != 0.)[0]
            network_tfs = [tf_scores.columns.values[i] for i in network_tf_index]
            network_genes = [tf_scores.index.values[i] for i in network_gene_index]

            print("Removing unidentifiability...")
            tf_scores = self.remove_unidentifiability(network_tfs, edge_matrix_low, tf_scores, self.args.min_targets)
            edge_matrix_high = self.binarize_scores(score_high, tf_scores)
            network_tf_index = np.where(utils.sum_col(edge_matrix_high.values) != 0.)[0]
            network_tfs = [tf_scores.columns.values[i] for i in network_tf_index]
            print(len(network_tfs), "TFs in network")


        labeled_binary_matrix = edge_matrix_low[network_tfs]
        network_genes_index = np.where(utils.sum_col(np.transpose(labeled_binary_matrix.values)) != 0.)[0]
        network_genes = [labeled_binary_matrix.index.values[i] for i in network_genes_index]
        binary_cs = labeled_binary_matrix.loc[network_genes]
        print("Total number of edges:", int(binary_cs.values.sum()))
        print("Final network:", binary_cs.shape[1], "TFs,", binary_cs.shape[0], "genes")
        print("End network construction.")
        self.construct_training_files(binary_cs, self.args.direction)
        self.construct_training_files(binary_cs, self.switch_direction())

    def find_tf_rank_threshold(self, tf_scores):

        """
        :param tf_scores: Labeled edge score matrix
        :return: Highest edge rank threshold that includes at least one edge for each TF
        """

        rank_high = self.args.num_tfs
        cur_num_tfs = self.calculate_cur_num_tfs(rank_high, tf_scores)
        tie = False
        max_rounds = 30
        counter = 0
        steps = self.args.num_tfs
        last_direction = 1
        while cur_num_tfs != self.args.num_tfs and not tie and counter < max_rounds:
            counter += 1
            if self.args.num_tfs > cur_num_tfs:
                if last_direction == -1:
                    if steps == 1:
                        tie = True
                    steps = np.max([np.floor(steps/2), 1])
                    last_direction = -1
                rank_high += steps
            else:
                if last_direction == 1:
                    if steps == 1:
                        tie = True
                    steps = np.max([np.floor(steps/2), 1])
                    last_direction = 1
                rank_high -= steps
            cur_num_tfs = self.calculate_cur_num_tfs(rank_high, tf_scores)

        return rank_high

    def calculate_cur_num_tfs(self, rank_high, tf_scores):
        """
        :param rank_high: Highest edge rank that includes at least one edge for each TF
        :param tf_scores: Labeled edge score matrix
        :return: Current number of TFs included in the selected edges
        """

        score_high = np.sort(tf_scores.to_numpy().flatten())[int(-1*rank_high)]
        binarized_scores = self.binarize_scores(score_high, tf_scores).values
        cur_num_targets = utils.sum_col(binarized_scores)
        cur_num_tfs = len([i for i in cur_num_targets if i != 0.])
        return cur_num_tfs

    def binarize_scores(self, abs_threshold, tf_scores):
        binary_scores = np.zeros(tf_scores.shape)
        for i in np.arange(tf_scores.shape[0]):
            for j in np.arange(tf_scores.shape[1]):
                if tf_scores.values[i][j] > abs_threshold:
                    binary_scores[i][j] = (1*np.sign(tf_scores.values[i][j]))

        binary_scores = self.remove_self_regulation(binary_scores, tf_scores.columns.values, tf_scores.index.values)
        return pd.DataFrame(binary_scores, index=tf_scores.index.values, columns=tf_scores.columns.values)

    @staticmethod
    def remove_self_regulation(binary_scores, tf_labels, gene_labels):
        for i in np.arange(binary_scores.shape[0]):
            for j in np.arange(binary_scores.shape[1]):
                if gene_labels[i] == tf_labels[j]:
                    print("Self regulation removed")
                    binary_scores[i][j] = 0
        return binary_scores

    def remove_unidentifiability(self, network_tfs, labeled_edge_matrix, labeled_score_matrix, min_targets):
        dict_targets = self.compile_target_gene_list(network_tfs, labeled_edge_matrix)
        updated_score_matrix = labeled_score_matrix.copy()

        for tf in dict_targets.keys():
            if len(dict_targets[tf]) < min_targets:
                updated_score_matrix = self.zero_out_genes(updated_score_matrix, dict_targets[tf])

            else:
                dict_targets_reduced = {key: val for key, val in dict_targets.items() if key != tf}
                current_tf_targets = dict_targets[tf]
                shared_genes = \
                    [t for t, targets in dict_targets_reduced.items() if sorted(targets) == sorted(current_tf_targets)]
                if len(shared_genes) > 0:
                    updated_score_matrix = self.zero_out_genes(updated_score_matrix, dict_targets[tf])

        return updated_score_matrix

    @staticmethod
    def compile_target_gene_list(network_tfs, labeled_edge_matrix):
        dict_targets = {tf: [] for tf in network_tfs}
        for tf in network_tfs:
            dict_targets[tf] = list(labeled_edge_matrix[tf].loc[labeled_edge_matrix[tf] != 0.].index)
        return dict_targets

    @staticmethod
    def zero_out_genes(score_matrix, targets):
        for gene in targets:
            score_matrix.loc[gene] = 0
        return score_matrix

    def construct_training_files(self, binary_cs, train_flag):
        """
        :param binary_cs: Binary network, gene x TF
        :param train_flag: Direction of perturbation
        :return: Files, unlabeled
        """
        print("Writing network to files...")
        chosen_tfs = binary_cs.columns.values
        chosen_genes = binary_cs.index.values
        chosen_tfs_common = self.get_common_names(chosen_tfs)
        if train_flag == 0:  # TFKO
            label = 'hk'
            mRNA = pd.read_csv('input/kemmerenKOexpressionMatrixTrimmed.csv', header=0, index_col=0)
        else:
            label = 'zev15'  # Induced
            mRNA = pd.read_csv('input/zev15minExpressionMatrixTrimmed.csv', header=0, index_col=0)

        chosen_tfs_wt = np.append(chosen_tfs.copy(), 'WT')
        sample_labels_nonrelevant = sorted(list(set(mRNA.columns.values) - set(chosen_tfs)), key=str.lower)
        sample_labels_nonrelevant.remove('WT') # Move WT to last column
        sample_labels_nonrelevant.append('WT')
        expression = mRNA.loc[chosen_genes]
        expression_relevant = expression[chosen_tfs_wt]
        expression_nonrelevant = expression[sample_labels_nonrelevant]

        binary_tfa = self.get_binary_tfa({'tfs': chosen_tfs}, mRNA.columns.values, train_flag)
        binary_tfa_relevant = binary_tfa[chosen_tfs_wt]
        binary_tfa_nonrelevant = binary_tfa[sample_labels_nonrelevant]

        binary_cs_signed = binary_cs.multiply(np.sign(expression_relevant[chosen_tfs].values))

        if train_flag == 0:
            binary_cs_signed = binary_cs_signed.multiply(-1)

        binary_cs_corr_signed, cs_corr = self.get_corr_binary_cs(binary_cs, chosen_tfs, chosen_genes, mRNA[sample_labels_nonrelevant])

        cur_dir = os.getcwd()
        if not os.path.exists(cur_dir + '/' + str(label)):
            os.makedirs(cur_dir + '/' + str(label))


        pd.DataFrame(chosen_tfs).to_csv(str(label)+'/'+str(label)+'_tf_list.csv', header=False, index=False)
        pd.DataFrame(chosen_genes).to_csv(str(label)+'/'+str(label)+'_gene_list.csv', header=False, index=False)
        pd.DataFrame(mRNA.columns.values).to_csv(str(label)+'/'+str(label)+'_sample_list_all.csv',
                                                 header=False, index=False)
        pd.DataFrame(chosen_tfs_wt).to_csv(str(label)+'/'+str(label)+'_sample_list_perturbed.csv',
                                           header=False, index=False)
        pd.DataFrame(sample_labels_nonrelevant).to_csv(str(label) + '/' + str(label) + '_sample_list_unperturbed.csv',
                                                       header=False, index=False)

        expression.to_csv(str(label)+'/'+str(label)+'_expression_all.csv', header=False, index=False)
        expression_relevant.to_csv(str(label)+'/'+str(label)+'_expression_perturbed.csv', header=False, index=False)
        expression_nonrelevant.to_csv(str(label) + '/' + str(label) + '_expression_unperturbed.csv',
                                      header=False, index=False)

        binary_tfa.to_csv(str(label) + '/' + str(label) + '_binary_tfa_all.csv', header=False, index=False)
        binary_tfa_relevant.to_csv(str(label) + '/' + str(label) + '_binary_tfa_perturbed.csv',
                                   header=False, index=False)
        binary_tfa_nonrelevant.to_csv(str(label) + '/' + str(label) + '_binary_tfa_unperturbed.csv',
                                      header=False, index=False)

        binary_cs_signed.to_csv(str(label) + '/' + str(label) + '_binary_cs_signed.csv', header=False, index=False)
        binary_cs_corr_signed.to_csv(str(label) + '/' + str(label) + '_binary_cs_corr_signed.csv',
                                     header=False, index=False)
        print("Done!")

    def switch_direction(self): #TODO: make constant vars
        if self.args.direction == 0:
            return 2
        if self.args.direction == 2:
            return 0

    @staticmethod
    def get_binary_tfa(tf_labels, sample_labels, perturbation):
        binary_tfa = np.ones((len(tf_labels['tfs']), len(sample_labels)))
        for i in np.arange(len(tf_labels['tfs'])):
            for j in np.arange(len(sample_labels)):
                if sample_labels[j] in tf_labels['tfs'][i]:
                    binary_tfa[i][j] = perturbation
        return pd.DataFrame(binary_tfa, index=tf_labels['tfs'], columns=sample_labels)

    @staticmethod
    def get_corr_binary_cs(cs, chosen_tfs, chosen_genes, expression):
        binary_cs_signed = cs.copy()
        cs_corr = cs.copy()
        for g in chosen_genes:
            gene_expression = expression.loc[g]
            for t in chosen_tfs:
                tf_expression = expression.loc[t]
                if cs[t].loc[g] != 0:
                    r_val, p_val = sy.pearsonr(gene_expression, tf_expression)
                    binary_cs_signed[t].loc[g] = np.sign(r_val)
                    cs_corr[t].loc[g] = p_val
        return binary_cs_signed, cs_corr

    @staticmethod
    def get_common_names(tfs):
        common_names = pd.read_csv('input/YeastGeneNames.tsv', delimiter='\t', header=None, index_col=0)
        common_names_chosen = [common_names.loc[i].values[0] for i in tfs]
        return common_names_chosen
