import RNA
from nupack import *
config.parallelism = True


class EnergyCalculator(object):
    def __init__(self):
        pass

    @staticmethod
    def calc_self_mfe_structure_and_fold_energy(seq):
        ss, mfe = RNA.fold(seq)
        return ss, mfe

    @staticmethod
    def calc_duplex_mfe_structure_and_fold_energy(seq1, seq2):
        # method 1 duplexfold
        # import RNA
        # energies = []
        # for seq in seqs:
        #     duplex = RNA.duplexfold(ref_seq, seq)
        #     energies.append(duplex.energy)

        # method 2 nupack
        # import nupack
        # energies = []
        # nup_model = nupack.Model(material='RNA')
        # for seq in seqs:
        #     mfe_structures = nupack.mfe(strands=[ref_seq, seq], model=nup_model)
        #     if len(mfe_structures) < 1:
        #         mfe = 0.0
        #     else:
        #         mfe = mfe_structures[0].energy
        #     energies.append(mfe)

        # method 3 nupack
        import nupack
        s1 = nupack.Strand(seq1, name='s1')
        s2 = nupack.Strand(seq2, name='s2')
        c = nupack.Complex([s1, s2], name='c')
        # c1 = nupack.Complex([s1], name='c1')
        # c2 = nupack.Complex([s2], name='c2')
        # c1_c2 = nupack.Complex([c1, c2], name='c1c2')
        # c2 = nupack.Complex([s, s], name='c2')

        t = nupack.Tube({s1: 1e-8, s2: 1e-8},
                        complexes=nupack.SetSpec(max_size=2, include=[c]),
                        name='t')

        # t = nupack.Tube({s: 1e-8},
        #                 name='t',
        #                 complexes=nupack.SetSpec(include=[c1, c2]))

        model = nupack.Model()
        anal_result = nupack.tube_analysis([t], model=model,
                                           compute=['pfunc', 'pairs', 'mfe', 'sample', 'subopt'],
                                           options={'num_sample': 2, 'energy_gap': 0.5})

        # print(anal_result)

        # result = anal_result['(s2+s1)']
        # result = anal_result['(s1+s2)']
        result = anal_result[c]

        # con_str = str(s2) + '+' + str(s1)
        concentrations = anal_result.tubes[t].complex_concentrations[c]
        return result, concentrations

        # method 4 fold_compound
        # fc = RNA.fold_compound(seq1 + "&" + seq2)
        # ss, mfe = fc.mfe_dimer()
        # return ss, mfe

    @staticmethod
    def calc_self_mfe_structure_and_fold_energy_and_probabilities(seq):

        import nupack

        s = nupack.Strand(seq, name='s')
        c1 = nupack.Complex([s], name='c1')
        # c2 = nupack.Complex([s, s], name='c2')

        t = nupack.Tube({s: 1e-8},
                        # complexes=nupack.SetSpec(include=[c1]),
                        name='t')

        # t = nupack.Tube({s: 1e-8},
        #                 name='t',
        #                 complexes=nupack.SetSpec(include=[c1, c2]))

        model = nupack.Model()
        anal_result = nupack.tube_analysis([t], model=model,
                                           compute=['pfunc', 'pairs', 'mfe', 'sample', 'subopt'],
                                           options={'num_sample': 2, 'energy_gap': 0.5})

        # print(anal_result)

        result = anal_result['(s)']
        # print('Physical quantities for complex c')
        # print('Complex free energy: %.2f kcal/mol' % result.free_energy)
        # print('Partition function: %.2e' % result.pfunc)
        # print('MFE proxy structure: %s' % result.mfe[0].structure)
        # print('Free energy of MFE proxy structure: %.2f kcal/mol' % result.mfe[0].energy)
        # print('Equilibrium pair probabilities: \n%s' % result.pairs)

        # aa_result = my_result['(s+s)']

        # result = anal_result[s]
        # print('MFE proxy structure:\n%s' % result.mfe[0].structure.matrix())

        return result

        # import matplotlib.pyplot as plt
        #
        # plt.imshow(my_result[c].pairs.to_array())
        # plt.xlabel('Base index')
        # plt.ylabel('Base index')
        # plt.title('Pair probabilities for complex c')
        # plt.colorbar()
        # plt.clim(0, 1)
        # # plt.savefig('my-figure.pdf')  # optionally, save a PDF of your figure
